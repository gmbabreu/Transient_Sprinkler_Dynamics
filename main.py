"""
    Written by Rachel Bertaud during PhD at Colorado School of Mines, 2026

    To analyze the transient behavior of the sprinkler, it is imperative to know
    the natural frequency (omega) and the damping coefficient (gamma) of the system.
    By assuming that the tail of the measured angular position data of the sprinkler
    fits the ODE for a damped, harmonic oscillator, we have found a method to find
    those values, and use them to extract the torque signal of the system.

    This code:
        0. Reads shared configuration and selects a Reynolds number/direction group
        1. Reads angular position CSVs and requests or reuses each trial's fit time
        2. Estimates omega (natural frequency) and gamma (damping coefficient)
        3. Keeps initial measurement noise and pads the recording for the FFT
        4. Uses Fourier transforms to extract the underlying torque signal
        5. Does a forward ODE solve using the extracted torque to verify the results
        6. Computes separate signed torque integrals around startup and shutdown
        7. Saves individual torque/impulse CSVs and plots the group's results

    All sections of this code are clearly labelled as above for
    process transparency. Sections 2-6 are narrated in scientific_pipeline;
    analyze_trial runs them around each trial's saved/interactive selections;
    section 7 saves and compares successful trials in run_group.

    Run, for example: python main.py 1000 f
    Named options also work: python main.py --re 1000 --direction f
    All calculations and saved selections use original video seconds;
    only comparison plots shift the time axis by the effective startup time.
"""

# SECTION ZERO - USER DEFINED INPUTS
###################################################################################################
# Shared inputs now live in experiment_config.json so MATLAB and Python use
# the same groups, physical constants and analysis settings. Command-line
# arguments select the group; per-trial metadata retains interactive choices.

# DEFINE DEPENDENCIES AND UDFs
import argparse
import csv
import json
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np

# DATA READ FUNCS
from Main_Functions.dataread_funcs import read_data, plot_data

# PROCESS FUNCS
from Main_Functions.process_funcs import validate_recording

# PLOT FUNCS
from Main_Functions.plot_funcs import plot_torque, plot_torque_int, plot_phi_gen

# IMPULSE FUNCS
from Main_Functions.impulse_funcs import (
    report_impulses, save_impulse_summary, relative_error_metrics,
    mark_impulse_windows, plot_impulse_diagnostics, select_impulse_windows, validate_window,
)

# SCIENTIFIC STAGES (pure calculations; workflow remains here)
from Main_Functions.scientific_pipeline import reconstruct_torque, measure_event_impulses


def parse_trial(trial_id):
    """Read Reynolds number and flow direction from an ID such as '1000f2'.

    The final number identifies the repeat; it stays in the full trial ID so
    files for repeats remain distinct, but is not part of group selection.
    """
    match = re.fullmatch(r"([1-9]\d*)([fr])([1-9]\d*)", trial_id)
    if match is None:
        raise ValueError(f"Invalid trial ID: {trial_id!r}; expected e.g. 1000f1")
    return int(match[1]), match[2]


def load_config(path=Path(__file__).with_name("experiment_config.json")):
    """Read the same JSON as MATLAB's loadBatchConfig and resolve its paths.

    Relative paths are based on the configuration file, not the shell's current
    directory. Expand each group into repeats 1..count, in group order. A scalar
    trials_per_group supplies every count; an array supplies one per group.
    The rest of the workflow still receives full trial IDs.
    Reject invalid group/count settings before processing any outputs.
    """
    path = Path(path).resolve()
    config = json.loads(path.read_text())
    config["paths"] = {key: (path.parent / value).resolve()
                       for key, value in config["paths"].items()}
    groups = config["trial_groups"]
    if (not isinstance(groups, list) or not groups
            or any(not isinstance(group, str) or re.fullmatch(r"[1-9]\d*[fr]", group) is None
                   for group in groups)):
        raise ValueError("trial_groups must be a nonempty list of group IDs such as 1000f.")
    if len(groups) != len(set(groups)):
        raise ValueError("Configured trial groups must be unique.")
    counts = config["trials_per_group"]
    if not isinstance(counts, list):
        counts = [counts] * len(groups)
    if len(counts) != len(groups):
        raise ValueError("trials_per_group must be a positive integer or an array "
                         "with one positive integer per trial group.")
    for count in counts:
        if (isinstance(count, bool) or not isinstance(count, (int, float))
                or not np.isfinite(count) or count < 1 or count != int(count)):
            raise ValueError("Each trials_per_group count must be a positive integer.")
    # Keep full IDs as a derived in-memory list for the existing per-trial loops.
    # Array counts follow trial_groups order; recordings are not discovered from disk.
    config["trials"] = [f"{group}{repeat}" for group, count in zip(groups, counts)
                        for repeat in range(1, int(count) + 1)]
    return config


def select_trials(config, reynolds, direction):
    """Return matching expanded trial IDs in configured group/repeat order."""
    return [trial_id for trial_id in config["trials"]
            if parse_trial(trial_id) == (int(reynolds), direction)]


def effective_event_time(metadata, event):
    """Use a manual event correction when set, otherwise the audio candidate.

    Python writes unset values as null; MATLAB may write []. Both mean that
    the candidate is used. Integration bounds never participate in this choice.
    """
    manual = metadata.get("manual_event_times_s", {}).get(event)
    value = metadata["audio_candidates_s"][event] if manual is None or manual == [] else manual
    if not np.isscalar(value) or not np.isfinite(value):
        raise ValueError(f"Missing/nonfinite {event} audio time or manual correction")
    return float(value)


def save_metadata(path, metadata):
    """Save the updated trial dictionary while retaining its other fields.

    Finite JSON values keep this file readable by both MATLAB and Python;
    callers change only their own selections, preserving manual audio edits.
    """
    path.write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")


def analyze_trial(trial_id, config, *, reselect=False, interactive=True):
    """Analyze one trial and return arrays/results for export and comparison.

    Reads {trial_id}_data.csv and {trial_id}_metadata.json from paths.data.
    Saved fit start and torque windows are reused unless reselect is True.
    Without interaction, a fit start must already be saved; missing integration
    windows use the effective event times plus configured window offsets.

    Failures raise to run_group, which reports the trial and continues. This
    function updates selections but leaves audio candidates/corrections intact.
    All numerical operations and selections stay in original video seconds.
    """
    # SECTION ONE - READ ANGULAR POSITION DATA AND REQUEST FIT LOCATION
    ###############################################################################################
    reynolds, direction = parse_trial(trial_id)
    data_dir = config["paths"]["data"]
    settings = config["analysis"]
    metadata_path = data_dir / f"{trial_id}_metadata.json"
    metadata = json.loads(metadata_path.read_text())
    # Generation marks a trial valid only after tracking and audio succeed.
    # Check that status even when an old angle CSV is still on disk.
    if metadata.get("valid") is not True:
        raise ValueError(f"Invalid generated trial: {metadata.get('failure_reason', '')}; "
                         f"frames {metadata.get('invalid_frames', [])}")

    fname = f"{trial_id}_data.csv"
    input_path = data_dir / fname
    print("------------------------------")
    print("Data coming from ", fname, " in directory ", data_dir)
    if settings["fit_enabled"]:
        print("Data is being fit!")
    print("Data incoming is in degrees and being converted to radians!!")

    # Read time t and angle phi, and find peaks of the data.
    full_t, full_phi, _, t_peaks, phi_peaks = read_data(input_path)
    # read_data converts degrees to radians. Report 1-based frame numbers to
    # match MATLAB and the numbered JPEGs; do not fill or discard bad samples.
    bad_frames = np.flatnonzero(~np.isfinite(full_phi)) + 1
    if bad_frames.size:
        raise ValueError(f"Nonfinite angles at frames {bad_frames.tolist()}")
    full_t, full_phi = validate_recording(full_t, full_phi)
    startup_time = effective_event_time(metadata, "startup")
    shutdown_time = effective_event_time(metadata, "shutdown")
    if not full_t[0] <= startup_time < shutdown_time <= full_t[-1]:
        raise ValueError(
            f"Invalid effective event times: startup={startup_time} s, "
            f"shutdown={shutdown_time} s; require "
            f"{full_t[0]} <= startup < shutdown <= {full_t[-1]} s."
        )

    # Reuse the user's free-decay selection before asking for another slider
    # choice. fit_segments keeps the existing nearest-peak fitting convention.
    t_fit = metadata.get("fit_start_s")
    if reselect or t_fit is None or t_fit == []:
        if not interactive:
            raise ValueError("Fit start is not saved; run interactively to select free decay.")
        # Reselection starts from the latest valid saved choice, not time zero.
        initial_time = (t_fit if np.isscalar(t_fit) and np.isfinite(t_fit)
                        and full_t[0] <= t_fit <= full_t[-1] else None)
        t_fit = float(plot_data(full_t, full_phi, 1, initial_time=initial_time))
    if not np.isscalar(t_fit) or not np.isfinite(t_fit) or not full_t[0] <= t_fit <= full_t[-1]:
        raise ValueError("Fit start must be within the original recording.")
    metadata["fit_start_s"] = float(t_fit)
    save_metadata(metadata_path, metadata)
    # SECTIONS TWO-SIX - RECONSTRUCT TORQUE BEFORE SELECTING INTEGRATION WINDOWS
    # The scientific narrative lives in scientific_pipeline; this function keeps
    # the original workflow order, metadata, interaction and parameter report.
    physics = config["physics"]
    reconstruction = reconstruct_torque(
        full_t, full_phi, t_peaks, phi_peaks, t_fit, settings, physics)
    fourier_t_seg = reconstruction["time"]
    true_torque = reconstruction["torque"]
    cummInt = reconstruction["cumulative"]
    phi_forward = reconstruction["phi_forward"]
    b, c = reconstruction["b"], reconstruction["c"]
    kappa, I = reconstruction["kappa"], reconstruction["inertia"]
    delta, xi = reconstruction["delta"], reconstruction["xi"]
    w0, T0, Td = reconstruction["w0"], reconstruction["T0"], reconstruction["Td"]
    print("Error of forward problem: ", reconstruction["error_forward"])

    # Retain the original parameter report and its existing damping convention.
    W = 28  # label column width
    print("=" * 52)
    print("SYSTEM PARAMETERS".center(52))
    print("=" * 52)
    print(f"{'Spring stiffness, κ':<{W}} {kappa:>8.2f}  dyn·cm")
    print(f"{'Moment of inertia, I':<{W}} {I:>8.2f}  g·cm²")
    print(f"{'Damping coefficient, δ':<{W}} {delta:>8.2f}  g·cm²·s⁻¹")
    print(f"{'Damping ratio, ξ':<{W}} {xi:>8.2f}  (dimensionless)")
    print(f"{'Natural frequency, ω₀':<{W}} {w0:>8.2f}  s⁻¹")
    print(f"{'Undamped natural period, T₀':<{W}} {T0:>8.2f}  s")
    print(f"{'Decay timescale, τ_d':<{W}} {Td:>8.2f}  s")
    print("=" * 52)

    # Audio initializes only missing windows. Existing integration selections
    # remain fixed even if a manual audio correction is edited later.
    saved_windows = metadata.get("integration_windows_s", {})
    windows = []
    needs_selection = reselect
    for event, event_time in (("startup", startup_time), ("shutdown", shutdown_time)):
        saved = saved_windows.get(event)
        if saved is None or saved == []:
            needs_selection = True
            saved = [event_time + offset
                     for offset in config["integration_window_offsets_s"][event]]
        windows.append(saved)
    # Invalid saved or audio-based suggestions can be repaired interactively.
    # Validation never clamps integration bounds or writes a replacement here.
    try:
        checked = [validate_window(fourier_t_seg, window) for window in windows]
        if checked[0][1] > checked[1][0]:
            raise ValueError("Startup and shutdown integration windows must not overlap.")
    except (TypeError, ValueError) as error:
        if not interactive:
            raise ValueError(f"Invalid integration window suggestions: {error}") from error
        needs_selection = True
    if interactive and needs_selection:
        windows = select_impulse_windows(fourier_t_seg, true_torque, *windows)

    # Finish the scientific calculation only after window selection/validation.
    impulse_results = measure_event_impulses(
        fourier_t_seg, true_torque, full_t, windows, physics, reynolds, direction)
    # Torque-window selection never modifies audio times or manual corrections.
    metadata["integration_windows_s"] = {
        event: [float(value) for value in window]
        for event, window in zip(("startup", "shutdown"), windows)
    }
    save_metadata(metadata_path, metadata)
    report_impulses(impulse_results)
    return dict(trial_id=trial_id, time=fourier_t_seg, torque=true_torque,
                cumulative=cummInt, startup_time=startup_time, impulses=impulse_results,
                phi=full_phi, phi_forward=phi_forward, b=b, c=c, inertia=I)


# SECTION SEVEN - SAVE SIGNAL AND PLOT RESULTS
###################################################################################################
# The per-trial analysis above feeds the following group plotting/export loop.

def plot_group_results(trials, output_dir, group):
    """Save angle, torque, cumulative-integral and individual impulse comparisons.

    Every curve uses its own sample grid shifted by effective startup time.
    CSV timestamps and integration selections are not changed by these plots.
    """
    for field, plotter, filename in (
        ("torque", plot_torque, "torque.png"),
        ("cumulative", plot_torque_int, "cumulative_integral.png"),
    ):
        fig, ax = plt.subplots()
        for trial in trials:
            plotter(trial["time"] - trial["startup_time"], trial[field],
                    ax=ax, label=trial["trial_id"])
            # Apply the same display-only time shift to the selected intervals.
            # Retain shaded windows and legend labels without numerical boxes.
            mark_impulse_windows(ax, trial["impulses"],
                                 time_offset_s=trial["startup_time"], trial_id=trial["trial_id"])
        ax.set_xlabel("Time from effective startup (s)")
        ax.set_title(f"{group}: {'cumulative torque' if field == 'cumulative' else field}")
        ax.legend(fontsize=7, loc="lower right")
        fig.tight_layout()
        fig.savefig(output_dir / filename)

    plot_jump_comparison(trials, output_dir, group)

    # Refresh this comparison too when a manual startup correction changes.
    fig, ax = plt.subplots()
    for trial in trials:
        ax.plot(trial["time"] - trial["startup_time"], np.rad2deg(trial["phi"]),
                label=trial["trial_id"])
    ax.set_xlabel("Time from effective startup (s)")
    ax.set_ylabel("Angular displacement (degrees)")
    ax.set_title(f"{group}: angular displacement")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "angular_displacement.png")


def plot_jump_comparison(trials, output_dir, group):
    """Compare saved startup/shutdown cumulative torque and annotate its errors."""
    # Paired observations share event colors; theory is a reference per event.
    # Keep sensitivity as the existing min/max range, not statistical uncertainty.
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(max(8, 1.25 * len(trials)), 6.5),
                             gridspec_kw={"height_ratios": [2, 1]})
    positions = np.arange(len(trials))
    net_impulses = np.zeros(len(trials))
    for event, color, marker, offset in (("startup", "tab:blue", "o", -.08),
                                          ("shutdown", "tab:orange", "s", .08)):
        rows = [next(row for row in trial["impulses"] if row["event"] == event)
                for trial in trials]
        measured = np.array([row["J_measured_dyn_cm_s"] for row in rows])
        net_impulses += measured
        x = positions + offset
        low = [row["sensitivity_min_dyn_cm_s"] for row in rows]
        high = [row["sensitivity_max_dyn_cm_s"] for row in rows]
        axes[0].vlines(x, low, high, color=color, linewidth=1.2, alpha=.55)
        axes[0].hlines(low, x - .025, x + .025, color=color, linewidth=1.2, alpha=.55)
        axes[0].hlines(high, x - .025, x + .025, color=color, linewidth=1.2, alpha=.55)
        axes[0].plot(x, measured, marker, color=color, markersize=6, zorder=3,
                     label=event.capitalize())
        for index, (trial, value) in enumerate(zip(trials, measured)):
            error = relative_error_metrics(trial["impulses"])[f"{event}_relative_error_pct"]
            label = "Error N/A" if error is None else f"Error {error:.2f}%"
            axes[0].annotate(label, (x[index], value), xytext=(0, -16 if event == "startup" else 10),
                             textcoords="offset points", ha="center",
                             va="top" if event == "startup" else "bottom", fontsize=9,
                             color=color, bbox=dict(facecolor="white", edgecolor="none", alpha=.8),
                             annotation_clip=False)
        theory = np.array([np.nan if row["J_theory_dyn_cm_s"] is None
                           else row["J_theory_dyn_cm_s"] for row in rows])
        if np.isfinite(theory).any():
            if np.isfinite(theory).all() and np.all(theory == theory[0]):
                axes[0].axhline(theory[0], color=color, linestyle="--", linewidth=1.3,
                               label=f"{event.capitalize()} theory")
            else:
                axes[0].plot(positions, theory, "--_", color=color, linewidth=1.3,
                             label=f"{event.capitalize()} theory")
    # This is only the sum of the two measured event impulses. Do not integrate
    # torque between events or infer uncertainty for their sum from sensitivity.
    axes[1].vlines(positions, 0, net_impulses, color="0.5", linewidth=1.2)
    axes[1].plot(positions, net_impulses, "o", color="0.2", markersize=6)
    for index, trial in enumerate(trials):
        mismatch = relative_error_metrics(trial["impulses"])["relative_magnitude_mismatch_pct"]
        label = "Mismatch N/A" if mismatch is None else f"Mismatch {mismatch:.3g}%"
        positive = net_impulses[index] >= 0
        axes[1].annotate(label, (positions[index], net_impulses[index]),
                         xytext=(0, -12 if positive else 10), textcoords="offset points",
                         ha="center", va="top" if positive else "bottom", fontsize=9,
                         bbox=dict(facecolor="white", edgecolor="none", alpha=.8),
                         annotation_clip=False)
    axes[1].set_title("Startup + shutdown (ideal: 0)", fontsize=10)
    for ax in axes:
        ax.axhline(0, color="0.35", linewidth=1, zorder=0)
        ax.set_axisbelow(True)
        ax.grid(axis="y", alpha=.18)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.margins(y=.2)
        # Center both panels on zero; the lower panel has its own scale.
        extent = max(abs(value) for value in ax.get_ylim())
        ax.set_ylim(-extent, extent)
    axes[0].set_ylabel("Cumulative torque (dyn·cm·s)")
    axes[1].set_ylabel("Net cumulative torque (dyn·cm·s)")
    axes[1].set_xticks(positions)
    axes[1].set_xticklabels([trial["trial_id"] for trial in trials],
                           rotation=45 if len(trials) > 8 else 0)
    axes[1].set_xlim(-.5, len(trials) - .5)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.5, .94),
               ncol=4, frameon=False, fontsize=9)
    fig.suptitle(f"{group}: cumulative jump comparison", fontsize=14)
    fig.text(.5, .015, "Error: relative error against theory. Mismatch: relative startup/shutdown magnitude mismatch.\n"
             "Thin ranges: window sensitivity (±0.25 s), not confidence intervals.",
             ha="center", fontsize=9, color="0.35")
    fig.tight_layout(rect=(0, .045, 1, .87))
    fig.savefig(output_dir / "cumulative_jump_comparison.png", dpi=180)



def run_group(config, reynolds, direction, *, reselect=False, interactive=True, show=True):
    """Process a selected group, export successes, and return (trials, failures).

    trials holds successful analysis dictionaries; failures maps each failed ID
    to its error message. One failed recording must not stop the other trials.
    Selection interaction and final diagnostic displays are separate options.
    """
    selected = select_trials(config, reynolds, direction)
    if not selected:
        raise ValueError(f"No configured trials for {reynolds}{direction}")
    group = f"{int(reynolds)}{direction}"
    output_dir = config["paths"]["output"] / group
    output_dir.mkdir(parents=True, exist_ok=True)
    # Replace only this group's analysis products, so failed trials cannot linger.
    for pattern in ("*_signal.csv", "*_jumps.csv", "jump_summary.csv",
                    "torque.png", "cumulative_integral.png", "cumulative_jump_comparison.png",
                    "angular_displacement.png", "*_jump_diagnostics.png",
                    "relative_errors.png", "relative_errors.csv",
                    "*_impulses.csv", "impulse_summary.csv", "impulse_comparison.png",
                    "*_impulse_diagnostics.png"):
        for path in output_dir.glob(pattern):
            path.unlink()
    trials, failures = [], {}
    for trial_id in selected:
        try:
            trial = analyze_trial(trial_id, config, reselect=reselect, interactive=interactive)
            np.savetxt(output_dir / f"{trial_id}_signal.csv",
                       np.column_stack([trial["time"], trial["torque"]]), delimiter=",",
                       header="time_seconds,torque_dyn_cm", comments="")
            save_impulse_summary(output_dir / f"{trial_id}_jumps.csv",
                                 trial_id, trial["impulses"])
        except Exception as exc:
            # An export may have failed after writing one file. Remove any
            # partial products before excluding this trial from the summary.
            for suffix in ("signal", "jumps"):
                (output_dir / f"{trial_id}_{suffix}.csv").unlink(missing_ok=True)
            failures[trial_id] = str(exc)
            print(f"FAILED {trial_id}: {exc}")
            continue
        trials.append(trial)
    if trials:
        # Keep two rows per successful trial, retaining its identity on both.
        # Passing None lets the existing exporter use those individual IDs.
        records = [dict(run_identity=trial["trial_id"], **row)
                   for trial in trials for row in trial["impulses"]]
        save_impulse_summary(output_dir / "jump_summary.csv", None, records)
        error_records = [dict(trial_id=trial["trial_id"], **relative_error_metrics(trial["impulses"]))
                         for trial in trials]
        with (output_dir / "relative_errors.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(error_records[0]))
            writer.writeheader()
            writer.writerows(error_records)
        plot_group_results(trials, output_dir, group)
        for trial in trials:
            if show:
                plot_phi_gen(trial["time"], trial["phi"], trial["phi_forward"])
            # Save the startup/shutdown close-ups and their local cumulative
            # impulses for every trial, including runs using --no-show.
            fig = plot_impulse_diagnostics(trial["time"], trial["torque"], trial["impulses"])
            fig.suptitle(f"{trial['trial_id']}: startup and shutdown cumulative torque")
            fig.savefig(output_dir / f"{trial['trial_id']}_jump_diagnostics.png", dpi=300)
            if not show:
                plt.close(fig)
        if show:
            plt.show()
        else:
            plt.close("all")
    print(f"{group}: {len(trials)} succeeded, {len(failures)} failed. Outputs: {output_dir}")
    return trials, failures


def main():
    """Translate command-line choices into one configured group analysis."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path(__file__).with_name("experiment_config.json"))
    parser.add_argument("reynolds_pos", type=int, nargs="?", help="Reynolds number, e.g. 1000")
    parser.add_argument("direction_pos", choices=("f", "r"), nargs="?", help="Flow direction")
    parser.add_argument("--re", type=int, dest="reynolds")
    parser.add_argument("--direction", choices=("f", "r"))
    parser.add_argument("--reselect", action="store_true", help="Pick fit start and torque windows again")
    parser.add_argument("--non-interactive", action="store_true", help="Require saved fit; use saved/default windows")
    parser.add_argument("--no-show", action="store_true", help="Save comparison plots without displaying diagnostics")
    args = parser.parse_args()
    # Accept the shorter positional form as well as the original named options.
    # Reject duplicate values so the requested group is always unambiguous.
    for name in ("reynolds", "direction"):
        positional = getattr(args, "reynolds_pos" if name == "reynolds" else "direction_pos")
        if positional is not None:
            if getattr(args, name) is not None:
                parser.error(f"Supply {name} only once, positionally or with its named option")
            setattr(args, name, positional)
    if args.reynolds is None or args.direction is None:
        parser.error("Supply Reynolds number and direction, e.g. python main.py 1000 f")
    # --no-show affects final diagnostics only. --non-interactive also disables
    # the selection pickers, so it cannot be combined with forced reselection.
    if args.reselect and args.non_interactive:
        parser.error("--reselect requires interactive selection")
    config = load_config(args.config)
    trials, failures = run_group(config, args.reynolds, args.direction,
                                reselect=args.reselect, interactive=not args.non_interactive,
                                show=config["analysis"]["display_plots"] and not args.no_show)
    # Report partial failures to shell callers even when valid outputs exist.
    return 0 if trials and not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
