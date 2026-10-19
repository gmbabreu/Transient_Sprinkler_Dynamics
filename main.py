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
    process transparency. Sections 1-6 run once per trial in analyze_trial;
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
import json
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
from scipy.signal import cheby2, sosfiltfilt

# DATA READ FUNCS
from Main_Functions.dataread_funcs import read_data, plot_data, fit_segments

# ESTIMATE FUNCS
from Main_Functions.estimate_funcs import est_a, est_b, get_constants, fit_phi

# PROCESS FUNCS
from Main_Functions.process_funcs import combine_data, validate_recording

# PLOT FUNCS
from Main_Functions.plot_funcs import plot_torque, plot_torque_int, plot_phi_gen

# IMPULSE FUNCS
from Main_Functions.impulse_funcs import (
    analyze_impulses, report_impulses, save_impulse_summary,
    mark_impulse_windows, plot_impulse_diagnostics, select_impulse_windows, validate_window,
)

# FOURIER TRANSFORM FUNCS
from Main_Functions.fft_funcs import torque_solver, phi_from_torque


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
        t_fit = float(plot_data(full_t, full_phi, 1))
    if not np.isscalar(t_fit) or not np.isfinite(t_fit) or not full_t[0] <= t_fit <= full_t[-1]:
        raise ValueError("Fit start must be within the original recording.")
    metadata["fit_start_s"] = float(t_fit)
    save_metadata(metadata_path, metadata)
    # Return the segment of t and phi after the chosen time for fitting.
    fit_index, fit_peaks_index, t_fit, t_seg, phi_seg = fit_segments(
        full_t, full_phi, t_peaks, t_fit)

    # SECTION TWO - ESTIMATE OMEGA AND GAMMA
    ###############################################################################################
    # Preserve the original peak estimates and optional damped-oscillator fit.
    # b is the fitted decay rate and c is the natural-frequency parameter.
    # Estimate a from data.
    a_est = est_a(fit_peaks_index, t_fit, t_peaks)
    # Estimate b from data.
    b_est = est_b(fit_peaks_index, t_peaks, phi_peaks)
    # Estimate c from a and b.
    c_est = np.sqrt((a_est*a_est) + (b_est*b_est))
    # Get the ODE constants using the estimates for gamma and omega.
    C1_est, C2_est = get_constants(b_est, c_est, full_phi[fit_index])
    if settings["fit_enabled"]:
        # Fit the data.
        b, c, C_1, C_2 = fit_phi(t_seg, phi_seg, b_est, c_est, C1_est, C2_est,
                                full_t[fit_index])
    else:
        # Set end values as estimates if no fit is requested.
        b = b_est
        c = c_est
        C_1 = C1_est
        C_2 = C2_est

    # Define the analytical solution in terms of the newly found values.
    def phi_an(t):
        """Retain the analytical fit accepted by the existing padding API.

        combine_data currently prepends constant initial-value padding; it
        does not append this analytical curve or replace measured samples.
        """
        t0 = t_peaks[fit_peaks_index]
        a = np.sqrt(c**2 - b**2)
        return np.exp(-b * (t - t0)) * (C_1 * np.cos(a * (t - t0))
                                        + C_2 * np.sin(a * (t - t0)))

    # SECTION THREE - PREPARE DATA WITH INITIAL MEASUREMENT NOISE RETAINED
    ###############################################################################################
    # The former remove_noise call is bypassed; keep the original FFT padding.
    fourier_t, fourier_phi, N_added = combine_data(full_t, full_phi, phi_an)

    # SECTION FOUR - FOURIER TRANSFORM
    ###############################################################################################
    N_f = len(fourier_phi)
    # Keep the existing FFT period and normalized Chebyshev cutoff convention;
    # batching changes which trials run, not the numerical reconstruction.
    L = fourier_t[-1] - fourier_t[0]
    # Output the extracted torque signal, then retain its real part.
    extracted_torque = torque_solver(N_f, b, c, L, fourier_phi)
    extracted_torque = np.real(extracted_torque)
    if settings["lowpass_enabled"]:
        sos = cheby2(N=settings["cheby2_order"], rs=settings["cheby2_rs"],
                     Wn=settings["cheby2_wn"], btype="low", output="sos")
        extracted_torque = sosfiltfilt(sos, extracted_torque)

    # Remove exactly the prepended padding. These timestamps then match the
    # input recording, so exports and integration use only measured samples.
    fourier_t_seg = fourier_t[N_added:]
    fourier_phi_seg = fourier_phi[N_added:]
    extracted_torque_seg = extracted_torque[N_added:]

    # SECTION FIVE - DO FORWARD PROBLEM WITH TORQUE SIGNAL
    ###############################################################################################
    phi_forward = phi_from_torque(N_f, fourier_t_seg, extracted_torque_seg, b, c)
    # SciPy's trapezoidal rule also supports environments using NumPy 1.x.
    error_forward = (np.sqrt(integrate.trapezoid((phi_forward - fourier_phi_seg)**2,
                                        x=fourier_t_seg))
                     / np.sqrt(integrate.trapezoid(fourier_phi_seg**2, x=fourier_t_seg)))
    print("Error of forward problem: ", error_forward)

    # SECTION SIX - COMPUTE TORQUE INTEGRALS
    ###############################################################################################
    physics = config["physics"]
    kappa = physics["kappa_dyn_cm_per_rad"]
    # The solver returns forcing per unit inertia. Scale by I to obtain
    # dimensional torque in dyn cm, using the same relation as before.
    I = kappa/c**2  # solve for inertia, in units of g cm^2
    delta = b*I  # solve for damping coefficient, in units of g cm^2 s^-1
    xi = delta/(2*np.sqrt(kappa*I))  # solve for damping ratio, dimensionless
    w0 = (delta/I)/(2*xi)  # solve for natural frequency, units of s^-1
    T0 = (2*np.pi)/w0  # solve for natural, undamped period, units of s
    Td = 1/(xi*w0)  # solve for decay timescale, units of s

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

    true_torque = extracted_torque_seg * I
    # This full-record cumulative curve is for plotting. The two measured
    # event impulses below are integrated separately within their own bounds.
    cummInt = integrate.cumulative_trapezoid(np.real(true_torque), fourier_t_seg, initial=0)

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

    flow_factor = physics["flow_rate_per_re_cm3_s"]
    # Reuse the existing signed integration, theory and sensitivity code for
    # startup and shutdown, with no new baseline or between-event integral.
    impulse_results = analyze_impulses(
        fourier_t_seg, true_torque, full_t,
        startup_window=windows[0], shutdown_window=windows[1],
        rho_g_cm3=physics["rho_g_cm3"],
        flow_rate_cm3_s=None if flow_factor is None else reynolds * flow_factor,
        geometry_factor_cm2=physics["geometry_factor_cm2"],
        sign_convention=physics["sign_convention"], spin_dir=direction,
    )
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
            # Multiple trials retain their shaded windows and legend labels;
            # numerical annotation boxes are shown only for a single trial.
            mark_impulse_windows(ax, trial["impulses"], annotate=len(trials) == 1,
                                 time_offset_s=trial["startup_time"], trial_id=trial["trial_id"])
        ax.set_xlabel("Time from effective startup (s)")
        ax.set_title(f"{group}: {field}")
        ax.legend(fontsize=7, loc="lower right")
        fig.tight_layout()
        fig.savefig(output_dir / filename)

    fig, ax = plt.subplots()
    for index, trial in enumerate(trials):
        # Small horizontal offsets keep individual trials visible at each
        # event category. The measured signed impulse values are unchanged.
        offset = (index - (len(trials) - 1) / 2) * .6 / max(1, len(trials))
        positions = np.array([0, 1]) + offset
        measured = [row["J_measured_dyn_cm_s"] for row in trial["impulses"]]
        line, = ax.plot(positions, measured, "o", label=trial["trial_id"])
        for position, row in zip(positions, trial["impulses"]):
            # Vertical lines show the existing integration-window sensitivity
            # range, not statistical uncertainty; black ticks mark theory.
            ax.vlines(position, row["sensitivity_min_dyn_cm_s"],
                      row["sensitivity_max_dyn_cm_s"], color=line.get_color())
            if row["J_theory_dyn_cm_s"] is not None:
                ax.plot(position, row["J_theory_dyn_cm_s"], "k_", markersize=14)
    ax.plot([], [], "k_", markersize=14, label="Theory")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Startup", "Shutdown"])
    ax.set_ylabel("Signed impulse (dyn·cm·s)")
    ax.set_title(f"{group}: individual impulses (lines: window sensitivity)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "impulse_comparison.png")

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
    for pattern in ("*_signal.csv", "*_impulses.csv", "impulse_summary.csv",
                    "torque.png", "cumulative_integral.png", "impulse_comparison.png",
                    "angular_displacement.png", "*_impulse_diagnostics.png"):
        for path in output_dir.glob(pattern):
            path.unlink()
    trials, failures = [], {}
    for trial_id in selected:
        try:
            trial = analyze_trial(trial_id, config, reselect=reselect, interactive=interactive)
            np.savetxt(output_dir / f"{trial_id}_signal.csv",
                       np.column_stack([trial["time"], trial["torque"]]), delimiter=",",
                       header="time_seconds,torque_dyn_cm", comments="")
            save_impulse_summary(output_dir / f"{trial_id}_impulses.csv",
                                 trial_id, trial["impulses"])
        except Exception as exc:
            # An export may have failed after writing one file. Remove any
            # partial products before excluding this trial from the summary.
            for suffix in ("signal", "impulses"):
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
        save_impulse_summary(output_dir / "impulse_summary.csv", None, records)
        plot_group_results(trials, output_dir, group)
        for trial in trials:
            if show:
                plot_phi_gen(trial["time"], trial["phi"], trial["phi_forward"])
            # Save the startup/shutdown close-ups and their local cumulative
            # impulses for every trial, including runs using --no-show.
            fig = plot_impulse_diagnostics(trial["time"], trial["torque"], trial["impulses"])
            fig.suptitle(f"{trial['trial_id']}: startup and shutdown impulses")
            fig.savefig(output_dir / f"{trial['trial_id']}_impulse_diagnostics.png", dpi=300)
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
