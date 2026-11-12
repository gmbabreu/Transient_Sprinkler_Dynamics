"""
Raw transient-impulse analysis for validating the sprinkler theory.

This module compares the measured angular impulse around pump startup and
shutdown with the theoretical prediction

    J_start = -s * rho * Q * G
    J_stop  = +s * rho * Q * G

where

    rho : fluid density [g/cm^3]
    Q   : signed total volume flow rate through both arms [cm^3/s]
    G   : signed geometry factor for one arm [cm^2]
    s   : sign-convention multiplier, either +1 or -1

The measured impulse is

    J = integral tau(t) dt

over a user-selected time window around each pump event.

For now, this file deliberately uses only the raw reconstructed torque.
No baseline subtraction or steady-flow correction is applied.
"""

import csv

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid


def validate_signal(time_s, torque_dyn_cm):
    """
    Validate and return the time and torque arrays.

    Requirements
    ------------
    - Both arrays must be one-dimensional and have the same length.
    - Values must be finite.
    - Time must be strictly increasing.
    - Torque must be real-valued.
    """
    time_s = np.asarray(time_s, dtype=float)
    torque_dyn_cm = np.asarray(torque_dyn_cm)

    valid = (
        time_s.ndim == 1
        and len(time_s) >= 2
        and torque_dyn_cm.shape == time_s.shape
        and not np.iscomplexobj(torque_dyn_cm)
        and np.all(np.isfinite(time_s))
        and np.all(np.isfinite(torque_dyn_cm))
        and np.all(np.diff(time_s) > 0)
    )

    if not valid:
        raise ValueError(
            "Time and torque must be finite, matching 1-D arrays, "
            "with strictly increasing time."
        )

    return time_s, torque_dyn_cm.astype(float)


def validate_window(time_s, window):
    """
    Check that an integration window [start, end] lies inside the data.
    """
    if len(window) != 2:
        raise ValueError("An integration window must contain exactly two endpoints.")

    start_s, end_s = window

    if (
        not np.all(np.isfinite([start_s, end_s]))
        or not time_s[0] <= start_s < end_s <= time_s[-1]
    ):
        raise ValueError(
            f"Invalid integration window {window}. "
            f"Data span [{time_s[0]}, {time_s[-1]}] s."
        )

    return float(start_s), float(end_s)


def window_samples(time_s, torque_dyn_cm, window):
    """
    Return the torque signal restricted to a time window.

    If a requested window endpoint lies between measured samples, the torque
    at that endpoint is obtained by linear interpolation. This lets the
    integration use the exact user-requested time limits without extrapolation.
    """
    time_s, torque_dyn_cm = validate_signal(time_s, torque_dyn_cm)
    start_s, end_s = validate_window(time_s, window)

    local_time = np.r_[
        start_s,
        time_s[(time_s > start_s) & (time_s < end_s)],
        end_s,
    ]
    local_torque = np.interp(local_time, time_s, torque_dyn_cm)

    return local_time, local_torque


def integrate_impulse(time_s, torque_dyn_cm, window):
    """
    Compute the signed angular impulse over one event window.

    Returns
    -------
    float
        Integral of torque over time in dyn*cm*s, equivalently g*cm^2/s.
    """
    local_time, local_torque = window_samples(time_s, torque_dyn_cm, window)
    return float(trapezoid(local_torque, local_time))


def theoretical_impulses(
    rho_g_cm3,
    flow_rate_cm3_s,
    geometry_factor_cm2,
    spin_dir,
    sign_convention,
):
    """
    Compute the theoretical startup and shutdown impulses.

    Theory
    ------
    Let Q be positive for forward/outward flow and negative for reverse/inward
    flow. Then

        J_start = -s * rho * Q * G
        J_stop  = +s * rho * Q * G

    The supplied flow_rate_cm3_s is the positive magnitude of the TOTAL flow
    through both arms. Therefore no extra factor of two is introduced here.

    Returns
    -------
    (J_start, J_stop, Q_signed)
        If either flow rate or geometry factor is unavailable, all theoretical
        impulse values are returned as None while Q_signed is still returned
        when possible.
    """
    if spin_dir not in ("f", "r"):
        raise ValueError(
            "spin_dir must be 'f' for forward/outward flow "
            "or 'r' for reverse/inward flow."
        )

    if sign_convention not in (-1, 1):
        raise ValueError("sign_convention must be explicitly +1 or -1.")

    if not np.isfinite(rho_g_cm3) or rho_g_cm3 <= 0:
        raise ValueError("Fluid density must be finite and positive.")

    if flow_rate_cm3_s is not None:
        if not np.isfinite(flow_rate_cm3_s) or flow_rate_cm3_s <= 0:
            raise ValueError(
                "flow_rate_cm3_s must be a positive total-flow magnitude or None."
            )

    if geometry_factor_cm2 is not None and not np.isfinite(geometry_factor_cm2):
        raise ValueError("geometry_factor_cm2 must be finite or None.")

    if flow_rate_cm3_s is None:
        q_signed = None
    else:
        q_signed = flow_rate_cm3_s if spin_dir == "f" else -flow_rate_cm3_s

    if q_signed is None or geometry_factor_cm2 is None:
        return None, None, q_signed

    theory_scale = (
        sign_convention
        * rho_g_cm3
        * q_signed
        * geometry_factor_cm2
    )

    return -theory_scale, theory_scale, q_signed


def _valid_sensitivity_windows(
    time_s,
    chosen_window,
    other_window,
    changes_s=(-0.25, 0.0, 0.25),
):
    """
    Generate nearby windows for checking sensitivity to integration limits.

    Each positive change expands the window equally on both sides.
    Each negative change contracts it equally on both sides.
    Defaults check a 0.25-second contraction, the original, and a 0.25-second expansion.
    Invalid windows are discarded.
    """
    candidates = []

    for change in changes_s:
        start_s = chosen_window[0] - change
        end_s = chosen_window[1] + change

        inside_data = time_s[0] <= start_s < end_s <= time_s[-1]
        avoids_other_event_window = (
            end_s <= other_window[0] or start_s >= other_window[1]
        )

        if inside_data and avoids_other_event_window:
            candidates.append((start_s, end_s))

    return candidates


def analyze_impulses(
    time_s,
    torque_dyn_cm,
    full_time_s,
    *,
    startup_window,
    shutdown_window,
    rho_g_cm3=1.0,
    flow_rate_cm3_s=None,
    geometry_factor_cm2=None,
    sign_convention=1,
    spin_dir="f",
):
    """
    Measure startup/shutdown impulses and compare them with theory.

    This function is intentionally limited to the RAW reconstructed torque.
    No background subtraction or steady-flow correction is performed.

    Parameters
    ----------
    time_s, torque_dyn_cm
        Reconstructed dimensional torque signal and its time coordinate.
    full_time_s
        Original experimental time array. It must match time_s exactly so the
        integration windows are interpreted in the original experiment coordinates.
    startup_window, shutdown_window
        Integration windows around the two transient events.
    rho_g_cm3
        Fluid density.
    flow_rate_cm3_s
        Positive magnitude of the TOTAL flow rate through both arms.
    geometry_factor_cm2
        Signed geometry factor G for ONE arm, traversed outward from the hub.
    sign_convention
        +1 or -1 mapping the experimental positive-rotation convention onto
        the positive z-direction used to define G.
    spin_dir
        'f' for forward/outward flow or 'r' for reverse/inward flow.

    Returns
    -------
    list of dict
        One result dictionary for startup and one for shutdown.
    """
    time_s, torque_dyn_cm = validate_signal(time_s, torque_dyn_cm)

    full_time_s = np.asarray(full_time_s, dtype=float)
    full_time_s, _ = validate_signal(full_time_s, np.zeros_like(full_time_s))

    if time_s.shape != full_time_s.shape or not np.array_equal(time_s, full_time_s):
        raise ValueError(
            "The torque time array must match the original experimental time array."
        )

    startup_window = validate_window(time_s, startup_window)
    shutdown_window = validate_window(time_s, shutdown_window)

    if startup_window[1] > shutdown_window[0]:
        raise ValueError("Startup and shutdown integration windows must not overlap.")

    theory_start, theory_stop, q_signed = theoretical_impulses(
        rho_g_cm3=rho_g_cm3,
        flow_rate_cm3_s=flow_rate_cm3_s,
        geometry_factor_cm2=geometry_factor_cm2,
        spin_dir=spin_dir,
        sign_convention=sign_convention,
    )

    event_specs = [
        ("startup", startup_window, shutdown_window, theory_start),
        ("shutdown", shutdown_window, startup_window, theory_stop),
    ]

    results = []

    for event, window, other_window, theory in event_specs:
        measured = integrate_impulse(time_s, torque_dyn_cm, window)

        sensitivity_windows = _valid_sensitivity_windows(
            time_s,
            window,
            other_window,
        )
        sensitivity_values = [
            integrate_impulse(time_s, torque_dyn_cm, candidate)
            for candidate in sensitivity_windows
        ]
        # Keep the two named changes so the report can show contraction first,
        # even when it produces a larger signed impulse than expansion does.
        contracted = next((value for candidate, value in zip(sensitivity_windows, sensitivity_values)
                           if candidate[0] > window[0] and candidate[1] < window[1]), None)
        expanded = next((value for candidate, value in zip(sensitivity_windows, sensitivity_values)
                         if candidate[0] < window[0] and candidate[1] > window[1]), None)

        residual = None if theory is None else measured - theory
        ratio = None if theory is None or theory == 0 else measured / theory

        results.append(
            {
                "event": event,
                "window_start_s": window[0],
                "window_end_s": window[1],
                "J_measured_dyn_cm_s": measured,
                "J_theory_dyn_cm_s": theory,
                "residual_dyn_cm_s": residual,
                "measured_over_theory": ratio,
                "sensitivity_contracted_dyn_cm_s": contracted,
                "sensitivity_expanded_dyn_cm_s": expanded,
                "sensitivity_min_dyn_cm_s": min(sensitivity_values),
                "sensitivity_max_dyn_cm_s": max(sensitivity_values),
                "sensitivity_windows_s": "; ".join(
                    f"{a:g}:{b:g}" for a, b in sensitivity_windows
                ),
                "rho_g_cm3": rho_g_cm3,
                "flow_rate_cm3_s": flow_rate_cm3_s,
                "Q_signed_cm3_s": q_signed,
                "geometry_factor_cm2": geometry_factor_cm2,
                "sign_convention": sign_convention,
                "spin_dir": spin_dir,
            }
        )

    return results


def relative_error_metrics(results):
    """Return the existing reported percentages for plots, tables and reports.

    Event error = 100*abs(J_measured-J_theory)/abs(J_theory).
    Magnitude mismatch = 100*abs(abs(J_start)-abs(J_stop))/max(abs(J_start),abs(J_stop)).
    Undefined denominators or missing theory return None, rather than zero.
    """
    metrics = {}
    events = {row["event"]: row for row in results}
    for event in ("startup", "shutdown"):
        row = events.get(event)
        theory = None if row is None else row["J_theory_dyn_cm_s"]
        metrics[f"{event}_relative_error_pct"] = (
            None if theory is None or theory == 0 else
            100 * (abs(row["residual_dyn_cm_s"]) / abs(theory)))
    metrics["relative_magnitude_mismatch_pct"] = None
    if "startup" in events and "shutdown" in events:
        j_start = events["startup"]["J_measured_dyn_cm_s"]
        j_stop = events["shutdown"]["J_measured_dyn_cm_s"]
        larger = max(abs(j_start), abs(j_stop))
        if larger > 0:
            metrics["relative_magnitude_mismatch_pct"] = (
                100 * (abs(abs(j_start) - abs(j_stop)) / larger))
    return metrics


def report_impulses(results):
    """
    Print a concise experimental-versus-theory comparison.
    """
    print("\nRAW TRANSIENT CUMULATIVE TORQUE COMPARISON")
    print("Units: dyn*cm*s = g*cm^2/s")
    print("Measured cumulative torque = integral of reconstructed torque over the selected window.")
    print("No baseline or steady-flow correction is applied.")

    metrics = relative_error_metrics(results)
    for row in results:
        print(
            f"\n{row['event'].upper()} "
            f"[{row['window_start_s']:g}, {row['window_end_s']:g}] s"
        )
        print(f"  measured J = {row['J_measured_dyn_cm_s']:.6g}")

        if row["J_theory_dyn_cm_s"] is None:
            print("  theory     = unavailable (flow rate and/or geometry factor missing)")
        else:
            theory = row["J_theory_dyn_cm_s"]
            print(f"  theory J   = {theory:.6g}")
            if theory == 0:
                print("  relative error = undefined (zero theoretical cumulative torque)")
            else:
                # Compare signed impulses against the theoretical magnitude,
                # so an opposite-sign measurement still counts as disagreement.
                relative_error = metrics[f"{row['event']}_relative_error_pct"]
                print(f"  relative error = {relative_error:.2f}%")

        contracted = row["sensitivity_contracted_dyn_cm_s"]
        expanded = row["sensitivity_expanded_dyn_cm_s"]
        contracted_text = "unavailable" if contracted is None else f"{contracted:.6g}"
        expanded_text = "unavailable" if expanded is None else f"{expanded:.6g}"
        print("  window sensitivity (contracted, expanded) = "
              f"[{contracted_text}, {expanded_text}]")

    if len(results) == 2:
        j_start = results[0]["J_measured_dyn_cm_s"]
        j_stop = results[1]["J_measured_dyn_cm_s"]

        print("\nSTARTUP/SHUTDOWN CONSISTENCY")
        print(
            "  J_start + J_stop = "
            f"{j_start + j_stop:.6g} "
            "(theory predicts approximately 0 for equal-and-opposite kicks)"
        )

        larger = max(abs(j_start), abs(j_stop))
        if larger > 0:
            relative_magnitude_mismatch = metrics["relative_magnitude_mismatch_pct"]
            print(
                "  relative magnitude mismatch = "
                f"{relative_magnitude_mismatch:.3g}%"
            )


def save_impulse_summary(path, run_identity, results):
    """
    Save the raw impulse/theory comparison to a CSV file.

    Supply a trial ID to attach it to that trial's startup/shutdown rows.
    For a group summary, supply None and rows that already contain each
    trial's run_identity; the individual result columns remain identical.
    """
    # Group summaries already carry an identity on every individual result.
    records = results if run_identity is None else [dict(run_identity=run_identity, **row) for row in results]

    with open(path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def mark_impulse_windows(ax, results, *, time_offset_s=0.0, trial_id=None):
    """
    Mark selected integration windows on an existing axis.

    Subtract time_offset_s for an aligned comparison plot while keeping the
    results in original recording seconds. trial_id identifies each window
    in a group legend. Existing individual plots keep their default behavior.
    """
    colors = ("tab:blue", "tab:orange")

    for row, color in zip(results, colors):
        event_label = f"{trial_id}: {row['event']}" if trial_id else row['event']
        ax.axvspan(
            row["window_start_s"] - time_offset_s,
            row["window_end_s"] - time_offset_s,
            color=color,
            alpha=0.15,
            label=f"{event_label} integration window",
        )

    ax.legend(fontsize=7, loc="lower right")


def select_impulse_windows(
    time_s, torque_dyn_cm, startup_window, shutdown_window, *, view_windows=None
):
    """Pick two endpoints per event on extracted torque close-ups.

    Close the figure, press Enter, or click Done to accept completed pairs.
    Events with no complete pair retain a valid supplied window. Invalid
    suggestions must be corrected before accepting; closing or cancelling
    without valid, nonoverlapping windows raises ValueError.
    Escape discards all clicks. view_windows controls only the initial views.
    """
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Button

    time_s, torque_dyn_cm = validate_signal(time_s, torque_dyn_cm)
    defaults = []
    for window in (startup_window, shutdown_window):
        try:
            defaults.append(validate_window(time_s, window))
        except (TypeError, ValueError):
            # No clamped fallback: this event requires a manual replacement.
            defaults.append(None)
    defaults_valid = (all(window is not None for window in defaults)
                      and defaults[0][1] <= defaults[1][0])
    if view_windows is None:
        # Include startup from 9.5 s and shutdown from 19.5 s, or earlier
        # when the 3 s lead-in requires it. Keep 2 s after each suggestion.
        # These are display limits only; integration bounds stay unchanged.
        view_windows = [(min(window[0] - 3, start), window[1] + 2)
                        if window is not None else (time_s[0], time_s[-1])
                        for window, start in zip(defaults, (9.5, 19.5))]
    if len(view_windows) != 2:
        raise ValueError("Supply one close-up view for each event.")

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    fig.subplots_adjust(top=0.76, bottom=0.25, wspace=0.3)
    fig.suptitle(
        "Extracted torque: click two interval endpoints in each panel (either order).\n"
        "Done / Enter / close: accept valid windows; invalid suggestions need correction.\n"
        "Reset / Escape: restore suggestions. Turn off zoom/pan to pick.",
        fontsize=10,
    )
    names = ("Startup", "Shutdown")
    colors = ("tab:blue", "tab:orange")
    chosen = list(defaults)
    pending = [None, None]
    manual = [False, False]
    markers = [[], []]
    status = fig.text(0.05, 0.14, "", fontsize=9)

    def checked_windows():
        for name, window in zip(names, chosen):
            if window is None:
                raise ValueError(f"{name} integration window is invalid; select two valid endpoints.")
            validate_window(time_s, window)
        if chosen[0][1] > chosen[1][0]:
            raise ValueError("Startup and shutdown integration windows must not overlap.")
        return tuple(chosen)

    for i, (ax, view) in enumerate(zip(axes, view_windows)):
        # Only the displayed context may be clipped to the recording. This
        # never changes a suggested or selected integration bound.
        if defaults[i] is not None:
            view = (min(view[0], defaults[i][0]), max(view[1], defaults[i][1]))
        try:
            context = validate_window(time_s, (max(time_s[0], view[0]), min(time_s[-1], view[1])))
        except ValueError:
            context = (time_s[0], time_s[-1])
        local_time, local_torque = window_samples(time_s, torque_dyn_cm, context)
        ax.plot(local_time, local_torque, color="hotpink")
        ax.set(xlim=context, xlabel="Experimental time (s)",
               ylabel="Extracted torque (dyn*cm)")
        ax.grid(alpha=0.25)

    def redraw(message=""):
        for i, ax in enumerate(axes):
            for artist in markers[i]:
                artist.remove()
            markers[i] = []
            if chosen[i] is None:
                ax.set_title(f"{names[i]}: select a valid integration window")
            else:
                a, b = chosen[i]
                markers[i] = [ax.axvspan(a, b, color=colors[i], alpha=0.15),
                              ax.axvline(a, color=colors[i], linestyle="--"),
                              ax.axvline(b, color=colors[i], linestyle="--")]
                source = "clicked" if manual[i] else "suggested"
                ax.set_title(f"{names[i]}: [{a:.3f}, {b:.3f}] s ({source})")
            if pending[i] is not None:
                markers[i].append(ax.axvline(pending[i], color=colors[i]))
        if not message:
            try:
                checked_windows()
                message = "Shaded windows will be used. Click Done when ready."
            except ValueError as error:
                message = str(error)
        status.set_text(message)
        fig.canvas.draw_idle()

    def on_click(event):
        if event.button != 1 or event.inaxes not in axes or event.xdata is None:
            return
        toolbar = getattr(fig.canvas, "toolbar", None)
        if toolbar is not None and toolbar.mode:
            return
        i = list(axes).index(event.inaxes)
        endpoint = float(event.xdata)
        if not np.isfinite(endpoint) or not time_s[0] <= endpoint <= time_s[-1]:
            return
        if pending[i] is None:
            # Keep the existing fallback behavior when suggestions are valid.
            trial = list(chosen)
            trial[i] = defaults[i]
            if defaults_valid and trial[0][1] > trial[1][0]:
                redraw("Reset first: the initial window would overlap the other clicked window.")
                return
            chosen[i], manual[i], pending[i] = defaults[i], False, endpoint
            redraw(f"{names[i]}: first endpoint {endpoint:.3f} s; click the other endpoint.")
            return
        candidate = tuple(sorted((pending[i], endpoint)))
        trial = list(chosen)
        trial[i] = candidate
        try:
            validate_window(time_s, candidate)
            # An invalid initial suggestion must not block repair of the other
            # event. Enforce overlap against accepted clicks or valid defaults.
            if (all(window is not None for window in trial)
                    and (defaults_valid or manual[1-i]) and trial[0][1] > trial[1][0]):
                raise ValueError("Startup must end before shutdown begins.")
        except ValueError as error:
            redraw(f"{error} Click a different second endpoint, or Reset.")
            return
        chosen[i], manual[i], pending[i] = candidate, True, None
        redraw()

    def reset(event=None):
        chosen[:] = defaults
        pending[:] = [None, None]
        manual[:] = [False, False]
        redraw()

    def done(event=None):
        try:
            checked_windows()
        except ValueError as error:
            redraw(str(error))
            return
        plt.close(fig)

    def on_key(event):
        if event.key == "enter":
            done()
        elif event.key == "escape":
            reset()
            plt.close(fig)

    reset_button = Button(fig.add_axes([0.61, 0.04, 0.16, 0.06]), "Reset selections")
    done_button = Button(fig.add_axes([0.79, 0.04, 0.16, 0.06]), "Done")
    reset_button.on_clicked(reset)
    done_button.on_clicked(done)
    connections = [fig.canvas.mpl_connect("button_press_event", on_click),
                   fig.canvas.mpl_connect("key_press_event", on_key)]
    redraw()
    try:
        plt.show(block=True)
    finally:
        for connection in connections:
            fig.canvas.mpl_disconnect(connection)
        plt.close(fig)
    # Closing/cancelling invalid suggestions fails rather than returning
    # clamped or overlapping bounds that main.py could save as selections.
    return checked_windows()


def plot_impulse_diagnostics(time_s, torque_dyn_cm, results):
    """
    Plot each transient and its locally accumulated RAW impulse.

    Top row
        Reconstructed torque near startup/shutdown.

    Bottom row
        Cumulative raw impulse, reset to zero at the left endpoint of each
        selected integration window.

    If the theoretical impulse is available, it is drawn as a horizontal
    reference line on the cumulative-impulse plot.
    """
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 6), constrained_layout=True)

    for col, row in enumerate(results):
        start_s = row["window_start_s"]
        end_s = row["window_end_s"]

        # Show extra context outside the integration window so missed pulse
        # tails are visually obvious.
        context = (
            max(time_s[0], start_s - 0.75),
            min(time_s[-1], end_s + 0.75),
        )

        context_time, context_torque = window_samples(
            time_s, torque_dyn_cm, context
        )
        axes[0, col].plot(context_time, context_torque, color="hotpink")
        axes[0, col].axvspan(
            start_s, end_s, color="tab:blue", alpha=0.15
        )
        axes[0, col].set(
            xlim=context,
            title=row["event"],
            ylabel="Torque (dyn*cm)",
        )

        local_time, local_torque = window_samples(
            time_s, torque_dyn_cm, (start_s, end_s)
        )
        cumulative_impulse = cumulative_trapezoid(
            local_torque, local_time, initial=0
        )

        axes[1, col].plot(
            local_time,
            cumulative_impulse,
            label=f"measured J = {cumulative_impulse[-1]:.4g}",
        )

        theory = row["J_theory_dyn_cm_s"]
        if theory is not None:
            axes[1, col].axhline(
                theory,
                linestyle="--",
                label=f"theory J = {theory:.4g}",
            )

        axes[1, col].axhline(0, color="gray", linewidth=0.5)
        axes[1, col].set(
            xlim=(start_s, end_s),
            xlabel="Experimental time (s)",
            ylabel="Cumulative torque (dyn*cm*s)",
        )
        axes[1, col].legend(fontsize=8)

    return fig
