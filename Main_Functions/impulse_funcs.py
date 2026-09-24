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
    event_time_s,
    chosen_window,
    other_window,
    changes_s=(0.0, 0.25, 0.5),
):
    """
    Generate nearby windows for checking sensitivity to integration limits.

    Each positive change expands the window equally on both sides.
    Each negative change contracts it equally on both sides.
    Invalid windows are discarded.
    """
    candidates = []

    for change in (-0.5, -0.25, 0.0, 0.25, 0.5):
        start_s = chosen_window[0] - change
        end_s = chosen_window[1] + change

        inside_data = time_s[0] <= start_s < end_s <= time_s[-1]
        contains_event = start_s <= event_time_s <= end_s
        avoids_other_event_window = (
            end_s <= other_window[0] or start_s >= other_window[1]
        )

        if inside_data and contains_event and avoids_other_event_window:
            candidates.append((start_s, end_s))

    return candidates


def analyze_impulses(
    time_s,
    torque_dyn_cm,
    full_time_s,
    *,
    pump_start_s,
    pump_stop_s,
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
        pump-event times are interpreted in the original experiment coordinates.
    pump_start_s, pump_stop_s
        Known pump startup and shutdown times.
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

    if (
        not np.all(np.isfinite([pump_start_s, pump_stop_s]))
        or not time_s[0] <= pump_start_s < pump_stop_s <= time_s[-1]
    ):
        raise ValueError(
            "Pump startup and shutdown times must be ordered and inside the data."
        )

    startup_window = validate_window(time_s, startup_window)
    shutdown_window = validate_window(time_s, shutdown_window)

    if startup_window[1] > shutdown_window[0]:
        raise ValueError("Startup and shutdown integration windows must not overlap.")

    if not startup_window[0] <= pump_start_s <= startup_window[1]:
        raise ValueError("The startup integration window must contain pump_start_s.")

    if not shutdown_window[0] <= pump_stop_s <= shutdown_window[1]:
        raise ValueError("The shutdown integration window must contain pump_stop_s.")

    theory_start, theory_stop, q_signed = theoretical_impulses(
        rho_g_cm3=rho_g_cm3,
        flow_rate_cm3_s=flow_rate_cm3_s,
        geometry_factor_cm2=geometry_factor_cm2,
        spin_dir=spin_dir,
        sign_convention=sign_convention,
    )

    event_specs = [
        ("startup", pump_start_s, startup_window, shutdown_window, theory_start),
        ("shutdown", pump_stop_s, shutdown_window, startup_window, theory_stop),
    ]

    results = []

    for event, event_time_s, window, other_window, theory in event_specs:
        measured = integrate_impulse(time_s, torque_dyn_cm, window)

        sensitivity_windows = _valid_sensitivity_windows(
            time_s,
            event_time_s,
            window,
            other_window,
        )
        sensitivity_values = [
            integrate_impulse(time_s, torque_dyn_cm, candidate)
            for candidate in sensitivity_windows
        ]

        residual = None if theory is None else measured - theory
        ratio = None if theory is None or theory == 0 else measured / theory

        results.append(
            {
                "event": event,
                "event_time_s": event_time_s,
                "window_start_s": window[0],
                "window_end_s": window[1],
                "J_measured_dyn_cm_s": measured,
                "J_theory_dyn_cm_s": theory,
                "residual_dyn_cm_s": residual,
                "measured_over_theory": ratio,
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


def report_impulses(results):
    """
    Print a concise experimental-versus-theory comparison.
    """
    print("\nRAW TRANSIENT IMPULSE COMPARISON")
    print("Units: dyn*cm*s = g*cm^2/s")
    print("Measured impulse = integral of reconstructed torque over the selected window.")
    print("No baseline or steady-flow correction is applied.")

    for row in results:
        print(
            f"\n{row['event'].upper()} "
            f"[{row['window_start_s']:g}, {row['window_end_s']:g}] s"
        )
        print(f"  measured J = {row['J_measured_dyn_cm_s']:.6g}")

        if row["J_theory_dyn_cm_s"] is None:
            print("  theory     = unavailable (flow rate and/or geometry factor missing)")
        else:
            print(f"  theory J   = {row['J_theory_dyn_cm_s']:.6g}")
            print(f"  residual   = {row['residual_dyn_cm_s']:.6g}")

            ratio = row["measured_over_theory"]
            if ratio is None:
                print("  measured/theory = undefined (zero theoretical impulse)")
            else:
                print(f"  measured/theory = {ratio:.6g}")

        print(
            "  window sensitivity = "
            f"[{row['sensitivity_min_dyn_cm_s']:.6g}, "
            f"{row['sensitivity_max_dyn_cm_s']:.6g}]"
        )

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
            relative_magnitude_mismatch = abs(abs(j_start) - abs(j_stop)) / larger
            print(
                "  relative magnitude mismatch = "
                f"{100 * relative_magnitude_mismatch:.3g}%"
            )


def save_impulse_summary(path, run_identity, results):
    """
    Save the raw impulse/theory comparison to a CSV file.
    """
    records = [dict(run_identity=run_identity, **row) for row in results]

    with open(path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def mark_impulse_windows(ax, results, annotate=True):
    """
    Mark pump events and selected integration windows on an existing axis.
    """
    colors = ("tab:blue", "tab:orange")

    for row, color in zip(results, colors):
        ax.axvline(
            row["event_time_s"],
            color=color,
            linestyle="--",
            label=f"Pump {row['event']}",
        )
        ax.axvspan(
            row["window_start_s"],
            row["window_end_s"],
            color=color,
            alpha=0.15,
            label=f"{row['event']} integration window",
        )

        if annotate:
            label = f"{row['event']} J = {row['J_measured_dyn_cm_s']:.4g} dyn*cm*s"
            if row["J_theory_dyn_cm_s"] is not None:
                label += f"\ntheory = {row['J_theory_dyn_cm_s']:.4g}"

            ax.text(
                0.02 if row["event"] == "startup" else 0.52,
                0.98,
                label,
                transform=ax.transAxes,
                va="top",
                fontsize=8,
                bbox=dict(facecolor="white", alpha=0.8, edgecolor=color),
            )

    ax.legend(fontsize=7, loc="lower right")


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
        axes[0, col].axvline(
            row["event_time_s"], color="tab:blue", linestyle="--"
        )
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

        axes[1, col].axvline(
            row["event_time_s"], color="gray", linestyle="--"
        )
        axes[1, col].axhline(0, color="gray", linewidth=0.5)
        axes[1, col].set(
            xlim=(start_s, end_s),
            xlabel="Experimental time (s)",
            ylabel="Cumulative impulse (dyn*cm*s)",
        )
        axes[1, col].legend(fontsize=8)

    return fig
