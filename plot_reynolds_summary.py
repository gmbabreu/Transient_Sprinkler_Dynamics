"""Plot configured, paired event integrals from existing group summaries only.

Regenerate summaries with main.py after changing analysis inputs or selections.
This command never reconstructs signals or updates group outputs/metadata.
"""

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from main import load_config, parse_trial


EVENTS = ("startup", "shutdown")
FIGURES = {
    "startup": ("startup_vs_reynolds.png", "Startup integral"),
    "shutdown": ("shutdown_vs_reynolds.png", "Shutdown integral"),
    "imbalance": ("mismatch_vs_reynolds.png", "Startup + shutdown imbalance"),
}
COMBINED_FIGURES = {
    "signed": "startup_shutdown_signed_vs_reynolds.png",
    "magnitude": "startup_shutdown_magnitude_vs_reynolds.png",
}
OUTPUT_FILENAMES = tuple(filename for filename, _ in FIGURES.values()) + tuple(COMBINED_FIGURES.values())
# Compare stored conventions, allowing only floating-point roundoff. Flow is
# normalized by Re because the configured flow conversion is shared across Re.
PHYSICS = ("rho_g_cm3", "geometry_factor_cm2", "sign_convention",
           "flow_rate_per_re_cm3_s", "Q_signed_per_re_cm3_s")


@dataclass
class TrialPair:
    trial_id: str
    measured: dict
    theory: dict
    conventions: tuple


@dataclass
class GroupSummary:
    reynolds: int
    trials: list
    theory: dict

    def values(self, event):
        return np.array([trial.measured[event] for trial in self.trials])


class IncompatibleSettings(ValueError):
    """Stored physical conventions cannot be combined in one experiment."""


def finite_number(value):
    number = float(value)
    if not np.isfinite(number):
        raise ValueError(f"nonfinite value {value!r}")
    return number


def optional_number(value):
    return None if value is None or value == "" else finite_number(value)


def agrees(first, second):
    if first is None or second is None:
        return first is second
    return bool(np.isclose(first, second, rtol=1e-10, atol=1e-12))


def stored_conventions(row, reynolds, direction):
    """Validate stored settings without recalculating any theoretical impulse."""
    if row.get("spin_dir") != direction:
        raise ValueError("spin_dir does not match the requested direction/identity")
    rho = finite_number(row["rho_g_cm3"])
    sign = finite_number(row["sign_convention"])
    geometry = optional_number(row["geometry_factor_cm2"])
    flow = optional_number(row["flow_rate_cm3_s"])
    signed_flow = optional_number(row["Q_signed_cm3_s"])
    if rho <= 0 or sign not in (-1, 1) or (flow is not None and flow <= 0):
        raise ValueError("invalid density, sign convention or total-flow magnitude")
    expected_flow = None if flow is None else flow * (1 if direction == "f" else -1)
    if not agrees(signed_flow, expected_flow):
        raise ValueError("Q_signed_cm3_s disagrees with flow magnitude/direction")
    return (rho, geometry, sign, None if flow is None else flow / reynolds,
            None if signed_flow is None else signed_flow / reynolds)


def group_prediction(trials, event, group, report):
    available = [trial.theory[event] for trial in trials if trial.theory[event] is not None]
    if not available:
        report(f"{group}: {event} theory unavailable; omitting theory marker.")
        return None
    if not all(agrees(value, available[0]) for value in available):
        report(f"{group}: inconsistent {event} theory predictions {available}; "
               "omitting theory marker (no averaging).")
        return None
    return available[0]


def read_group(path, selected, reynolds, direction, report):
    """Read once, then pair by exact identity and event, never CSV row order."""
    group = f"{reynolds}{direction}"
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream, strict=True)
            required = {"run_identity", "event", "J_measured_dyn_cm_s", "spin_dir",
                        "rho_g_cm3", "geometry_factor_cm2", "sign_convention",
                        "flow_rate_cm3_s", "Q_signed_cm3_s"}
            if not required.issubset(reader.fieldnames or []):
                report(f"{group}: invalid summary; missing columns "
                       f"{sorted(required - set(reader.fieldnames or []))}; skipped.")
                return None
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        report(f"{group}: missing/unreadable summary {path}: {exc}; skipped.")
        return None

    by_trial = {trial_id: [] for trial_id in selected}
    for row in rows:
        trial_id = row.get("run_identity")
        try:
            identity = parse_trial(trial_id or "")
        except ValueError as exc:
            report(f"{group}: excluded malformed identity: {exc}")
            continue
        if identity != (reynolds, direction):
            report(f"{group}: excluded {trial_id}: wrong group/direction.")
        elif trial_id not in by_trial:
            report(f"{group}: excluded deselected trial {trial_id}.")
        else:
            by_trial[trial_id].append(row)

    trials = []
    for trial_id, trial_rows in by_trial.items():
        events = [row.get("event") for row in trial_rows]
        if len(events) != 2 or any(events.count(event) != 1 for event in EVENTS):
            report(f"{group}: excluded {trial_id}: duplicate, invalid or incomplete "
                   f"event pair {events}; require exactly one startup and one shutdown.")
            continue
        try:
            measured, theory, conventions = {}, {}, []
            for row in trial_rows:
                if None in row or any(value is None for value in row.values()):
                    raise ValueError("malformed CSV row")
                event = row["event"]
                measured[event] = finite_number(row["J_measured_dyn_cm_s"])
                conventions.append(stored_conventions(row, reynolds, direction))
                try:
                    theory[event] = optional_number(row.get("J_theory_dyn_cm_s"))
                except (ValueError, TypeError):
                    theory[event] = None
                if theory[event] is None:
                    report(f"{trial_id}: {event} theory unavailable/invalid; keeping measurement.")
            if not all(agrees(a, b) for a, b in zip(*conventions)):
                raise ValueError("incompatible stored physical settings between events")
            measured["imbalance"] = finite_number(measured["startup"] + measured["shutdown"])
            theory["imbalance"] = (theory["startup"] + theory["shutdown"]
                                   if all(theory[e] is not None for e in EVENTS) else None)
            if theory["imbalance"] is not None and not np.isfinite(theory["imbalance"]):
                theory["imbalance"] = None
            trials.append(TrialPair(trial_id, measured, theory, conventions[0]))
        except (ValueError, TypeError, KeyError) as exc:
            report(f"{group}: excluded {trial_id}: {exc}.")
    if not trials:
        report(f"{group}: no valid paired trials; skipped.")
        return None
    theory = {event: group_prediction(trials, event, group, report) for event in EVENTS}
    # Even cancelling disagreements in the two event predictions must not
    # manufacture a group imbalance prediction.
    theory["imbalance"] = (group_prediction(trials, "imbalance", group, report)
                           if all(theory[e] is not None for e in EVENTS) else None)
    return GroupSummary(reynolds, trials, theory)


def collect_summaries(config, direction, report=print):
    if direction not in ("f", "r"):
        raise ValueError("direction must be f or r")
    selected = {}
    for trial_id in config["trials"]:
        reynolds, trial_direction = parse_trial(trial_id)
        if trial_direction == direction:
            selected.setdefault(reynolds, []).append(trial_id)
    groups = []
    for reynolds, trial_ids in sorted(selected.items()):
        path = config["paths"]["output"] / f"{reynolds}{direction}" / "jump_summary.csv"
        group = read_group(path, trial_ids, reynolds, direction, report)
        if group is not None:
            groups.append(group)
    settings = [(trial.trial_id, trial.conventions) for group in groups for trial in group.trials]
    if settings:
        conflicts = []
        for index, field in enumerate(PHYSICS):
            if any(not agrees(values[index], settings[0][1][index]) for _, values in settings):
                conflicts.append(f"{field}: " + ", ".join(
                    f"{trial_id}={values[index]}" for trial_id, values in settings))
        if conflicts:
            raise IncompatibleSettings("Incompatible stored physical settings; cannot combine "
                                       "these trials/groups:\n" + "\n".join(conflicts))
    return groups


def mean_annotation(group, event, mean):
    """Label means with event error or the existing measured magnitude mismatch."""
    if event == "imbalance":
        start, stop = (np.mean(group.values(name)) for name in EVENTS)
        larger = max(abs(start), abs(stop))
        mismatch = None if larger == 0 else 100 * abs(abs(start) - abs(stop)) / larger
        return "Mismatch N/A" if mismatch is None else f"Mismatch {mismatch:.2f}%"
    theory = group.theory[event]
    error = None if theory is None or theory == 0 else 100 * abs(mean - theory) / abs(theory)
    return "Error N/A" if error is None else f"Error {error:.2f}%"


def plot_summary(groups, direction):
    """Return trial/mean/theory figures, annotating only the means."""
    figures = {}
    groups = sorted(groups, key=lambda group: group.reynolds)
    for event, (filename, title) in FIGURES.items():
        fig, ax = plt.subplots(figsize=(7, 4.5))
        x = [group.reynolds for group in groups]
        ax.scatter([group.reynolds for group in groups for _ in group.trials],
                   [value for group in groups for value in group.values(event)],
                   s=18, alpha=.5, color="tab:blue", label="Trials", zorder=3)
        means = [np.mean(group.values(event)) for group in groups]
        ax.scatter(x, means, s=55, marker="o", color="tab:orange",
                   alpha=.65, linewidths=0, label="Mean", zorder=4)
        predicted = [group for group in groups if group.theory[event] is not None]
        if predicted:
            ax.scatter([group.reynolds for group in predicted],
                       [group.theory[event] for group in predicted],
                       s=85, marker="_", color="black", alpha=.7, linewidths=2,
                       label="Theory", zorder=5)
        midpoint = (x[0] + x[-1]) / 2
        for group, mean in zip(groups, means):
            label = mean_annotation(group, event, mean)
            right = group.reynolds <= midpoint
            ax.annotate(label, (group.reynolds, mean),
                        xytext=(12 if right else -12, 0), textcoords="offset points",
                        ha="left" if right else "right", va="center", fontsize=9,
                        bbox=dict(facecolor="white", edgecolor="none", alpha=.85), zorder=6)
        ax.axhline(0, color="0.4", linewidth=.9, zorder=0)
        flow_label = "Reverse" if direction == "r" else "Forward"
        ylabel = (f"{title}\n(dyn·cm·s)" if event == "imbalance"
                  else f"{title} (dyn·cm·s)")
        ax.set(xlabel="Reynolds number", ylabel=ylabel,
               title=f"{title} vs Reynolds number\n{flow_label} flow")
        ax.set_xticks(x)
        ax.ticklabel_format(axis="x", style="plain", useOffset=False)
        ax.grid(alpha=.2)
        ax.set_axisbelow(True)
        ax.margins(x=.1, y=.15)
        ax.legend(frameon=False, ncol=3)
        fig.tight_layout()
        figures[filename] = fig
    return figures


def combined_theory_slopes(groups, report=print):
    """Derive origin-line slopes from stored predictions, never measured data.

    Missing event predictions permit a line from the remaining consistent
    predictions. Conflicts suppress the affected event line; a shared magnitude
    line additionally requires no conflicting predictions in either event.
    """
    slopes = {}
    inconsistent = False
    for event in EVENTS:
        available = [(group.reynolds, group.theory[event] / group.reynolds)
                     for group in groups if group.theory[event] is not None]
        # A rejected group prediction is distinguishable from missing theory:
        # it has available trial predictions, but they failed group validation.
        rejected = [group.reynolds for group in groups if group.theory[event] is None
                    and any(trial.theory[event] is not None for trial in group.trials)]
        if rejected or (available and any(not agrees(value, available[0][1])
                                         for _, value in available)):
            report(f"Combined plots: inconsistent {event} theory slopes J_theory/Re "
                   f"{available}; rejected group predictions at Re={rejected}. "
                   "Omitting this theoretical line.")
            slopes[event] = None
            inconsistent = True
        elif available:
            slopes[event] = available[0][1]
        else:
            report(f"Combined plots: {event} theory unavailable; omitting this theoretical line.")
            slopes[event] = None
    asymmetric = [group.reynolds for group in groups
                  if all(group.theory[event] is not None for event in EVENTS)
                  and not agrees(group.theory['startup'] / group.reynolds,
                                 -group.theory['shutdown'] / group.reynolds)]
    if asymmetric or (all(slopes[event] is not None for event in EVENTS)
                      and not agrees(slopes['startup'], -slopes['shutdown'])):
        report(f"Combined plots: startup/shutdown theory slopes are not equal and opposite "
               f"(Re={asymmetric}, slopes={slopes}); omitting theoretical lines.")
        slopes = dict.fromkeys(EVENTS)
        inconsistent = True
    available = [slopes[event] for event in EVENTS if slopes[event] is not None]
    slopes['magnitude'] = abs(available[0]) if available and not inconsistent else None
    if inconsistent:
        report("Combined magnitude plot: inconsistent predictions; omitting shared theory line.")
    return slopes


def plot_combined_summary(groups, direction, report=print):
    """Add signed and mean-of-magnitudes views without changing stored values."""
    groups = sorted(groups, key=lambda group: group.reynolds)
    slopes = combined_theory_slopes(groups, report)
    figures = {}
    reynolds = [group.reynolds for group in groups]
    line_x = np.array([0., max(reynolds) * 1.1])
    for mode, filename in COMBINED_FIGURES.items():
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for event, color, marker in (("startup", "tab:blue", "o"),
                                     ("shutdown", "tab:orange", "o")):
            values = [group.values(event) for group in groups]
            if mode == 'magnitude':
                values = [np.abs(value) for value in values]
            style = dict(facecolors=color, edgecolors=color)
            ax.scatter([group.reynolds for group in groups for _ in group.trials],
                       np.concatenate(values), s=18, alpha=.5, marker=marker,
                       label=f"{event.capitalize()} trials", zorder=3, **style)
            ax.scatter(reynolds, [np.mean(value) for value in values],
                       s=55, alpha=.65, marker=marker, linewidths=1.3,
                       label=f"{event.capitalize()} mean", zorder=4, **style)
        if mode == 'signed':
            for event, color in (("startup", "tab:blue"), ("shutdown", "tab:orange")):
                if slopes[event] is not None:
                    ax.plot(line_x, slopes[event] * line_x, color=color,
                            linestyle='-', linewidth=1.4,
                            label=f"{event.capitalize()} theory", zorder=2)
        elif slopes['magnitude'] is not None:
            ax.plot(line_x, slopes['magnitude'] * line_x, color='black',
                    linestyle='-', linewidth=1.4, label='Theory magnitude', zorder=2)
        ax.axhline(0, color='0.4', linewidth=.9, zorder=0)
        flow = 'Reverse' if direction == 'r' else 'Forward'
        ax.set(xlabel='Reynolds number',
               ylabel=('Signed impulse' if mode == 'signed' else 'Impulse magnitude') + ' (dyn·cm·s)',
               title=f"Startup and shutdown: {mode} impulses\n{flow} flow",
               xlim=(0, line_x[-1]))
        ax.set_xticks([0] + reynolds)
        ax.ticklabel_format(axis='x', style='plain', useOffset=False)
        ax.grid(alpha=.2)
        ax.set_axisbelow(True)
        ax.margins(y=.15)
        if mode == 'magnitude':
            ax.set_ylim(bottom=0)
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        order = [f"{event} {series}" for event in ("Startup", "Shutdown")
                 for series in ("trials", "mean", "theory")]
        if mode == 'magnitude':
            order.insert(2, 'Theory magnitude')
        order = [label for label in order if label in by_label]
        ax.legend([by_label[label] for label in order], order, frameon=False, ncol=2, fontsize=8)
        fig.tight_layout()
        figures[filename] = fig
    return figures


def run_summary(config, direction, *, show=True, report=print):
    output = config["paths"]["output"] / "global" / direction
    try:
        groups = collect_summaries(config, direction, report)
    except IncompatibleSettings as exc:
        report(str(exc))
        groups = []
    if not groups:
        for filename in OUTPUT_FILENAMES:
            (output / filename).unlink(missing_ok=True)
        report(f"No comparable valid paired data for {direction}; removed old global figures.")
        return []
    figures = plot_summary(groups, direction)
    figures.update(plot_combined_summary(groups, direction, report))
    output.mkdir(parents=True, exist_ok=True)
    for filename, fig in figures.items():
        fig.savefig(output / filename, dpi=300)
    for group in groups:
        report(f"{group.reynolds}{direction}: {len(group.trials)} trials; means (dyn·cm·s): "
               + ", ".join(f"{event}={np.mean(group.values(event)):.10f}" for event in FIGURES))
    report(f"Saved {len(figures)} global figures to {output}")
    if show:
        plt.show()
    for fig in figures.values():
        plt.close(fig)
    return groups


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path,
                        default=Path(__file__).with_name("experiment_config.json"))
    parser.add_argument("--direction", required=True, choices=("f", "r"))
    parser.add_argument("--no-show", action="store_true", help="Save figures without displaying them")
    args = parser.parse_args()
    if args.no_show:
        plt.switch_backend("Agg")
    config = load_config(args.config)
    groups = run_summary(config, args.direction,
                         show=config["analysis"]["display_plots"] and not args.no_show)
    return 0 if groups else 1


if __name__ == "__main__":
    raise SystemExit(main())
