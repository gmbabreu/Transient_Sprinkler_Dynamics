"""Synthetic summary-only checks: selection, pairing, values and plot artists."""

import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection
import numpy as np

import main
import plot_reynolds_summary as summary


class ReynoldsSummaryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.addCleanup(lambda: plt.close("all"))
        self.root = Path(self.temp.name)
        self.config = main.load_config()
        self.config["paths"]["output"] = self.root / "results"
        self.config["trials"] = ["1500r1", "500r1", "1000r1", "1000r2", "1000f1"]
        self.messages = []

    def pair(self, trial_id, start=-10., stop=14., theory=(-8., 8.)):
        reynolds, direction = main.parse_trial(trial_id)
        flow = .01 * reynolds
        return [dict(run_identity=trial_id, event=event, J_measured_dyn_cm_s=measured,
                     J_theory_dyn_cm_s=predicted, rho_g_cm3=1., geometry_factor_cm2=-2.,
                     sign_convention=1, flow_rate_cm3_s=flow,
                     Q_signed_cm3_s=flow * (1 if direction == "f" else -1), spin_dir=direction)
                for event, measured, predicted in zip(summary.EVENTS, (start, stop), theory)]

    def write(self, group, rows):
        path = self.config["paths"]["output"] / group / "jump_summary.csv"
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        return path

    def collect(self, direction="r"):
        return summary.collect_summaries(self.config, direction, self.messages.append)

    def test_numeric_order_pairing_exact_series_and_no_uncertainty_artists(self):
        self.write("500r", self.pair("500r1", -5, 7, (-4, 4))[::-1])
        first = self.pair("1000r1", -10, 14)
        second = self.pair("1000r2", -20, 30)
        # Events interleaved and reversed: adjacent rows are not a trial pair.
        self.write("1000r", [second[1], first[0], second[0], first[1]])
        self.write("1500r", self.pair("1500r1", -30, 50, (-12, 12)))
        groups = self.collect()
        self.assertEqual([g.reynolds for g in groups], [500, 1000, 1500])
        self.assertEqual([t.trial_id for t in groups[1].trials], ["1000r1", "1000r2"])
        np.testing.assert_array_equal(groups[1].values("imbalance"), [4, 10])
        figures = summary.plot_summary(groups[::-1], "r")
        expected = {
            "startup_vs_reynolds.png": ([-5, -10, -20, -30], [-5, -15, -30], [-4, -8, -12]),
            "shutdown_vs_reynolds.png": ([7, 14, 30, 50], [7, 22, 50], [4, 8, 12]),
            "mismatch_vs_reynolds.png": ([2, 4, 10, 20], [2, 7, 20], [0, 0, 0]),
        }
        for filename, (trials, means, theory) in expected.items():
            ax = figures[filename].axes[0]
            artists = {artist.get_label(): artist for artist in ax.collections}
            self.assertEqual(set(artists), {"Trials", "Mean", "Theory"})
            for label, x, y in (("Trials", [500, 1000, 1000, 1500], trials),
                                ("Mean", [500, 1000, 1500], means),
                                ("Theory", [500, 1000, 1500], theory)):
                np.testing.assert_array_equal(artists[label].get_offsets(), np.column_stack([x, y]))
            for artist in artists.values():
                self.assertTrue(0 < artist.get_alpha() < 1)
            self.assertGreater(artists["Mean"].get_sizes()[0], artists["Trials"].get_sizes()[0])
            self.assertEqual(artists["Mean"].get_facecolors().shape, (1, 4))
            self.assertFalse(np.array_equal(artists["Mean"].get_facecolors()[0, :3],
                                           artists["Trials"].get_facecolors()[0, :3]))
            np.testing.assert_array_equal(ax.get_xticks(), [500, 1000, 1500])
            self.assertEqual([t.get_text() for t in ax.get_xticklabels()], ["500", "1000", "1500"])
            self.assertIn("dyn·cm·s", ax.get_ylabel())
            self.assertEqual([t.get_text() for t in ax.get_legend().get_texts()], ["Trials", "Mean", "Theory"])
            # Exactly the three scatter series and a zero line: no error bars,
            # line segments, bands, regressions, trial labels or extra axes.
            self.assertTrue(all(isinstance(a, PathCollection) for a in ax.collections))
            self.assertEqual(len(ax.lines), 1)
            np.testing.assert_array_equal(ax.lines[0].get_ydata(), [0, 0])
            self.assertEqual(len(ax.containers), 0)
            self.assertEqual(len(ax.patches), 0)
            self.assertEqual(len(ax.texts), len(means))
            for index, (annotation, x, mean, prediction) in enumerate(zip(ax.texts, [500, 1000, 1500], means, theory)):
                self.assertEqual(annotation.xy, (x, mean))
                if filename == "mismatch_vs_reynolds.png":
                    expected_error = ["Mismatch 28.57%", "Mismatch 31.82%", "Mismatch 40.00%"][index]
                else:
                    expected_error = f"Error {100 * abs(mean - prediction) / abs(prediction):.2f}%"
                self.assertEqual(annotation.get_text(), expected_error)
            self.assertEqual(len(figures[filename].axes), 1)
        self.assertIn("Startup + shutdown imbalance", figures["mismatch_vs_reynolds.png"].axes[0].get_title())

    def test_error_is_for_the_mean_not_average_of_trial_errors(self):
        self.config["trials"] = ["1000r1", "1000r2"]
        self.write("1000r", self.pair("1000r1", start=-4) + self.pair("1000r2", start=-12))
        groups = self.collect()
        figures = summary.plot_summary(groups, "r")
        annotations = figures["startup_vs_reynolds.png"].axes[0].texts
        self.assertEqual(len(annotations), 1)
        self.assertEqual(annotations[0].xy, (1000, -8))
        self.assertEqual(annotations[0].get_text(), "Error 0.00%")
        # Missing and inconsistent theory both yield N/A, without inventing a
        # percentage for event errors.
        groups[0].theory["startup"] = None
        figures = summary.plot_summary(groups, "r")
        self.assertEqual(figures["startup_vs_reynolds.png"].axes[0].texts[0].get_text(), "Error N/A")

    def test_imbalance_label_uses_mismatch_of_means_independent_of_theory(self):
        self.config["trials"] = ["1000r1", "1000r2"]
        self.write("1000r", self.pair("1000r1", -10, 20, (None, None))
                   + self.pair("1000r2", -30, 30, (None, None)))
        group, = self.collect()
        figures = summary.plot_summary([group], "r")
        annotation, = figures["mismatch_vs_reynolds.png"].axes[0].texts
        self.assertEqual(annotation.xy, (1000, 5))
        # Mean start=-20 and stop=25: 20%, not the 25% average trial mismatch.
        self.assertEqual(annotation.get_text(), "Mismatch 20.00%")
        for start, stop, expected in ((0, 0, "Mismatch N/A"),
                                      (-10, 10, "Mismatch 0.00%"),
                                      (10, 10, "Mismatch 0.00%"),
                                      (0, 10, "Mismatch 100.00%")):
            for trial in group.trials:
                trial.measured.update(startup=start, shutdown=stop, imbalance=start + stop)
            self.assertEqual(summary.mean_annotation(group, "imbalance", start + stop), expected)

    def test_invalid_duplicate_incomplete_and_deselected_pairs_are_excluded(self):
        self.config["trials"] = [f"1000r{i}" for i in range(1, 8)]
        rows = self.pair("1000r1")
        rows += self.pair("1000r2") + [self.pair("1000r2")[0]]  # Duplicate startup
        rows += [self.pair("1000r3")[1]]  # Missing startup
        bad = self.pair("1000r4")
        bad[0]["J_measured_dyn_cm_s"] = "nan"
        rows += bad
        bad = self.pair("1000r5")
        bad[1]["event"] = "steady"
        rows += bad
        bad = self.pair("1000r6")
        bad[0]["spin_dir"] = "f"
        rows += bad
        bad = self.pair("1000r7")
        bad[1]["J_measured_dyn_cm_s"] = "bad"
        rows += bad
        rows += self.pair("1000r8")  # No longer selected
        rows += self.pair("1000f1") + self.pair("500r1")  # Wrong direction/group
        malformed = self.pair("1000r1")
        for row in malformed:
            row["run_identity"] = "garbage"
        rows += malformed
        self.write("1000r", rows)
        groups = self.collect()
        self.assertEqual([t.trial_id for t in groups[0].trials], ["1000r1"])
        np.testing.assert_array_equal(groups[0].values("startup"), [-10])
        for fragment in ("1000r2", "1000r3", "1000r4", "1000r5", "1000r6", "1000r7",
                         "deselected trial 1000r8", "wrong group/direction", "malformed identity"):
            self.assertTrue(any(fragment in message for message in self.messages), fragment)

    def test_missing_groups_and_partial_or_invalid_theory_keep_measurements(self):
        first = self.pair("1000r1", theory=(None, 8.))
        second = self.pair("1000r2", -20, 30, (-8., "inf"))
        self.write("1000r", first + second)
        groups = self.collect()
        self.assertEqual([g.reynolds for g in groups], [1000])
        self.assertEqual(groups[0].theory, dict(startup=-8., shutdown=8., imbalance=None))
        figures = summary.plot_summary(groups, "r")
        ax = figures["mismatch_vs_reynolds.png"].axes[0]
        self.assertEqual([a.get_label() for a in ax.collections], ["Trials", "Mean"])
        np.testing.assert_array_equal(ax.collections[0].get_offsets()[:, 1], [4, 10])
        self.assertTrue(any("500r: missing/unreadable" in m for m in self.messages))
        self.assertTrue(any("1500r: missing/unreadable" in m for m in self.messages))
        self.assertTrue(any("1000r1: startup theory unavailable" in m for m in self.messages))
        self.assertTrue(any("imbalance theory unavailable" in m for m in self.messages))

    def test_absent_theory_column_is_supported(self):
        self.config["trials"] = ["1000r1"]
        rows = self.pair("1000r1")
        for row in rows:
            del row["J_theory_dyn_cm_s"]
        self.write("1000r", rows)
        group, = self.collect()
        self.assertEqual(group.theory, dict(startup=None, shutdown=None, imbalance=None))
        np.testing.assert_array_equal(group.values("shutdown"), [14])

    def test_inconsistent_theory_is_not_averaged_even_if_sums_cancel(self):
        self.config["trials"] = ["1000r1", "1000r2"]
        self.write("1000r", self.pair("1000r1") + self.pair("1000r2", theory=(-9, 9)))
        group, = self.collect()
        self.assertTrue(all(value is None for value in group.theory.values()))
        self.assertTrue(any("inconsistent startup theory" in m for m in self.messages))
        self.assertTrue(any("inconsistent shutdown theory" in m for m in self.messages))
        self.assertEqual(len(group.trials), 2)

    def test_direction_separation_and_single_read_per_selected_group(self):
        self.write("1000r", self.pair("1000r1") + self.pair("1000r2"))
        self.write("1000f", self.pair("1000f1", 10, -14, (8, -8)))
        original_open = Path.open
        opened = []

        def track_open(path, *args, **kwargs):
            opened.append(path)
            return original_open(path, *args, **kwargs)

        with patch.object(Path, "open", track_open):
            reverse = self.collect("r")
        self.assertEqual([p.parent.name for p in opened], ["500r", "1000r", "1500r"])
        self.assertTrue(all(p.name == "jump_summary.csv" for p in opened))
        forward, = self.collect("f")
        self.assertEqual([t.trial_id for t in forward.trials], ["1000f1"])
        np.testing.assert_array_equal(forward.values("startup"), [10])
        np.testing.assert_array_equal(reverse[0].values("startup"), [-10, -10])
        self.assertEqual(forward.theory["startup"], 8)

    def test_physical_conventions_must_agree_across_trials_and_groups(self):
        for same_group in (False, True):
            other = "1000r2" if same_group else "1500r1"
            for field, value in (("rho_g_cm3", 2.), ("geometry_factor_cm2", 3.),
                                 ("sign_convention", -1), ("flow_rate_cm3_s", 99.)):
                with self.subTest(same_group=same_group, field=field):
                    self.config["trials"] = ["1000r1", other]
                    rows = self.pair(other)
                    for row in rows:
                        row[field] = value
                        if field == "flow_rate_cm3_s":
                            row["Q_signed_cm3_s"] = -value
                    if same_group:
                        self.write("1000r", self.pair("1000r1") + rows)
                    else:
                        self.write("1000r", self.pair("1000r1"))
                        self.write("1500r", rows)
                    with self.assertRaisesRegex(summary.IncompatibleSettings, field.split("_cm3_s")[0]):
                        self.collect()

    def test_missing_optional_physics_and_roundoff_are_comparable(self):
        self.config["trials"] = ["1000r1", "1500r1"]
        for trial in self.config["trials"]:
            rows = self.pair(trial, theory=(None, None))
            for row in rows:
                row["geometry_factor_cm2"] = None
                if trial == "1500r1":
                    row["rho_g_cm3"] += 1e-12
            self.write(trial[:-1], rows)
        self.assertEqual(len(self.collect()), 2)

    def test_invalid_stored_settings_exclude_whole_pair(self):
        self.config["trials"] = ["1000r1", "1000r2"]
        for field, value in (("rho_g_cm3", 0), ("geometry_factor_cm2", "nan"),
                             ("sign_convention", 0), ("Q_signed_cm3_s", 10),
                             ("flow_rate_cm3_s", -10), ("rho_g_cm3", 2)):
            with self.subTest(field=field, value=value):
                bad = self.pair("1000r2")
                bad[0][field] = value
                self.write("1000r", self.pair("1000r1") + bad)
                group, = self.collect()
                self.assertEqual([t.trial_id for t in group.trials], ["1000r1"])

    def stale_outputs(self):
        output = self.config["paths"]["output"] / "global" / "r"
        output.mkdir(parents=True, exist_ok=True)
        for filename, _ in summary.FIGURES.values():
            (output / filename).write_bytes(b"old")
        (output / "keep.txt").write_text("unrelated")
        return output

    def test_no_valid_data_removes_only_own_old_figures(self):
        for rows in (None, self.pair("1000r1")[:1]):
            with self.subTest(rows=rows):
                output = self.stale_outputs()
                self.config["trials"] = ["1000r1"]
                if rows is not None:
                    self.write("1000r", rows)
                result = summary.run_summary(self.config, "r", show=False, report=self.messages.append)
                self.assertEqual(result, [])
                self.assertEqual([p.name for p in output.iterdir()], ["keep.txt"])

    def test_incompatible_groups_remove_old_figures_and_report_conflicts(self):
        output = self.stale_outputs()
        self.write("1000r", self.pair("1000r1") + self.pair("1000r2"))
        rows = self.pair("1500r1")
        for row in rows:
            row["sign_convention"] = -1
        self.write("1500r", rows)
        self.assertEqual(summary.run_summary(self.config, "r", show=False,
                                            report=self.messages.append), [])
        self.assertEqual([p.name for p in output.iterdir()], ["keep.txt"])
        self.assertTrue(any("Incompatible stored physical settings" in m for m in self.messages))

    def test_cli_config_paths_no_show_300dpi_and_preservation(self):
        path = self.write("1000r", self.pair("1000r1") + self.pair("1000r2"))
        # These unrelated outputs must neither be read nor updated.
        jump = path.with_name("1000r1_jumps.csv")
        jump.write_text("not an input")
        metadata = self.root / "1000r1_metadata.json"
        metadata.write_text("not an input")
        per_group = path.with_name("torque.png")
        per_group.write_bytes(b"keep existing figure")
        before = {p: p.read_bytes() for p in (path, jump, metadata, per_group)}
        raw_config = json.loads(Path(main.__file__).with_name("experiment_config.json").read_text())
        raw_config["trial_groups"] = ["1000r"]
        raw_config["trials_per_group"] = 2
        raw_config["paths"]["output"] = "results"
        config_path = self.root / "config.json"
        config_path.write_text(json.dumps(raw_config))
        with patch("sys.argv", ["plot_reynolds_summary.py", "--config", str(config_path),
                                "--direction", "r", "--no-show"]), \
             patch.object(plt, "show", side_effect=AssertionError("displayed")), \
             patch.object(main, "analyze_trial", side_effect=AssertionError("reanalyzed")), \
             patch.object(main, "save_metadata", side_effect=AssertionError("updated metadata")), \
             patch.object(plt.Figure, "savefig", autospec=True) as save, \
             contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(summary.main(), 0)
        self.assertEqual(save.call_count, 3)
        for call in save.call_args_list:
            self.assertEqual(call.kwargs["dpi"], 300)
            self.assertEqual(call.args[1].parent, self.root / "results" / "global" / "r")
        self.assertEqual({call.args[1].name for call in save.call_args_list},
                         {filename for filename, _ in summary.FIGURES.values()})
        self.assertEqual({p: p.read_bytes() for p in before}, before)
        self.assertEqual(plt.get_fignums(), [])


if __name__ == "__main__":
    unittest.main()
