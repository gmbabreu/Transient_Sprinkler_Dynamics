import contextlib
import csv
import io
import tempfile
import unittest
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Main_Functions.impulse_funcs import (
    analyze_impulses, integrate_impulse, report_impulses,
    save_impulse_summary, mark_impulse_windows, plot_impulse_diagnostics,
)


class ImpulseTests(unittest.TestCase):
    def setUp(self):
        self.t = np.arange(0., 30.01, .125)
        # Triangular pulses with signed areas +6 and -6 on a constant background.
        start = 12 * np.maximum(0., 1 - abs(self.t - 10) / .5)
        stop = -12 * np.maximum(0., 1 - abs(self.t - 20) / .5)
        self.y = 3 + start + stop
        self.settings = dict(startup_window=(9.03, 12.07),
                             shutdown_window=(19.03, 22.07))

    def analyze(self, **overrides):
        return analyze_impulses(self.t, self.y, self.t, **(self.settings | overrides))

    def test_interpolated_endpoints_and_signed_areas(self):
        self.assertAlmostEqual(integrate_impulse([0, 1, 2, 3], [1, 3, 5, 7], (.25, 2.75)), 10.)
        for sign in (-1, 1):
            self.assertAlmostEqual(integrate_impulse([0, 1, 2], np.array([0, 4, 0])*sign,
                                                    (.25, 1.75)), sign*3.75)

    def test_windows_alone_determine_raw_impulses_and_export(self):
        rows = self.analyze()
        for row, kick in zip(rows, (6., -6.)):
            duration = row['window_end_s'] - row['window_start_s']
            self.assertAlmostEqual(row['J_measured_dyn_cm_s'], kick + 3*duration)
            self.assertIsNone(row['J_theory_dyn_cm_s'])
            self.assertNotIn('event_time_s', row)
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            report_impulses(rows)
        self.assertIn('theory     = unavailable', output.getvalue())
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'summary.csv'
            save_impulse_summary(path, 'synthetic', rows)
            with path.open() as stream:
                saved = list(csv.DictReader(stream))
            self.assertEqual(len(saved), 2)
            self.assertEqual(saved[0]['run_identity'], 'synthetic')
            self.assertEqual(saved[0]['J_theory_dyn_cm_s'], '')
            self.assertNotIn('event_time_s', saved[0])
        fig = plot_impulse_diagnostics(self.t, self.y, rows)
        mark_impulse_windows(fig.axes[0], rows)
        fig.canvas.draw()
        for ax, row in zip(fig.axes[2:], rows):
            self.assertEqual(ax.lines[0].get_ydata()[0], 0)
            self.assertAlmostEqual(ax.lines[0].get_ydata()[-1], row['J_measured_dyn_cm_s'])
        plt.close(fig)

    def test_window_markings_follow_plot_offset_without_changing_results(self):
        rows = self.analyze()
        original = [dict(row) for row in rows]
        fig, ax = plt.subplots()
        try:
            mark_impulse_windows(ax, rows, time_offset_s=9., trial_id='1000r1')
            self.assertEqual(len(ax.patches), 2)
            for span, row in zip(ax.patches, rows):
                # Polygon and Rectangle spans differ by Matplotlib version;
                # transform their vertices back to the displayed data axis.
                vertices = ax.transData.inverted().transform(
                    span.get_transform().transform(span.get_path().vertices))
                np.testing.assert_allclose(
                    [vertices[:, 0].min(), vertices[:, 0].max()],
                    [row['window_start_s'] - 9., row['window_end_s'] - 9.])
                self.assertIn('1000r1', span.get_label())
            self.assertEqual(rows, original)
        finally:
            plt.close(fig)

    def test_theory_signs_and_zero(self):
        for direction, flow_sign in [('f', 1), ('r', -1)]:
            for convention in (-1, 1):
                rows = self.analyze(flow_rate_cm3_s=2, geometry_factor_cm2=3,
                                    spin_dir=direction, sign_convention=convention)
                for row, event_sign in zip(rows, (-1, 1)):
                    theory = event_sign*flow_sign*convention*6
                    self.assertEqual(row['J_theory_dyn_cm_s'], theory)
                    self.assertEqual(row['residual_dyn_cm_s'], row['J_measured_dyn_cm_s']-theory)
                    self.assertEqual(row['measured_over_theory'], row['J_measured_dyn_cm_s']/theory)
        for row in self.analyze(flow_rate_cm3_s=2, geometry_factor_cm2=0):
            self.assertEqual(row['J_theory_dyn_cm_s'], 0)
            self.assertIsNone(row['measured_over_theory'])

    def test_invalid_inputs(self):
        for t in ([0, 0, 1], [0, 2, 1], [0, float('nan'), 2]):
            with self.assertRaises(ValueError):
                integrate_impulse(t, [0, 1, 0], (0, 1))
        for window in ((-1, 1), (1, 4), (1, 1)):
            with self.assertRaises(ValueError):
                integrate_impulse([0, 1, 2], [0, 1, 0], window)
        for override in [dict(startup_window=(9, 20)), dict(startup_window=(0, 31)),
                         dict(shutdown_window=(22, 21)), dict(sign_convention=0),
                         dict(flow_rate_cm3_s=-1)]:
            with self.assertRaises(ValueError):
                self.analyze(**override)
        with self.assertRaises(ValueError):
            analyze_impulses(self.t+.1, self.y, self.t, **self.settings)

    def test_windows_need_no_separate_event_time(self):
        rows = self.analyze(shutdown_window=(21, 23))
        self.assertAlmostEqual(rows[1]['J_measured_dyn_cm_s'], 6.)

    def test_sensitivity_reports_contraction_before_expansion(self):
        rows = analyze_impulses(self.t, np.full_like(self.t, -2.), self.t,
                                startup_window=(9., 12.), shutdown_window=(19., 22.))
        for row in rows:
            self.assertEqual(row['sensitivity_contracted_dyn_cm_s'], -5.)
            self.assertEqual(row['sensitivity_expanded_dyn_cm_s'], -7.)
            # The plotted range still uses numerical minimum/maximum.
            self.assertEqual(row['sensitivity_min_dyn_cm_s'], -7.)
            self.assertEqual(row['sensitivity_max_dyn_cm_s'], -5.)
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            report_impulses(rows)
        self.assertEqual(output.getvalue().count(
            'window sensitivity (contracted, expanded) = [-5, -7]'), 2)

    def test_sensitivity_filters_invalid_variants(self):
        rows = self.analyze()
        self.assertEqual(rows[0]['sensitivity_windows_s'], '9.28:11.82; 9.03:12.07; 8.78:12.32')
        rows = self.analyze(startup_window=(9.9, 19.0), shutdown_window=(19.1, 20.1))
        self.assertEqual(rows[0]['sensitivity_windows_s'], '10.15:18.75; 9.9:19')
        self.assertEqual(rows[1]['sensitivity_windows_s'], '19.35:19.85; 19.1:20.1')
        self.assertIsNone(rows[0]['sensitivity_expanded_dyn_cm_s'])
        rows = self.analyze(startup_window=(9.9, 10.2))
        self.assertIsNone(rows[0]['sensitivity_contracted_dyn_cm_s'])
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            report_impulses(rows)
        self.assertIn('window sensitivity (contracted, expanded) = [unavailable, ', output.getvalue())


if __name__ == '__main__':
    unittest.main()
