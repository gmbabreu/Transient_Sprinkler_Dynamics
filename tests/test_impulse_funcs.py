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
    analyze_impulses, integrate_window, background_area, report_impulses,
    save_impulse_summary, plot_impulse_diagnostics,
)


class ImpulseTests(unittest.TestCase):
    def setUp(self):
        self.t = np.arange(0., 30.01, .125)
        # Offset 3, steady increment 4, triangular kicks of areas +6 and -6.
        start = 12 * np.maximum(0., 1 - abs(self.t-10)/.5)
        stop = -12 * np.maximum(0., 1 - abs(self.t-20)/.5)
        steady = 4 * ((self.t > 10) & (self.t < 20))
        steady[(self.t == 10) | (self.t == 20)] = 2
        self.y = 3 + steady + start + stop
        self.settings = dict(pump_start_s=10., pump_stop_s=20.,
                             startup_window=(9.03, 12.07), shutdown_window=(19.03, 22.07),
                             baseline_window=(5., 8.), steady_window=(14., 18.))

    def analyze(self, **overrides):
        return analyze_impulses(self.t, self.y, self.t, **(self.settings | overrides))

    def test_interpolated_endpoints_and_signed_areas(self):
        # Linear function: exact integral even when both endpoints fall between samples.
        self.assertAlmostEqual(integrate_window([0, 1, 2, 3], [1, 3, 5, 7], (.25, 2.75)), 10.)
        for sign in (-1, 1):
            self.assertAlmostEqual(integrate_window([0, 1, 2], np.array([0, 4, 0])*sign, (.25, 1.75)), sign*3.75)

    def test_backgrounds_missing_theory_and_export(self):
        rows = self.analyze()
        for row, kick in zip(rows, (6., -6.)):
            a, b = row['window_start_s'], row['window_end_s']
            self.assertAlmostEqual(row['tau_off_dyn_cm'], 3.)
            self.assertAlmostEqual(row['tau_on_dyn_cm'], 7.)
            self.assertAlmostEqual(row['J_raw_dyn_cm_s'], kick+background_area(a, b, 3, 7, 10, 20))
            self.assertAlmostEqual(row['J_offset_dyn_cm_s'], row['J_raw_dyn_cm_s']-3*(b-a))
            self.assertAlmostEqual(row['J_step_dyn_cm_s'], kick)
            self.assertIsNone(row['theory_dyn_cm_s'])
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            report_impulses(rows)
        self.assertIn('theory unavailable', output.getvalue())
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'summary.csv'
            save_impulse_summary(path, 'synthetic', rows)
            with path.open() as f:
                saved = list(csv.DictReader(f))
            self.assertEqual(len(saved), 2)
            self.assertEqual(saved[0]['run_identity'], 'synthetic')
            self.assertEqual(saved[0]['theory_dyn_cm_s'], '')
        fig = plot_impulse_diagnostics(self.t, self.y, rows)
        fig.canvas.draw()
        for ax in fig.axes[2:]:
            for line in ax.lines[:3]:
                self.assertEqual(line.get_ydata()[0], 0)
        plt.close(fig)

    def test_theory_signs_and_zero(self):
        for direction, flow_sign in [('f', 1), ('r', -1)]:
            for convention in (-1, 1):
                rows = self.analyze(flow_rate_cm3_s=2, geometry_factor_cm2=3,
                                    spin_dir=direction, sign_convention=convention)
                for row, event_sign in zip(rows, (-1, 1)):
                    theory = event_sign*flow_sign*convention*6
                    self.assertEqual(row['theory_dyn_cm_s'], theory)
                    for k in ('raw', 'offset', 'step'):
                        self.assertEqual(row[f'residual_{k}_dyn_cm_s'], row[f'J_{k}_dyn_cm_s']-theory)
                        self.assertEqual(row[f'ratio_{k}'], row[f'J_{k}_dyn_cm_s']/theory)
        for row in self.analyze(flow_rate_cm3_s=2, geometry_factor_cm2=0):
            self.assertEqual(row['theory_dyn_cm_s'], 0)
            self.assertIsNone(row['ratio_raw'])

    def test_invalid_inputs(self):
        for t in ([0, 0, 1], [0, 2, 1], [0, float('nan'), 2]):
            with self.assertRaises(ValueError):
                integrate_window(t, [0, 1, 0], (0, 1))
        for w in ((-1, 1), (1, 4), (1, 1)):
            with self.assertRaises(ValueError):
                integrate_window([0, 1, 2], [0, 1, 0], w)
        for override in [dict(startup_window=(9, 20)), dict(startup_window=(0, 31)),
                         dict(shutdown_window=(21, 23)), dict(sign_convention=0),
                         dict(flow_rate_cm3_s=-1)]:
            with self.assertRaises(ValueError):
                self.analyze(**override)
        with self.assertRaises(ValueError):
            analyze_impulses(self.t+.1, self.y, self.t, **self.settings)

    def test_sensitivity_filters_invalid_variants(self):
        rows = self.analyze(startup_window=(9.9, 19.0), shutdown_window=(19.1, 20.1))
        self.assertEqual(rows[0]['sensitivity_windows_s'], '9.9:19')
        self.assertEqual(rows[1]['sensitivity_windows_s'], '19.1:20.1')


if __name__ == '__main__':
    unittest.main()
