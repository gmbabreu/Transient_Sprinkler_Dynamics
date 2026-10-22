"""Run main's real reconstruction/analysis/export with deterministic CSV input.

Saved metadata replaces interactive selections. These integration checks run
with both NumPy 1.x and NumPy 2.x.
"""

import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import main
from Main_Functions import scientific_pipeline


class MainReconstructionTests(unittest.TestCase):
    def test_export_matches_original_time_and_impulse_analysis_torque(self):
        for count in (1530, 2048, 2049, 5001):
            with self.subTest(count=count), tempfile.TemporaryDirectory() as directory:
                directory = Path(directory)
                config = main.load_config()
                config['paths']['data'] = directory
                config['paths']['output'] = directory / 'results'
                config['trials'] = ['1000r1']
                time = 5. + np.arange(count) * .02
                angle = .025 * np.exp(-.03 * (time - 5.)) * np.cos(2 * (time - 5.))
                np.savetxt(directory / '1000r1_data.csv',
                           np.column_stack([time, np.rad2deg(angle)]),
                           delimiter=',', header='time_seconds,angle_degrees', comments='')
                metadata = dict(valid=True, fit_start_s=10.,
                                audio_candidates_s=dict(startup=10.149896, shutdown=20.149896),
                                integration_windows_s=dict(startup=[11.149896, 12.149896],
                                                           shutdown=[21.149896, 22.149896]))
                (directory / '1000r1_metadata.json').write_text(json.dumps(metadata))
                with patch.object(scientific_pipeline, 'analyze_impulses', wraps=scientific_pipeline.analyze_impulses) as analyze, \
                     patch.object(scientific_pipeline, 'phi_from_torque', wraps=scientific_pipeline.phi_from_torque) as forward, \
                     contextlib.redirect_stdout(io.StringIO()):
                    trials, failures = main.run_group(config, 1000, 'r', interactive=False, show=False)
                try:
                    self.assertFalse(failures)
                    trial = trials[0]
                    export = directory / 'results/1000r/1000r1_signal.csv'
                    self.assertEqual(export.read_text().splitlines()[0], 'time_seconds,torque_dyn_cm')
                    saved = np.loadtxt(export, delimiter=',', skiprows=1)
                    self.assertEqual(saved.shape, (count, 2))
                    analyze.assert_called_once()
                    analyzed_time, analyzed_torque, original_time = analyze.call_args.args
                    np.testing.assert_array_equal(saved[:, 0], time)
                    np.testing.assert_array_equal(saved[:, 0], original_time)
                    np.testing.assert_array_equal(saved[:, 0], analyzed_time)
                    np.testing.assert_array_equal(saved[:, 1], analyzed_torque)
                    np.testing.assert_array_equal(saved[:, 1], trial['torque'])
                    # Retain the original check that export uses dimensional
                    # torque, rather than the unscaled forcing sent to the ODE.
                    forward.assert_called_once()
                    unscaled_forcing = forward.call_args.args[2]
                    np.testing.assert_array_equal(saved[:, 1], unscaled_forcing * trial['inertia'])
                    self.assertFalse(np.allclose(saved[:, 1], unscaled_forcing))
                    self.assertEqual(len(trial['phi_forward']), count)
                    self.assertTrue(np.isfinite(trial['phi_forward']).all())
                finally:
                    plt.close('all')


if __name__ == '__main__':
    unittest.main()
