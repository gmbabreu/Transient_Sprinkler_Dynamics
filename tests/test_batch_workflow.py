"""Focused checks for group selection, timing reuse and invalid-trial exclusion."""

import contextlib
import csv
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


class BatchConfigTests(unittest.TestCase):
    """Both workflow stages need the same explicit group/repeat expansion."""

    def load(self, groups, count):
        config = json.loads(Path(main.__file__).with_name('experiment_config.json').read_text())
        config['trial_groups'] = groups
        config['trials_per_group'] = count
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'config.json'
            path.write_text(json.dumps(config))
            return main.load_config(path)

    def test_expands_same_repeat_count_for_each_group_in_order(self):
        config = self.load(['1000f', '500r', '1500f'], 3)
        self.assertEqual(config['trials'], [
            '1000f1', '1000f2', '1000f3',
            '500r1', '500r2', '500r3',
            '1500f1', '1500f2', '1500f3',
        ])
        self.assertEqual(main.select_trials(config, 500, 'r'), ['500r1', '500r2', '500r3'])
        self.assertEqual(main.select_trials(config, 1000, 'r'), [])
        self.assertEqual(self.load(['1000f', '500r'], 1)['trials'], ['1000f1', '500r1'])
        # Counts are numeric integers in MATLAB, so accept JSON 3.0 as well.
        self.assertEqual(self.load(['500r'], 3.0)['trials'], ['500r1', '500r2', '500r3'])

    def test_rejects_malformed_or_duplicate_groups(self):
        for groups in ([], ['1000f1'], ['1000x'], ['0f'], ['1000f', '1000f'], '1000f', [1000]):
            with self.subTest(groups=groups), self.assertRaisesRegex(ValueError, 'group'):
                self.load(groups, 3)

    def test_rejects_nonpositive_or_noninteger_repeat_counts(self):
        for count in (0, -1, 1.5, True, '3', None, []):
            with self.subTest(count=count), self.assertRaisesRegex(ValueError, 'positive integer'):
                self.load(['1000f'], count)


class BatchWorkflowTests(unittest.TestCase):
    def setUp(self):
        # Use a finite, uniformly sampled damped oscillation so the real fitting
        # and torque code can run without recordings or a graphical interaction.
        # All generated CSVs, metadata and plots stay in this temporary folder.
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.addCleanup(lambda: plt.close('all'))
        self.root = Path(self.temp.name)
        self.config = main.load_config()
        self.config['paths']['data'] = self.root
        self.config['paths']['output'] = self.root / 'results'
        self.config['trials'] = ['1000r1']
        self.time = np.arange(1530) * .02
        self.angle = np.rad2deg(.025 * np.exp(-.03 * self.time) * np.cos(2 * self.time))
        self.metadata = dict(valid=True, invalid_frames=[], fit_start_s=10.,
                             audio_candidates_s=dict(startup=8., shutdown=18.),
                             manual_event_times_s=dict(startup=9., shutdown=None),
                             integration_windows_s=dict(startup=[10., 11.], shutdown=[20., 21.]))
        self.write_trial('1000r1', self.metadata)

    def write_trial(self, trial, metadata, angle=None):
        """Create the same CSV/metadata pair that MATLAB generation supplies."""
        (self.root / f'{trial}_metadata.json').write_text(json.dumps(metadata))
        np.savetxt(self.root / f'{trial}_data.csv',
                   np.column_stack([self.time, self.angle if angle is None else angle]),
                   delimiter=',', header='time_seconds,angle_degrees', comments='')

    def read_metadata(self):
        return json.loads((self.root / '1000r1_metadata.json').read_text())

    def analyze(self, **kwargs):
        """Run the actual per-trial analysis while keeping test output concise."""
        with contextlib.redirect_stdout(io.StringIO()):
            return main.analyze_trial('1000r1', self.config, **kwargs)

    def test_selects_explicit_group_in_configured_order(self):
        # Include similar prefixes to catch substring-based group selection.
        self.config['trials'] = ['1000r3', '500r1', '1000f1', '1000r1', '10000r1']
        self.assertEqual(main.select_trials(self.config, 1000, 'r'), ['1000r3', '1000r1'])
        self.assertEqual(main.select_trials(self.config, 500, 'f'), [])
        with self.assertRaisesRegex(ValueError, 'Invalid trial'):
            main.parse_trial('1000x1')

    def test_saved_selections_reuse_manual_timing_without_prompting(self):
        with patch.object(main, 'plot_data', side_effect=AssertionError('Unexpected fit picker')), \
             patch.object(main, 'select_impulse_windows', side_effect=AssertionError('Unexpected window picker')):
            first = self.analyze()
            metadata = self.read_metadata()
            # Correcting a cue should change plot alignment while leaving saved
            # torque integration bounds, measured impulses and torque intact.
            metadata['manual_event_times_s']['startup'] = 7.5
            self.write_trial('1000r1', metadata)
            second = self.analyze()
        self.assertEqual(first['startup_time'], 9.)
        self.assertEqual(second['startup_time'], 7.5)
        self.assertEqual(first['impulses'], second['impulses'])
        np.testing.assert_array_equal(first['torque'], second['torque'])
        self.assertEqual(self.read_metadata()['integration_windows_s'], self.metadata['integration_windows_s'])

    def test_first_selection_starts_from_audio_and_persists_for_rerun(self):
        self.metadata['fit_start_s'] = []  # MATLAB's unset JSON values
        self.metadata['manual_event_times_s']['shutdown'] = []
        self.metadata['integration_windows_s'] = dict(startup=[], shutdown=None)
        self.write_trial('1000r1', self.metadata)
        with patch.object(main, 'plot_data', return_value=10.) as fit, \
             patch.object(main, 'select_impulse_windows', return_value=([9.5, 11.5], [19.5, 21.5])) as picker:
            self.analyze()
        fit.assert_called_once()
        # Manual startup=9 overrides candidate=8; shutdown uses candidate=18.
        # Each default window then adds that event's configured [1, 2] offset.
        self.assertEqual(picker.call_args.args[2:], ([10., 11.], [19., 20.]))
        self.assertEqual(self.read_metadata()['fit_start_s'], 10.)
        with patch.object(main, 'plot_data', side_effect=AssertionError('fit not reused')), \
             patch.object(main, 'select_impulse_windows', side_effect=AssertionError('windows not reused')):
            self.analyze()

    def test_reselect_changes_windows_without_changing_event_times(self):
        with patch.object(main, 'plot_data', return_value=12.), \
             patch.object(main, 'select_impulse_windows', return_value=([9.5, 11.5], [19.5, 21.5])):
            result = self.analyze(reselect=True)
        metadata = self.read_metadata()
        self.assertEqual(metadata['fit_start_s'], 12.)
        self.assertEqual(metadata['audio_candidates_s'], self.metadata['audio_candidates_s'])
        self.assertEqual(metadata['manual_event_times_s'], self.metadata['manual_event_times_s'])
        self.assertEqual(result['startup_time'], 9.)
        self.assertEqual(metadata['integration_windows_s']['startup'], [9.5, 11.5])

    def test_invalid_trials_are_reported_excluded_and_stale_outputs_removed(self):
        # Cover an explicit generation failure, a stale CSV with a NaN angle,
        # and missing inputs before the valid trial to verify continuation.
        self.config['trials'] = ['1000r2', '1000r3', '1000r4', '1000r1']
        invalid = dict(self.metadata, valid=False, invalid_frames=[1, 7], failure_reason='Missing laser')
        self.write_trial('1000r2', invalid)
        bad_angle = self.angle.copy()
        bad_angle[6] = np.nan
        self.write_trial('1000r3', self.metadata, bad_angle)
        # r4 deliberately has no metadata or data.
        output = self.root / 'results/1000r'
        output.mkdir(parents=True)
        stale = ['1000r2_signal.csv', '1000r2_impulses.csv', '1000r3_signal.csv',
                 '1000r2_impulse_diagnostics.png',
                 'impulse_summary.csv', 'torque.png', 'cumulative_integral.png', 'impulse_comparison.png']
        for name in stale:
            (output / name).write_text('stale')
        with contextlib.redirect_stdout(io.StringIO()):
            trials, failures = main.run_group(self.config, 1000, 'r', interactive=False, show=False)
        self.assertEqual([trial['trial_id'] for trial in trials], ['1000r1'])
        self.assertEqual(set(failures), {'1000r2', '1000r3', '1000r4'})
        self.assertIn('[1, 7]', failures['1000r2'])
        self.assertIn('frames [7]', failures['1000r3'])
        self.assertFalse((output / '1000r2_signal.csv').exists())
        self.assertFalse((output / '1000r2_impulses.csv').exists())
        self.assertFalse((output / '1000r3_signal.csv').exists())
        self.assertFalse((output / '1000r2_impulse_diagnostics.png').exists())
        with (output / 'impulse_summary.csv').open() as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 2)
        self.assertEqual({row['run_identity'] for row in rows}, {'1000r1'})
        self.assertEqual({row['event'] for row in rows}, {'startup', 'shutdown'})
        for name in ('torque.png', 'cumulative_integral.png', 'impulse_comparison.png',
                     '1000r1_impulse_diagnostics.png'):
            self.assertGreater((output / name).stat().st_size, 100)
        # A subsequent all-failed run must not leave the earlier valid products.
        self.write_trial('1000r1', invalid)
        with contextlib.redirect_stdout(io.StringIO()):
            trials, failures = main.run_group(self.config, 1000, 'r', interactive=False, show=False)
        self.assertFalse(trials)
        self.assertEqual(len(failures), 4)
        self.assertEqual(list(output.iterdir()), [])


if __name__ == '__main__':
    unittest.main()
