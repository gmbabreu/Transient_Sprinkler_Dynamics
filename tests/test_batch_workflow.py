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

    def test_expands_individual_repeat_counts_in_group_order(self):
        config = self.load(['1000f', '1000r', '500f', '500r'], [1, 2, 4, 5])
        self.assertEqual(config['trials'], [
            '1000f1',
            '1000r1', '1000r2',
            '500f1', '500f2', '500f3', '500f4',
            '500r1', '500r2', '500r3', '500r4', '500r5',
        ])
        self.assertEqual(main.select_trials(config, 1000, 'r'), ['1000r1', '1000r2'])
        self.assertEqual(self.load(['1000f'], [3])['trials'], ['1000f1', '1000f2', '1000f3'])
        self.assertEqual(self.load(['1000f', '500r'], [1., 2.])['trials'],
                         ['1000f1', '500r1', '500r2'])

    def test_rejects_mismatched_or_invalid_count_arrays(self):
        for counts in ([], [1], [1, 2, 3], [0, 2], [-1, 2], [1.5, 2],
                       [True, 2], ['1', 2], [None, 2], [[1], 2]):
            with self.subTest(counts=counts), self.assertRaisesRegex(ValueError, 'trials_per_group'):
                self.load(['1000f', '500r'], counts)

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

    def test_event_times_must_be_ordered_within_recording(self):
        for source in ('audio_candidates_s', 'manual_event_times_s'):
            for startup, shutdown in ((-1., 18.), (18., 18.), (19., 18.),
                                      (8., self.time[-1] + 1), (8., -1.)):
                with self.subTest(source=source, startup=startup, shutdown=shutdown):
                    metadata = json.loads(json.dumps(self.metadata))
                    metadata['manual_event_times_s'] = dict(startup=None, shutdown=None)
                    metadata[source] = dict(startup=startup, shutdown=shutdown)
                    self.write_trial('1000r1', metadata)
                    path = self.root / '1000r1_metadata.json'
                    before = path.read_text()
                    with self.assertRaisesRegex(ValueError, 'Invalid effective event times') as caught:
                        self.analyze()
                    self.assertIn(f'startup={startup:g}', str(caught.exception))
                    self.assertIn(f'shutdown={shutdown:g}', str(caught.exception))
                    self.assertEqual(path.read_text(), before)

    def test_event_times_at_recording_endpoints_are_valid(self):
        self.metadata['manual_event_times_s'] = dict(startup=self.time[0], shutdown=self.time[-1])
        self.write_trial('1000r1', self.metadata)
        self.assertEqual(self.analyze(interactive=False)['startup_time'], self.time[0])

    def test_invalid_windows_are_not_saved_without_valid_manual_selection(self):
        windows = dict(startup=[-1., 1.], shutdown=[29., 35.])
        self.metadata['integration_windows_s'] = windows
        self.write_trial('1000r1', self.metadata)
        for interactive in (False, True):
            with self.subTest(interactive=interactive), \
                 patch.object(plt, 'show', return_value=None), \
                 patch.object(main, 'select_impulse_windows', wraps=main.select_impulse_windows) as picker:
                with self.assertRaisesRegex(ValueError, 'integration window'):
                    self.analyze(interactive=interactive)
                self.assertEqual(picker.called, interactive)
                self.assertEqual(self.read_metadata()['integration_windows_s'], windows)

    def test_invalid_saved_windows_prompt_for_repair(self):
        self.metadata['integration_windows_s'] = dict(startup=[9., 23.], shutdown=[19., 22.])
        self.write_trial('1000r1', self.metadata)
        with patch.object(main, 'select_impulse_windows', return_value=([10., 11.], [20., 21.])) as picker:
            self.analyze()
        picker.assert_called_once()
        self.assertEqual(self.read_metadata()['integration_windows_s'],
                         dict(startup=[10., 11.], shutdown=[20., 21.]))

    def test_angular_plot_refreshes_after_manual_timing_correction(self):
        first = self.analyze(interactive=False)
        metadata = self.read_metadata()
        metadata['manual_event_times_s']['startup'] = 7.5
        self.write_trial('1000r1', metadata)
        second = self.analyze(interactive=False)
        for trial in (first, second):
            with patch.object(plt.Figure, 'savefig', autospec=True) as save:
                main.plot_group_results([trial], self.root, '1000r')
            fig = next(call.args[0] for call in save.call_args_list
                       if Path(call.args[1]).name == 'angular_displacement.png')
            self.assertEqual(len(fig.axes[0].lines), 1)
            curve = fig.axes[0].lines[0]
            np.testing.assert_array_equal(curve.get_xdata(), self.time - trial['startup_time'])
            np.testing.assert_array_equal(curve.get_ydata(), np.rad2deg(trial['phi']))
            self.assertEqual(curve.get_label(), '1000r1')
            plt.close('all')
        np.testing.assert_array_equal(first['torque'], second['torque'])
        self.assertEqual(first['impulses'], second['impulses'])

    def test_reselect_changes_windows_without_changing_event_times(self):
        with patch.object(main, 'plot_data', return_value=12.) as fit, \
             patch.object(main, 'select_impulse_windows', return_value=([9.5, 11.5], [19.5, 21.5])) as picker:
            result = self.analyze(reselect=True)
        self.assertEqual(fit.call_args.kwargs['initial_time'], 10.)
        self.assertEqual(picker.call_args.args[2:], ([10., 11.], [20., 21.]))
        metadata = self.read_metadata()
        self.assertEqual(metadata['fit_start_s'], 12.)
        self.assertEqual(metadata['audio_candidates_s'], self.metadata['audio_candidates_s'])
        self.assertEqual(metadata['manual_event_times_s'], self.metadata['manual_event_times_s'])
        self.assertEqual(result['startup_time'], 9.)
        self.assertEqual(metadata['integration_windows_s']['startup'], [9.5, 11.5])
        # A second reselection must use the most recently accepted choices.
        with patch.object(main, 'plot_data', return_value=12.) as fit, \
             patch.object(main, 'select_impulse_windows', return_value=([9.5, 11.5], [19.5, 21.5])) as picker:
            self.analyze(reselect=True)
        self.assertEqual(fit.call_args.kwargs['initial_time'], 12.)
        self.assertEqual(picker.call_args.args[2:], ([9.5, 11.5], [19.5, 21.5]))

    def test_fit_picker_opens_at_saved_time(self):
        def inspect():
            marker = plt.gcf().axes[0].lines[1]
            np.testing.assert_array_equal(marker.get_xdata(), [17.25, 17.25])
        with patch.object(plt, 'show', side_effect=inspect):
            chosen = main.plot_data(self.time, np.deg2rad(self.angle), 1, initial_time=17.25)
        self.assertEqual(chosen, 17.25)

    def test_impulse_comparison_preserves_values_and_labels_each_trial(self):
        first = self.analyze(interactive=False)
        trials = [first, dict(first, trial_id='1000r2')]
        with patch.object(plt.Figure, 'savefig', autospec=True) as save:
            main.plot_group_results(trials, self.root, '1000r')
        fig = next(call.args[0] for call in save.call_args_list
                   if Path(call.args[1]).name == 'cumulative_jump_comparison.png')
        self.assertEqual(len(fig.axes), 2)
        self.assertEqual([text.get_text() for text in fig.axes[1].get_xticklabels()],
                         ['1000r1', '1000r2'])
        for event_index, event in enumerate(('startup', 'shutdown')):
            measured = next(line for line in fig.axes[0].lines
                            if line.get_label() == event.capitalize())
            theory = next(line for line in fig.axes[0].lines
                          if line.get_label() == event.capitalize() + ' theory')
            rows = [next(row for row in trial['impulses'] if row['event'] == event)
                    for trial in trials]
            np.testing.assert_array_equal(measured.get_ydata(),
                                          [row['J_measured_dyn_cm_s'] for row in rows])
            np.testing.assert_array_equal(theory.get_ydata(),
                                          [rows[0]['J_theory_dyn_cm_s']] * 2)
            segments = fig.axes[0].collections[3 * event_index].get_segments()
            for segment, row in zip(segments, rows):
                np.testing.assert_array_equal(segment[:, 1],
                    [row['sensitivity_min_dyn_cm_s'], row['sensitivity_max_dyn_cm_s']])
        net = next(line for line in fig.axes[1].lines if line.get_marker() == 'o')
        np.testing.assert_array_equal(net.get_ydata(),
            [sum(row['J_measured_dyn_cm_s'] for row in trial['impulses']) for trial in trials])
        for ax in fig.axes:
            self.assertTrue(any(np.array_equal(line.get_ydata(), [0, 0]) for line in ax.lines))
        self.assertEqual(len(fig.legends[0].texts), 4)
        labels = [text.get_text() for text in fig.axes[0].texts]
        metrics = main.relative_error_metrics(first['impulses'])
        self.assertIn(f"Error {metrics['startup_relative_error_pct']:.2f}%", labels)
        self.assertIn(f"Error {metrics['shutdown_relative_error_pct']:.2f}%", labels)
        self.assertIn(f"Mismatch {metrics['relative_magnitude_mismatch_pct']:.3g}%",
                      [text.get_text() for text in fig.axes[1].texts])
        startup = next(row for row in first['impulses'] if row['event'] == 'startup')
        np.testing.assert_array_equal(fig.axes[0].texts[0].xy,
                                      [-.08, startup['J_measured_dyn_cm_s']])
        for figure_call in save.call_args_list:
            figure = figure_call.args[0]
            labels = [figure._suptitle.get_text()] if figure._suptitle is not None else []
            for axis in figure.axes:
                labels.extend([axis.get_title(), axis.get_xlabel(), axis.get_ylabel()])
                labels.extend(text.get_text() for text in axis.texts)
            self.assertNotIn('impulse', ' '.join(labels).lower())
        self.assertNotIn('relative_errors.png', [Path(call.args[1]).name for call in save.call_args_list])

    def test_invalid_trials_are_reported_excluded_and_stale_outputs_removed(self):
        # Cover an explicit generation failure, a stale CSV with a NaN angle,
        # and missing inputs before the valid trial to verify continuation.
        self.config['trials'] = ['1000r2', '1000r3', '1000r4', '1000r5', '1000r1']
        invalid = dict(self.metadata, valid=False, invalid_frames=[1, 7], failure_reason='Missing laser')
        self.write_trial('1000r2', invalid)
        bad_angle = self.angle.copy()
        bad_angle[6] = np.nan
        self.write_trial('1000r3', self.metadata, bad_angle)
        invalid_events = dict(self.metadata, manual_event_times_s=dict(startup=20., shutdown=10.))
        self.write_trial('1000r5', invalid_events)
        # r4 deliberately has no metadata or data.
        output = self.root / 'results/1000r'
        output.mkdir(parents=True)
        stale = ['1000r2_signal.csv', '1000r2_jumps.csv', '1000r3_signal.csv',
                 '1000r2_jump_diagnostics.png',
                 'jump_summary.csv', 'torque.png', 'cumulative_integral.png', 'cumulative_jump_comparison.png',
                 'angular_displacement.png', 'relative_errors.png', 'relative_errors.csv',
                 'impulse_comparison.png', 'impulse_summary.csv', '1000r2_impulses.csv',
                 '1000r2_impulse_diagnostics.png']
        for name in stale:
            (output / name).write_text('stale')
        with contextlib.redirect_stdout(io.StringIO()):
            trials, failures = main.run_group(self.config, 1000, 'r', interactive=False, show=False)
        self.assertEqual([trial['trial_id'] for trial in trials], ['1000r1'])
        self.assertEqual(set(failures), {'1000r2', '1000r3', '1000r4', '1000r5'})
        self.assertIn('Invalid effective event times', failures['1000r5'])
        self.assertIn('[1, 7]', failures['1000r2'])
        self.assertIn('frames [7]', failures['1000r3'])
        self.assertFalse((output / '1000r2_signal.csv').exists())
        self.assertFalse((output / '1000r2_jumps.csv').exists())
        self.assertFalse((output / '1000r3_signal.csv').exists())
        self.assertFalse((output / '1000r2_jump_diagnostics.png').exists())
        with (output / 'jump_summary.csv').open() as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 2)
        self.assertEqual({row['run_identity'] for row in rows}, {'1000r1'})
        self.assertEqual({row['event'] for row in rows}, {'startup', 'shutdown'})
        self.assertFalse(any('impulse' in path.name for path in output.iterdir()))
        self.assertFalse((output / 'relative_errors.png').exists())
        with (output / 'relative_errors.csv').open() as stream:
            errors = list(csv.DictReader(stream))
        self.assertEqual(len(errors), 1)
        self.assertEqual(errors[0]['trial_id'], '1000r1')
        for field, value in main.relative_error_metrics(trials[0]['impulses']).items():
            self.assertEqual(float(errors[0][field]), value)
        for name in ('torque.png', 'cumulative_integral.png', 'cumulative_jump_comparison.png',
                     '1000r1_jump_diagnostics.png', 'angular_displacement.png'):
            self.assertGreater((output / name).stat().st_size, 100)
        # A subsequent all-failed run must not leave the earlier valid products.
        self.write_trial('1000r1', invalid)
        with contextlib.redirect_stdout(io.StringIO()):
            trials, failures = main.run_group(self.config, 1000, 'r', interactive=False, show=False)
        self.assertFalse(trials)
        self.assertEqual(len(failures), 5)
        self.assertEqual(list(output.iterdir()), [])


if __name__ == '__main__':
    unittest.main()
