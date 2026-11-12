"""Exercise the torque picker through Matplotlib's mouse/key callbacks."""

import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backend_bases import KeyEvent, MouseEvent
import numpy as np

from Main_Functions.impulse_funcs import analyze_impulses, select_impulse_windows


class ImpulseSelectionTests(unittest.TestCase):
    def setUp(self):
        self.time = np.linspace(0, 30, 241)
        self.torque = np.full_like(self.time, 2.)
        self.defaults = ((11., 12.), (21., 22.))

    def pick(self, actions, view_windows=((9, 13), (19, 23))):
        def interact(*, block):
            self.assertTrue(block)
            fig = plt.gcf()
            fig.canvas.draw()
            # Verify the selectable curves are torque, rather than its integral.
            for ax in fig.axes[:2]:
                curve = ax.lines[0]
                np.testing.assert_allclose(
                    curve.get_ydata(),
                    np.interp(curve.get_xdata(), self.time, self.torque),
                )
                self.assertEqual(ax.get_ylabel(), "Extracted torque (dyn*cm)")
            for action in actions:
                if action[0] == "open":
                    self.assertTrue(plt.fignum_exists(fig.number))
                    continue
                elif action[0] == "key":
                    event = KeyEvent("key_press_event", fig.canvas, key=action[1])
                else:
                    _, panel, time = action
                    ax = fig.axes[panel]
                    ax.set_xlim(0, 30)
                    x_pixel, y_pixel = ax.transData.transform((time, np.mean(ax.get_ylim())))
                    event = MouseEvent("button_press_event", fig.canvas,
                                       x_pixel, y_pixel, button=1)
                fig.canvas.callbacks.process(event.name, event)

        with patch.object(plt, "show", side_effect=interact):
            return select_impulse_windows(
                self.time, self.torque, *self.defaults,
                view_windows=view_windows,
            )

    def assertWindowsEqual(self, actual, expected):
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)

    def test_close_without_clicks_and_incomplete_pairs_keep_sound_windows(self):
        self.assertWindowsEqual(self.pick([]), self.defaults)
        self.assertWindowsEqual(
            self.pick([("click", 0, 10), ("click", 1, 20)]), self.defaults,
        )

    def test_completed_pairs_in_either_order_feed_analysis(self):
        windows = self.pick([
            ("click", 0, 12.5), ("click", 0, 9.5),
            ("click", 1, 19.5), ("click", 1, 22.5), ("key", "enter"),
        ])
        self.assertWindowsEqual(windows, ((9.5, 12.5), (19.5, 22.5)))
        results = analyze_impulses(
            self.time, self.torque, self.time,
            startup_window=windows[0], shutdown_window=windows[1],
        )
        for row in results:
            self.assertAlmostEqual(row["J_measured_dyn_cm_s"], 6.)

    def test_completed_event_survives_other_events_incomplete_pair(self):
        windows = self.pick([
            ("click", 0, 9.5), ("click", 0, 12.5), ("click", 1, 20),
        ])
        self.assertWindowsEqual(windows, ((9.5, 12.5), self.defaults[1]))

    def test_escape_restores_both_sound_windows(self):
        windows = self.pick([
            ("click", 0, 9.5), ("click", 0, 12.5), ("key", "escape"),
        ])
        self.assertWindowsEqual(windows, self.defaults)

    def test_overlapping_or_zero_width_pair_can_be_retried(self):
        windows = self.pick([
            ("click", 0, 10), ("click", 0, 10),
            ("click", 0, 23), ("click", 0, 13),
        ])
        self.assertWindowsEqual(windows, ((10, 13), self.defaults[1]))

    def test_incomplete_replacement_falls_back_to_sound(self):
        windows = self.pick([
            ("click", 0, 9.5), ("click", 0, 12.5), ("click", 0, 10),
        ])
        self.assertWindowsEqual(windows, self.defaults)

    def test_invalid_suggestions_can_be_replaced(self):
        for defaults in (
            ((-5, 2), (29, 35)),
            ((9, 23), (19, 22)),
            ((12, 11), (21, 21)),
            ((np.nan, 12), (21, np.inf)),
        ):
            with self.subTest(defaults=defaults):
                self.defaults = defaults
                # Repair shutdown first even when it overlaps the uncorrected
                # startup suggestion. Final acceptance still needs both valid.
                windows = self.pick([
                    ("click", 1, 19.5), ("click", 1, 22.5),
                    ("click", 0, 9.5), ("click", 0, 12.5), ("key", "enter"),
                ], view_windows=None)
                self.assertWindowsEqual(windows, ((9.5, 12.5), (19.5, 22.5)))

    def test_valid_suggestion_can_be_retained_while_other_is_repaired(self):
        self.defaults = ((11., 12.), (29., 35.))
        windows = self.pick([("click", 1, 20), ("click", 1, 23)])
        self.assertWindowsEqual(windows, ((11., 12.), (20., 23.)))

    def test_invalid_suggestions_without_valid_replacements_raise(self):
        for defaults in (((-5, 2), (29, 35)), ((9, 23), (19, 22))):
            for actions in ([], [("click", 0, 10)], [("key", "enter")], [
                ("click", 0, 9.5), ("click", 0, 12.5),
                ("click", 1, 19.5), ("click", 1, 22.5), ("key", "escape"),
            ]):
                with self.subTest(defaults=defaults, actions=actions):
                    self.defaults = defaults
                    with self.assertRaisesRegex(ValueError, 'integration window'):
                        self.pick(actions)
                    self.assertEqual(plt.get_fignums(), [])

    def test_invalid_done_keeps_picker_open_for_repair(self):
        self.defaults = ((-5, 2), (29, 35))
        windows = self.pick([
            ("key", "enter"), ("open",),
            ("click", 0, 9.5), ("click", 0, 12.5),
            ("click", 1, 19.5), ("click", 1, 22.5), ("key", "enter"),
        ], view_windows=((40, 50), (40, 50)))
        self.assertWindowsEqual(windows, ((9.5, 12.5), (19.5, 22.5)))


if __name__ == "__main__":
    unittest.main()
