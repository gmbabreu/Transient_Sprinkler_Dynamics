"""Regression checks for recording validation and FFT padding boundaries."""

import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import cheby2, sosfiltfilt

from Main_Functions.fft_funcs import torque_solver
from Main_Functions.plot_funcs import plot_franken
from Main_Functions.process_funcs import combine_data


class PaddingTests(unittest.TestCase):
    def test_lengths_preserve_every_original_sample_and_timestamp(self):
        cases = ((17, 4096), (2047, 4096), (2048, 4096),
                 (2049, 8192), (4096, 8192), (4097, 16384),
                 (5001, 16384), (8192, 16384))
        for count, total in cases:
            with self.subTest(count=count):
                time = 7.25 + np.arange(count) * .03125
                signal = .3 + np.sin(time)
                padded_t, padded_y, added = combine_data(time, signal, None)
                self.assertEqual(added, total - count)
                self.assertEqual(len(padded_t), total)
                self.assertEqual(len(padded_y), total)
                np.testing.assert_array_equal(padded_t[added:], time)
                np.testing.assert_array_equal(padded_y[added:], signal)
                np.testing.assert_array_equal(padded_y[:added], signal[0])
                np.testing.assert_array_equal(
                    padded_t[:added], time[0] + np.arange(-added, 0) * .03125)
                self.assertEqual(padded_t[added - 1], time[0] - .03125)
                self.assertTrue(np.all(np.diff(padded_t) > 0))
                forcing = torque_solver(total, .2, 2., padded_t[-1] - padded_t[0], padded_y)
                self.assertEqual(forcing.shape, padded_y.shape)
                self.assertTrue(np.isfinite(forcing).all())

    def test_legacy_short_recording_padding_and_filtered_reconstruction(self):
        time = np.arange(1530) * .125
        signal = .4 + np.exp(-.02 * time) * np.cos(2 * time)
        # Independent reference to the previously supported <2048 code path.
        legacy_added = 2048 + (2048 - len(signal))
        legacy_t = np.concatenate([
            np.arange(-legacy_added * .125, 0, .125), time[:2048]])
        legacy_y = np.concatenate([
            np.full(legacy_added, signal[0]), signal[:2048]])
        padded_t, padded_y, added = combine_data(time, signal, None)
        self.assertEqual(added, legacy_added)
        np.testing.assert_array_equal(padded_t, legacy_t)
        np.testing.assert_array_equal(padded_y, legacy_y)
        # Keep the legacy frequency ordering, period and transfer convention.
        frequencies = (2 * np.pi / (legacy_t[-1] - legacy_t[0])) * np.r_[
            np.arange(2048), [0], np.arange(-2047, 0)]
        response = -1 / ((-frequencies)**2 + 2j * .2 * (-frequencies) - 2.**2)
        legacy_torque = np.fft.ifft(np.fft.fft(legacy_y) / response).real
        actual = torque_solver(len(padded_y), .2, 2., padded_t[-1] - padded_t[0], padded_y).real
        np.testing.assert_allclose(actual, legacy_torque, rtol=1e-13, atol=1e-13)
        sos = cheby2(N=4, rs=60, Wn=.5, btype="low", output="sos")
        np.testing.assert_allclose(sosfiltfilt(sos, actual),
                                   sosfiltfilt(sos, legacy_torque),
                                   rtol=1e-13, atol=1e-13)

    def test_floating_interval_has_exact_padding_count_and_preserves_time(self):
        time = 12.3 + np.arange(1515) * (1 / 29.97)
        padded_t, _, added = combine_data(time, np.sin(time), None)
        self.assertEqual(len(padded_t), 4096)
        self.assertEqual(padded_t[added - 1], time[0] - np.mean(np.diff(time)))
        np.testing.assert_array_equal(padded_t[added:], time)

    def test_small_timestamp_roundoff_is_supported_without_resampling(self):
        time = np.round(3. + np.arange(100) / 29.97, 8)
        padded_t, _, added = combine_data(time, np.sin(time), None)
        np.testing.assert_array_equal(padded_t[added:], time)

    def test_invalid_recordings_report_clear_errors(self):
        cases = (
            ([0, 1], [0], "matching 1-D"),
            ([[0, 1]], [[0, 1]], "matching 1-D"),
            ([], [], "at least two"),
            ([0], [1], "at least two"),
            ([0, np.nan], [0, 1], "finite"),
            ([0, 1], [0, np.inf], "finite"),
            ([0, 0, 1], [0, 1, 0], "strictly increasing"),
            ([0, 2, 1], [0, 1, 0], "strictly increasing"),
            ([0, 1, 2.1], [0, 1, 0], "uniform sampling"),
            ([0, 1j], [0, 1], "real-valued"),
            ([0, 1], [0, 1j], "real-valued"),
            ([0, "invalid"], [0, 1], "numeric"),
            ([[0], [1, 2]], [0, 1], "numeric"),
        )
        for time, signal, message in cases:
            with self.subTest(time=time, signal=signal):
                with self.assertRaisesRegex(ValueError, message):
                    combine_data(time, signal, None)

    def test_padding_that_cannot_be_represented_is_rejected(self):
        # Original timestamps are distinct, but extending backward at this
        # offset crosses into a float range where that interval is too small.
        time = np.array([-2.**53, -2.**53 + 1])
        with self.assertRaisesRegex(ValueError, "padded timestamps"):
            combine_data(time, [1, 2], None)

    def test_padding_plot_uses_actual_boundary(self):
        for count in (1530, 2048, 5001):
            with self.subTest(count=count):
                time = np.arange(count) * .125
                signal = .4 + np.sin(time)
                padded_t, padded_y, added = combine_data(time, signal, None)
                try:
                    plot_franken(padded_t, padded_y, time, 0, None)
                    padding, original = plt.gca().lines
                    self.assertEqual(len(padding.get_xdata()), added)
                    np.testing.assert_array_equal(original.get_xdata(), time)
                    np.testing.assert_array_equal(original.get_ydata(), signal)
                finally:
                    plt.close("all")


if __name__ == "__main__":
    unittest.main()
