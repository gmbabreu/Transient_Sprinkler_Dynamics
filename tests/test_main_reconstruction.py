"""Run main's real reconstruction/analysis/export with deterministic CSV input.

Only user configuration and the interactive fit-time selector are replaced.
Like main.py, these integration checks require NumPy 2 or newer.
"""

import ast
import contextlib
import importlib
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
MAIN = REPO / "main.py"


class MainReconstructionTests(unittest.TestCase):
    def test_export_matches_original_time_and_impulse_analysis_torque(self):
        for count in (1530, 2048, 2049, 5001):
            with self.subTest(count=count), tempfile.TemporaryDirectory() as directory:
                time = 5. + np.arange(count) * .02
                angle = .025 * np.exp(-.03 * (time - 5.)) * np.cos(2 * (time - 5.))
                np.savetxt(Path(directory) / "1000r1_data.csv",
                           np.column_stack([time, np.rad2deg(angle)]),
                           delimiter=",", header="time_seconds,angle_degrees", comments="")
                # Override editable settings, then execute the actual main.py
                # body, including fitting, filtering, forward solve and impulses.
                tree = ast.parse(MAIN.read_text())
                settings = dict(data_dir=directory, write_t=1, plot_switch=0,
                                proc_data_switch=0, write_impulse_summary=0)
                for node in tree.body:
                    if isinstance(node, ast.Assign) and len(node.targets) == 1:
                        target = node.targets[0]
                        if isinstance(target, ast.Name) and target.id in settings:
                            node.value = ast.Constant(settings[target.id])
                ast.fix_missing_locations(tree)
                namespace = {"__file__": str(MAIN), "__name__": "__main__"}
                with patch.object(sys, "path", sys.path + [str(REPO / "Main_Functions")]):
                    reader = importlib.import_module("dataread_funcs")
                    impulses = importlib.import_module("impulse_funcs")
                    with patch.object(reader, "plot_data", return_value=10.), \
                         patch.object(impulses, "analyze_impulses", wraps=impulses.analyze_impulses) as analyze, \
                         contextlib.redirect_stdout(io.StringIO()):
                        exec(compile(tree, str(MAIN), "exec"), namespace)
                try:
                    export = Path(directory) / "1000r1_data_signal.csv"
                    self.assertEqual(export.read_text().splitlines()[0], "time_seconds,torque_dyn_cm")
                    saved = np.loadtxt(export, delimiter=",", skiprows=1)
                    self.assertEqual(saved.shape, (count, 2))
                    analyze.assert_called_once()
                    analyzed_time, analyzed_torque, original_time = analyze.call_args.args
                    np.testing.assert_array_equal(saved[:, 0], time)
                    np.testing.assert_array_equal(saved[:, 0], original_time)
                    np.testing.assert_array_equal(saved[:, 0], analyzed_time)
                    np.testing.assert_array_equal(saved[:, 1], analyzed_torque)
                    np.testing.assert_array_equal(saved[:, 1], namespace["true_torque"])
                    np.testing.assert_array_equal(
                        saved[:, 1], namespace["extracted_torque_seg"] * namespace["I"])
                    self.assertFalse(np.allclose(saved[:, 1], namespace["extracted_torque_seg"]))
                    self.assertEqual(namespace["N_f"], len(namespace["fourier_phi"]))
                    self.assertEqual(len(namespace["phi_forward"]), count)
                    self.assertTrue(np.isfinite(namespace["phi_forward"]).all())
                finally:
                    plt.close("all")


if __name__ == "__main__":
    unittest.main()
