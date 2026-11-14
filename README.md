# Transient Sprinkler Dynamics

Work by Rachel Bertaud, Gabriel Abreu, and Kelly Yu.

Reconstruct torque from each trial's angular displacement by fitting its free-decay
tail, applying the existing Fourier solver and torque filter, and checking a
forward ODE reconstruction. Measure two separate signed cumulative torque values around pump
startup and shutdown, and compare individual trials against theory.

## Four-step workflow

Run these commands from the repository root.

1. **Configure ingestion.** Keep recordings directly under `paths.trials` on the
   existing external drive, with exact names such as `1000r5.MTS`. Edit
   [`experiment_config.json`](experiment_config.json):

   ```json
   "trial_ranges": {
     "500r": [], "500f": [], "1000r": [5, 9],
     "1000f": [], "1500r": [], "1500f": []
   }
   ```

   Endpoints are inclusive: this ingests five recordings, `1000r5`–`1000r9`.
   `[]` skips ingestion for that group and keeps its history available for analysis.
   Endpoints must be positive integers in ascending order; malformed groups,
   duplicate keys, nested arrays and invalid ranges are rejected. Separate batches
   may leave gaps in numbering: historical IDs are never renumbered.

   Paths without a leading `/` are relative to the configuration file. Keep
   `paths.trials` and `paths.frames` on their configured drives. Scripts never
   rename, move or modify original MTS files. Generic camera names such as
   `00001.MTS` require an explicit source-identification mechanism; file ordering
   never assigns a trial ID. The default extension remains `.MTS`.

2. **Extract frames and calibrate (MATLAB).** Run extraction, then click each
   ruler mark requested at the top of the calibration window. Points stay
   labeled and can be dragged to correct them. After placing all marks, click a
   point to select it (cyan): arrow keys move it 1 pixel; Shift+arrow moves it
   10 pixels. Click **Accept points** after review, or close the window to
   cancel without saving calibration parameters.
   Set `calibration.reference_trial` explicitly (for example `1000r1`, even if
   its ingestion range is empty). Alternatively, omit `reference_trial` and
   provide an image at `paths.calibration_reference`. Completing the calibration
   replaces the MAT file at `paths.ruler_calibration`; cancelling keeps the
   previous calibration. The reference image is not replaced.

   ```matlab
   run('Data_Generation/save_frames.m')
   run('Data_Generation/calibration/calibrateRulerAxis.m')
   ```

   Skip the calibration command when reusing the saved calibration.
   Completed frame extractions retain their `.complete` validation checks.

3. **Generate trial data (MATLAB).**

   ```matlab
   run('Data_Generation/generateData.m')
   ```

   This processes only ingestion trials, skips completed CSV/metadata pairs,
   and saves angle CSVs, audio-event candidates and per-trial metadata.
   Use `audio.inspect` to review the cues; put corrections in the trial's
   `manual_event_times_s` metadata fields. Failed trials are reported and skipped.

4. **Analyze a group (Python).**

   ```bash
   python main.py 1000 r
   ```

   On the first run, choose the free-decay fit start and two integration
   endpoints for each event. Selections are saved and reused next time.
   Results and plots go to `results/1000r/` by default.

| Optional argument | Use |
| --- | --- |
| `--reselect` | Reopen the fit start and integration windows at the latest saved selections |
| `--no-show` | Save plots without showing final diagnostics; selection remains interactive |
| `--non-interactive` | Use a saved fit start and valid saved/audio-based windows without prompting |
| `--config path/to/config.json` | Use a different batch configuration |

The named form `python main.py --re 1000 --direction r` also works.

MATLAB scripts add their own helper directories. Run them by full path or add
`Data_Generation` and `Data_Generation/calibration` to MATLAB's search path.
MATLAB with Image Processing Toolbox and supported video/audio codecs is required.
Python requires `numpy`, `scipy`, and `matplotlib`. NumPy 1.x and 2.x are
supported; trapezoidal integration uses SciPy.

```bash
pip install numpy scipy matplotlib
python -m unittest discover -s tests -v
```

## Settings and saved selections

All shared batch settings live in `experiment_config.json`:

| Setting | Purpose |
| --- | --- |
| `paths` | Inputs, generated data, calibration files and group results |
| `trial_ranges` | Groups such as `1000f`: inclusive `[first, last]` ingestion ranges, or `[]` |
| `calibration.reference_trial` | Optional explicit reference trial, independent of ingestion ranges |
| `video` | Extension, mirror distance in inches, laser YCbCr limits, extraction regeneration and tracking diagnostics |
| `physics` | Spring stiffness, density, flow conversion, signed one-arm geometry factor and sign convention |
| `analysis` | Existing fit switch, Chebyshev torque-filter settings and diagnostic display |
| `audio` | Separate search intervals in **audio seconds**, audio-to-video offset, inspection and playback |
| `integration_window_offsets_s` | Initial bounds relative to each effective event time in **video seconds** |

To deliberately regenerate a trial, include it in its ingestion range and in
`video.regenerate_trials`, for example `["1000r5"]`. This allowlist permits
`generateData.m` to replace that trial's CSV. To also rebuild its frames, set
`video.regenerate_frames` to `true`. Remove the allowlist entry and restore the
flag afterward. Regeneration invalidates the old CSV and retains manual audio
corrections, fit starts, integration windows and other metadata. A failed
regeneration leaves that trial invalid; successful generation restores validity.
Unrelated trials and group results are retained. Make a separate copy of any
trial you want to keep before deliberately regenerating it. Incomplete extraction
folders can be retried, but replacing existing trial data still requires the
allowlist. Run ingestion and migration serially.

The existing YCbCr laser thresholds are unchanged and editable under `video.laser_ycbcr_limits`.
[`Data_Generation/ycbcrMask.m`](Data_Generation/ycbcrMask.m) applies these limits
when isolating the laser in each frame.
Spring calibration is a different physical setup: its geometry and measured
weights remain in `Data_Generation/calibration/calibrateSpring.m`. Copy its
resulting stiffness into `physics.kappa_dyn_cm_per_rad`; do not substitute the
video geometry into that script.

`audio.inspect` shows the original waveform, 20 ms RMS and angular-position
panels for each event, sharing the original video time axis.
`audio.playback_enabled` auditions selected channels at `audio.playback_gain`,
with a moving cursor on the waveform.
Enable these settings in `experiment_config.json` and run
`Data_Generation/generateData.m` to inspect or replay the audio. It calls
[`Data_Generation/estimateSoundTimes.m`](Data_Generation/estimateSoundTimes.m)
and saves candidates separately from manual event corrections, as video seconds:
`video_time = audio_time + audio_to_video_offset_s`.

Each `data/{group}/{trial}_metadata.json` contains:

- `valid`, `failure_reason`, `invalid_frames`: tracking status and 1-based frame numbers.
- `audio_candidates_s`: independently estimated `startup` and `shutdown` times.
- `manual_event_times_s`: optional `startup` and `shutdown` corrections in original video seconds.
- `fit_start_s`: selected free-decay start in original seconds.
- `integration_windows_s`: selected `startup` and `shutdown` pairs in original seconds.

The picker can repair invalid suggested windows. Accepting requires two
positive-duration, nonoverlapping windows inside the recording. Incomplete
pairs retain valid suggestions; closing/cancelling without valid bounds fails
without overwriting saved windows. Noninteractive runs reject invalid bounds.

Unset selections may be JSON `null` or MATLAB's empty array `[]`. For example,
edit only the appropriate manual correction field to confirm an audio cue:

```json
"manual_event_times_s": {"startup": 10.2, "shutdown": null}
```

The effective event time is the manual correction when present, otherwise the
candidate. Ordinary regeneration preserves manual corrections, fit start and
integration selections. Regenerating frames invalidates the old tracking data
until generation succeeds again. Missing/invalid metadata and nonfinite input
angles are rejected by Python, even if a stale CSV exists.
Both MATLAB generation and Python analysis require effective times satisfying
`recording_start <= startup < shutdown <= recording_end`; invalid trials are reported and skipped.

Comparison curves use `aligned_time = original_time - effective_startup_time`
on each trial's own samples. Saved data and integration selections retain original
timestamps. Changing an integration window never changes audio timing or plot
alignment. Python analysis refreshes all comparison plots, including angular
displacement, from successful trials using the latest manual startup corrections.
Previously saved integration bounds remain unchanged until edited or reselected.
There is no active motion alignment, noise removal, averaging, shared-grid
interpolation, or mean-data export. `plotData.m` has been retired.

## Standard group outputs

The default group folder is `results/{Re}{direction}/`, for example `results/1000f/`:

| File | Contents |
| --- | --- |
| `angular_displacement.png` | Python analysis: one audio-aligned line per successful saved trial |
| `{trial}_signal.csv` | Original timestamps and dimensional torque: `time_seconds,torque_dyn_cm` |
| `{trial}_jumps.csv` | Separate signed startup/shutdown results, predictions and window sensitivity |
| `jump_summary.csv` | All successful individual event rows with trial identity |
| `torque.png` | Overlaid reconstructed torque with each trial's selected integration windows |
| `cumulative_integral.png` | Overlaid raw cumulative torque integrals with each trial's selected integration windows |
| `cumulative_jump_comparison.png` | Main cumulative jump comparison: startup/shutdown cumulative torque against theory and their sum relative to zero; point labels show relative errors and magnitude mismatch; thin ranges show window sensitivity |
| `relative_errors.csv` | The same three percentages at full precision; undefined values are blank |
| `{trial}_jump_diagnostics.png` | Per-trial startup/shutdown torque close-ups and local cumulative torque against theory, saved at 300 dpi |

The torque and cumulative-integral plots shade every trial's startup/shutdown
windows on the same audio-aligned axis as its curve. Legends identify the trial
and event.
The per-trial diagnostic figure uses original recording seconds and is saved
even with `--no-show`; its filename is independent of figure numbering.

MATLAB generation leaves group analysis products intact. Python analysis rebuilds
the requested group's products using every saved trial, reports failures,
and removes stale analysis products for failed trials. Results retain their
existing names and layout under `results/{group}/`.

Generated data uses `paths.data/{group}/{trial}_data.csv` and
`{trial}_metadata.json`. Optional masked/overlay MP4s and overlay-frame folders
use that same group directory. Extracted frames use
`paths.frames/{group}/{trial}_frames/`, including when `paths.frames` is on an
external drive. Legacy flat paths remain readable; ambiguous duplicate flat and
grouped artifacts cause a collision error. CSV columns and units are unchanged.

## Existing data and migration

Python discovers saved `{trial}_data.csv` / `{trial}_metadata.json` pairs in
both flat and grouped folders. Ingestion ranges only select MATLAB preprocessing;
Python analysis includes all saved trials, including nonconsecutive IDs.

For example, if `1000r1`–`1000r4` already exist, set `"1000r": [5, 9]` and run
extraction and generation for `1000r5`–`1000r9`. `python main.py 1000 r` then
analyzes all nine. Setting `"1000r": []` skips new ingestion and keeps all nine
available. Original recordings are unnecessary for Python analysis and plots.

To organize existing generated files, preview and then apply the migration:

```bash
python main.py --migrate-data                         # dry run; data only
python main.py --migrate-data --apply                 # move generated data
python main.py --migrate-data --include-frames         # preview external frame moves
python main.py --migrate-data --include-frames --apply
```

The command preflights destination conflicts, never overwrites destinations,
keeps CSV/metadata bytes unchanged, and safely reruns. Frame migration requires
`--include-frames`; no external-drive migration has been performed. Unrecognized
legacy names and original MTS recordings are left untouched. Stop other writers
during migration. Flat paths remain readable until migration is completed.

To rebuild plots with saved selections and no connected recordings:

```bash
python main.py 1000 r --non-interactive --no-show
python main.py 1500 r --non-interactive --no-show
python plot_reynolds_summary.py --direction r --no-show
```

## Global integral plots versus Reynolds number

Use the same `experiment_config.json` as generation and analysis. Global plots
select all saved IDs for the requested direction, regardless of empty or
new ingestion ranges. Generate group summaries first, then run:

```bash
python plot_reynolds_summary.py --direction r
python plot_reynolds_summary.py --direction r --no-show
python plot_reynolds_summary.py --config path/to/config.json --direction f --no-show
```

The command uses the same default configuration and config-relative path
resolution as `main.load_config`. It reads saved groups' existing
`jump_summary.csv` files and selects saved trial IDs for the requested
direction. It does not rerun fitting, open selection pickers, update metadata or
rewrite per-group results. **Regenerate summaries with `main.py` after changing
analysis inputs, physical constants or saved selections**; this plotting command
cannot detect stale summaries from those changes.

Five separate figures are saved at 300 dpi under
`{configured output path}/global/{direction}/` (default `results/global/r/`):

| File | Plotted value (dyn·cm·s) |
| --- | --- |
| `startup_vs_reynolds.png` | Startup `J_measured_dyn_cm_s` |
| `shutdown_vs_reynolds.png` | Shutdown `J_measured_dyn_cm_s` |
| `mismatch_vs_reynolds.png` | **Startup + shutdown imbalance:** `J_start + J_stop` |
| `startup_shutdown_signed_vs_reynolds.png` | Startup and shutdown signed integrals together, with a mean for each event/group |
| `startup_shutdown_magnitude_vs_reynolds.png` | Individual absolute integrals and the mean of absolute integrals for each event/group |

The imbalance is the signed sum of the two event integrals, matching the lower
panel in `main.plot_jump_comparison`. It is neither the percentage magnitude
mismatch nor an integral between the events. Reverse-flow startup remains
negative and shutdown positive in the signed views; directions are never pooled.
Only the combined magnitude view applies absolute values for visualization.

In the original three figures, each numeric Reynolds-number cluster shows
individual trials, their unweighted arithmetic mean as a larger orange dot,
and a distinct marker for the stored
theoretical prediction. All markers are partially transparent to reveal overlap.
Only means receive percentage labels. Startup and shutdown use relative error,
`100 * abs(mean - theory) / abs(theory)`, rather than averaging individual trial
errors. These labels show `Error N/A` for missing, inconsistent or zero theory.
The imbalance plot instead labels the **magnitude mismatch of the mean event
integrals**, `100 * abs(abs(mean_start) - abs(mean_stop)) /
max(abs(mean_start), abs(mean_stop))`. This uses the same definition as the
per-trial `Mismatch` labels in `cumulative_jump_comparison.png`, applied to the
group's mean startup and shutdown, rather than averaging the trial percentages.
It is independent of theory, and is undefined only when both mean event integrals
are zero (`Mismatch N/A`). The imbalance y-values remain signed `J_start + J_stop`
in dyn·cm·s; the percentage is only a label.
The imbalance theory uses paired stored startup-plus-shutdown predictions
(zero for the implemented model). Missing or inconsistent theory is reported
and its marker omitted; valid measurements and means remain. Available theory
must agree within a group (roundoff tolerance `rtol=1e-10`, `atol=1e-12`), and
predictions are never averaged. There are no uncertainty ranges, error bars,
regression lines or connections between trial points in these global figures.

The two combined views use filled blue circles for startup and filled orange circles
for shutdown, with larger mean markers. Partial transparency and distinct colors
help reveal overlapping events. All trials remain at their actual Reynolds
number. The magnitude mean is `mean(abs(J))` over the selected trials, which can
differ from `abs(mean(J))`; the underlying signed measurements are unchanged.
These views omit percentage annotations to keep the event series readable.

In the combined signed plot, startup and shutdown each have a continuous straight
theory line through the origin. The magnitude plot has one shared theoretical
magnitude line. Slopes come from the validated stored predictions `J_theory / Re`,
with the same roundoff tolerance as group predictions. Available slopes must
agree across Reynolds groups, and startup/shutdown slopes must be equal and
opposite. No measured-data regression or hardcoded theory coefficient is used.
Lines span zero to the upper displayed Reynolds limit, including single-group
plots. Missing theory omits its event line if no other group supplies a valid
prediction; consistent available groups can supply the model slope across the
displayed range. A shared magnitude line can use either available event.
Inconsistent predictions suppress the affected event line and the shared
magnitude line; failures of the equal-and-opposite check suppress both event lines.
Diagnostics report these conflicts while retaining all valid measurements and
means. Validation for these new lines does not alter the original three plots'
stored-theory markers or annotations.

Events are paired by exact `run_identity` and event name, independent of row
order. Malformed identities, wrong directions, nonfinite measurements, duplicate
or incomplete event pairs and deselected trials are reported and excluded;
missing groups are reported and skipped. All five figures use the same valid
trial pairs. Stored density, geometry, rotation sign and flow conversion
(`flow_rate_cm3_s / Re` and `Q_signed_cm3_s / Re`) must be comparable across
included trials/groups, with signed flow consistent with direction. Optional
flow/geometry values may be uniformly unavailable; mixed known/unknown settings
are not assumed comparable. Conflicting physical settings are reported and
prevent combined figures. If no comparable valid data remain, only this
command's five old figures for the requested direction are removed, and the
command exits with status 1.

## Numerical conventions

`main.py` retains Rachel Bertaud's introductory explanation and numbered analysis
sections. The per-trial system-parameter report also retains its original formulas
for inertia, damping, natural period and decay timescale.

Initial measurement noise remains in the input. FFT padding, period/frequency
conventions, filtering, reconstruction and the forward solver retain their existing
behavior. Each recording must contain at least two finite samples with increasing,
approximately uniform timestamps (`rtol=1e-4`, `atol=1e-12` s). Irregular or missing
samples invalidate this workflow; repeat the experiment rather than filling gaps.
The FFT length is twice the smallest power of two at least as large as the input,
with a minimum base length of 2048. Constant initial-value padding precedes the
full recording. Exports contain every original sample, without padding.

`I = kappa / c**2` scales reconstructed forcing to torque in dyn·cm.
`impulse_funcs.py` integrates the raw signed torque separately around startup and
shutdown with exact, linearly interpolated integration endpoints, preserving its
existing units (**dyn·cm·s = g·cm²/s**), signs and sensitivity calculations. There
is no baseline subtraction and no integral between events. Window sensitivity varies
both ends by 0.25 s where valid, including the original window in the reported
range; it is not a statistical confidence interval.
The terminal reports `[contracted, expanded]` in that order, rather than sorting
those two cumulative torque values. Invalid variants appear as `unavailable`. CSVs include
the corresponding `sensitivity_contracted_dyn_cm_s` and
`sensitivity_expanded_dyn_cm_s` fields. The numerical minimum/maximum across
valid variants and the original window remain available for plotting the range.

The terminal compares each measured cumulative torque with theory using
`relative error = 100 * abs(measured - theory) / abs(theory)` percent, displayed
to two decimal places. A zero theoretical cumulative torque gives an undefined percentage.
Residuals and measured/theory ratios remain available in the CSV exports.

The current flow setting preserves the existing relation `Q = Re * 0.0094` cm³/s,
a positive **total flow through both arms**. Verify this conversion for your batch.
`G` is the signed integral `x dy - y dx` along **one arm outward from the hub**;
`sign_convention` maps measured positive rotation to the positive z direction.
With signed outward-positive flow, predictions are `J_start = -rho * Q_signed * G`
and `J_stop = +rho * Q_signed * G`, multiplied by `sign_convention`. There is no
extra factor of two. Set the flow conversion or geometry factor to `null` to leave
theory unavailable. Forward-flow comparisons remain exploratory because the
arm-only model may omit external-jet contributions. Damping-factor investigation
and other numerical changes are deferred.

## Verification

The original test suite and focused batch checks exercise reconstruction/export,
group selection, saved timing/selections, reselection and invalid-trial exclusion.
Run `python -m unittest discover -s tests -v` to include the synthetic global-plot
checks for numeric Reynolds ordering, exact trial/mean/theory series, event
pairing, exclusions, missing theory/groups, physical comparability and direction
separation. Plot-artist checks also verify that no uncertainty ranges appear.
The combined-plot checks also verify exact signed and magnitude series, means of
individual magnitudes, theory slopes and symmetry, missing/inconsistent theory,
single-group inputs, all five outputs and byte-identical original PNGs.
The current 64-test suite passes in the default Python 3.9 environment. The
earlier 42-test suite passed with both Python 3.12/NumPy 2 and the default
Python 3.9/NumPy 1.23.5 environment, including shared group/repeat expansion
and invalid configuration checks. MATLAB R2026a syntax
checks and synthetic AVI/JPEG/WAV checks also passed: extraction completion and
reuse, stale-frame cleanup, reference preservation, independent audio candidates
and offsets, metadata preservation, missing-detection exclusion and continuation.
For the existing `1000r1_data.csv`, the pre-refactor and refactored calculations
were compared with motion preprocessing disabled, fit start 25 s and fixed windows
`[11.149896, 12.149896]` and `[21.149896, 22.149896]`. Fitted constants, every torque
sample, forward reconstruction and both signed cumulative torque values were identical. The
cumulative torque values were -247.57673961727286 and 276.5888171556742 dyn·cm·s. These fixed windows
are a numerical regression check, not confirmation of the audio event times.

For the storage/range change, all six existing trials (`1000r1`–`1000r4`,
`1500r1`–`1500r2`) were recomputed from isolated copies before and after the
change. Fitted parameters, torque arrays, cumulative integrals, event integrals,
theoretical predictions, forward reconstruction and reconstruction errors were
exactly equal. Torque arrays and saved event-result fields also matched the
existing result CSVs exactly. Original data, metadata and result CSVs remained
byte-identical. These checks validate numerical preservation for the available
CSV inputs, including their existing saved selections.

Actual MTS recordings and a saved ruler calibration are absent from this checkout,
so real video decoding, laser tracking, audio cue quality, interactive calibration
and playback require verification with the experiment recordings.

## Mathematical Details of Method 

In the experiments, the sprinkler undergoes forcing by the syringe pump for ~10 seconds, then is left to ring without any forcing. During the ringing with no forcing, we assume the sprinkler's motion can be modelled as a damped harmonic oscillator,

$$
    \ddot{\phi} + 2\gamma\dot{\phi} + \omega^2 \phi = 0
$$

with damping coefficient $\gamma$ and natural frequency $\omega$. To extract the torque using Fourier transforms, we need to know what $\gamma$ and $\omega$ are. We go about finding them two ways:
1. Using simple estimates from data manipuation (NO FIT)
   
   <img width="642" height="402" alt="Screenshot 2026-06-03 141036" src="https://github.com/user-attachments/assets/5f5ccc62-f70b-4177-8661-80cf5a65668f" />
      
   Here, the blue, solid curve is spring response of sprinkler. Orange, circular data points are where peaks occur (change in direction of curve).​ This excellent data shown         was provided by Kelly and Jesse, Dec. 2025​.

   <img width="284" height="305.5" alt="Screenshot 2026-06-03 141212" src="https://github.com/user-attachments/assets/83be9d5d-7133-4ec4-9390-364e1e822899" />

   ### Finding an estimate for $\gamma$
   The black curve in the image above will fit the form of $y(t) = \omega(t_0) e ^{-\gamma t}$. Our goal here will be to fit the data to $y(t)$ to retrieve $\gamma$. We can take all the positive peaks that lie on that black curve, normalize them by dividing by $\omega(t_0)$, then shift the curve to $t = 0$. In this modified for, the data is simple to fit for $\gamma$ using standard scipy/numpy functions.

   
   ### Finding an estimate for $\omega$
   The period of this data is the time between every other peak. Thus, for all $t_i > t_0$, we define
   $$
        b_i = t_{i + 2} - t_i.
   $$
   It follows that the estimate for $\Omega$ is equal to $\frac{2\pi}{\text{mean}(b_i)}$ and $\omega = \sqrt{\Omega^2 + \gamma^2}$.
      
2. Fitting the experimental data to the ODE form using non-linear least squares ([scipy.curve_fit](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html)) (FIT)

   The code refers to this as the fitting step. After step 1 ( which is performed no matter the fit flag ), we have good estimates for $\omega, \gamma$. We take these estimates, set them as initial conditions for the parameters of interest in scipy.curve, and fit the data to the ODE form using non-linear least squares. The output of this curve fit is then used for $\gamma, \omega$ for the rest of the code.

\noindent Once we have values for $\gamma, \omega$, we can use a Fourier transform to solve 

$$
    \ddot{\phi} + 2\gamma\dot{\phi} + \omega^2 \phi = \tau(t).
$$

 We rewrite our ODE as a convolution in the time domain as
 
$$
\phi(t) = \int G(t - s)\tau(s)ds,
$$

where $G$ is the Green's function of the system. In Fourier space, this convolution becomes multiplication, or

$$
\hat{\phi} = \hat{G} \cdot \hat{\tau}. 
$$

\noindent We take the Fourier transform of our ODE directly, where $k$ is the frequency variable, to solve for $\hat{G}$.

$$
\begin{aligned}
-k^2\hat{\phi} - 2i\gamma k \hat{\phi} + \omega^2 \hat{\phi} &= \hat{\tau}\\
(-k^2 - 2i\gamma k + \omega^2)\hat{\phi} &= \hat{\tau} \\
\hat{\phi} &= \left( \frac{-1}{k^2 + 2i\gamma k - \omega^2} \right) \hat{\tau} \\
&= \hat{G} \cdot \hat{\tau}
\end{aligned}
$$

 Therefore, 
 
$$
\hat{G} = \frac{-1}{k^2 + 2i\gamma k - \omega^2}.
$$

Since we know the Green's function and $\phi$, we can solve for the torque signal of the system as:

$$
\begin{aligned}
\hat{\tau} &= \frac{\hat{\phi}}{\hat{G}}\\
\hat{\tau} &= \frac{\text{fft}(\phi)}{\hat{G}} \\
\text{ifft}(\hat{\tau}) &= \text{ifft} \left( \frac{\text{fft}(\phi)}{\hat{G}}\right)  \\
\tau &= \text{ifft} \left( \frac{\text{fft}(\phi)}{\hat{G}}\right) 
\end{aligned}
$$

Now, the question becomes how well this torque signal represents our data. We can now take $\tau$ and solve the ODE form forward for $\phi$. Thus, we are generating analytical data using the extracted $\tau$ from the experimental data. We then calculate the norm 2 error between the analytical data and the experimental data as a reference for the accuracy of our method.
