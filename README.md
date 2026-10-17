# Transient Sprinkler Dynamics

Work by Rachel Bertaud, Gabriel Abreu, and Kelly Yu.

Reconstruct torque from each trial's angular displacement by fitting its free-decay
tail, applying the existing Fourier solver and torque filter, and checking a
forward ODE reconstruction. Measure two separate signed impulses around pump
startup and shutdown, and compare individual trials against theory.

## Four-step workflow

1. **Configure the batch.** Edit [`experiment_config.json`](experiment_config.json).
   Put recordings in the configured `paths.trials` folder. List the run groups
   under `trial_groups` and set one shared `trials_per_group` count. For example,
   `"trial_groups": ["1000f", "500r", "1500f"]` with `"trials_per_group": 3`
   processes repeats 1, 2 and 3 for every group: `1000f1` through `1000f3`,
   `500r1` through `500r3`, and `1500f1` through `1500f3`. The default count of 1
   retains the four existing example trial names. Keep only groups from one
   compatible camera/ruler setup. Names follow `{Re}{direction}{trial}`, where
   `f` is forward/outward flow and `r` is reverse/inward flow. Paths are relative
   to the JSON file, independent of the working directory. Python also accepts
   `--config /path/to/experiment_config.json`. Both loaders expand the groups
   and repeat count into full trial IDs in memory; individual CSV and metadata
   filenames stay unchanged.
2. **Extract and calibrate in MATLAB.** Run `Data_Generation/save_frames.m`, then
   `Data_Generation/calibration/calibrateRulerAxis.m` to click the ruler marks.
   Complete extractions are skipped. Set `video.regenerate_frames` to `true` to
   replace generated frames, then restore it to `false`. A `.complete` marker
   records frame count and frame rate only after extraction succeeds; incomplete
   folders are rebuilt and old generated JPEGs are removed first, so shorter
   replacement recordings cannot leave stale trailing frames. The first
   group's repeat 1 supplies the calibration reference only if it is absent.
   Ordinary reruns, including frame regeneration, preserve that image. For a new
   setup, choose new reference/calibration paths or deliberately remove both old
   files and recalibrate. The ruler projection method is unchanged and its saved
   result is reused for the entire batch.
3. **Generate individual data in MATLAB.** Run `Data_Generation/generateData.m`.
   It reuses extracted frames, tracks the weighted laser centroid, projects onto
   the ruler, and keeps the original first-position angle zero. It writes
   `data/{trial}_data.csv` with `time_seconds,angle_degrees`, starting at time zero,
   plus one `{trial}_metadata.json`. It independently searches startup and shutdown
   audio intervals for the loudest sample across all channels. These are candidate
   cues; inspect and correct them when necessary. Missing laser detections or
   nonfinite angles invalidate the entire trial: affected frame numbers are
   reported and saved, the old angle CSV is removed, and other trials continue.
   No missing positions are interpolated. One angular-displacement comparison
   figure per Reynolds number and direction is saved from valid trials.
4. **Analyze a selected group in Python.** For example:

   ```bash
   python main.py 1000 f
   ```

   The original form, `python main.py --re 1000 --direction f`, also works.
   Both forms accept the same additional options below.

   Every matching configured trial is analyzed separately. On its first run,
   select the free-decay fit start with the existing slider and choose two torque
   endpoints around each event. Close the picker or press Enter to accept;
   incomplete pairs retain the supplied windows. Fit start and both windows are
   saved in the same per-trial metadata and reused on subsequent runs. Add
   `--reselect` to select them again. `--no-show` suppresses final diagnostic
   displays while still allowing required selection; `--non-interactive` requires
   a saved fit start and uses saved windows or audio-initialized defaults.
   Failed trials are reported while valid trials continue. A run with any failure
   exits with status 1.

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
| `trial_groups`, `trials_per_group` | Groups such as `1000f` and one positive integer repeat count shared by every group |
| `video` | Extension, mirror distance in inches, laser YCbCr limits, extraction regeneration and tracking diagnostics |
| `physics` | Spring stiffness, density, flow conversion, signed one-arm geometry factor and sign convention |
| `analysis` | Existing fit switch, Chebyshev torque-filter settings and diagnostic display |
| `audio` | Separate search intervals in **audio seconds**, audio-to-video offset, inspection and playback |
| `integration_window_offsets_s` | Initial bounds relative to each effective event time in **video seconds** |

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

Each `data/{trial}_metadata.json` contains:

- `valid`, `failure_reason`, `invalid_frames`: tracking status and 1-based frame numbers.
- `audio_candidates_s`: independently estimated `startup` and `shutdown` times.
- `manual_event_times_s`: optional `startup` and `shutdown` corrections in original video seconds.
- `fit_start_s`: selected free-decay start in original seconds.
- `integration_windows_s`: selected `startup` and `shutdown` pairs in original seconds.

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

Comparison curves use `aligned_time = original_time - effective_startup_time`
on each trial's own samples. Saved data and integration selections retain original
timestamps. Changing an integration window never changes audio timing or plot
alignment. After editing manual times, rerun generation to refresh its angular
comparison plot; Python comparison plots refresh on the next analysis run.
Previously saved integration bounds remain unchanged until edited or reselected.
There is no active motion alignment, noise removal, averaging, shared-grid
interpolation, or mean-data export. `plotData.m` has been retired.

## Standard group outputs

The default group folder is `results/{Re}{direction}/`, for example `results/1000f/`:

| File | Contents |
| --- | --- |
| `angular_displacement.png` | MATLAB generation: one line per valid trial, audio-aligned |
| `{trial}_signal.csv` | Original timestamps and dimensional torque: `time_seconds,torque_dyn_cm` |
| `{trial}_impulses.csv` | Separate signed startup/shutdown results, predictions and window sensitivity |
| `impulse_summary.csv` | All successful individual impulse rows with trial identity |
| `torque.png` | Overlaid reconstructed torque with each trial's selected integration windows |
| `cumulative_integral.png` | Overlaid raw cumulative torque integrals with each trial's selected integration windows |
| `impulse_comparison.png` | Each trial's startup/shutdown impulses against theory; lines show window sensitivity |
| `{trial}_impulse_diagnostics.png` | Per-trial startup/shutdown torque close-ups and local cumulative impulses against theory, saved at 300 dpi |

The torque and cumulative-integral plots shade every trial's startup/shutdown
windows on the same audio-aligned axis as its curve. Legends identify the trial
and event; a single-trial plot also includes numerical impulse annotations.
The per-trial diagnostic figure uses original recording seconds and is saved
even with `--no-show`. For a single trial with saved selections, this is the
figure displayed as Figure 5; its filename is independent of figure numbering.

Generation removes the affected groups' previous generated results before
rebuilding. Analysis replaces its group products and removes stale per-trial
analysis CSVs, including those belonging to failed or deselected trials. It
preserves the angular figure produced by generation. Enable `video.masked_video`,
`video.overlay_video`, or `video.overlay_frames` to save masked videos, tracking
overlay videos, or individual overlay JPEGs under the configured data folder.
Older example data and impulse CSVs in `Data_Generation/data` are not selected as
results; Python requires the configured trial's valid metadata from generation.

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
is no baseline subtraction and no jump-to-jump impulse. Window sensitivity varies
both ends by 0.25 s where valid, including the original window in the reported
range; it is not a statistical confidence interval.
The terminal reports `[contracted, expanded]` in that order, rather than sorting
those two impulse values. Invalid variants appear as `unavailable`. CSVs include
the corresponding `sensitivity_contracted_dyn_cm_s` and
`sensitivity_expanded_dyn_cm_s` fields. The numerical minimum/maximum across
valid variants and the original window remain available for plotting the range.

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
All 30 Python tests passed with both Python 3.12/NumPy 2 and the default
Python 3.9/NumPy 1.23.5 environment, including shared group/repeat expansion
and invalid configuration checks. MATLAB R2026a syntax
checks and synthetic AVI/JPEG/WAV checks also passed: extraction completion and
reuse, stale-frame cleanup, reference preservation, independent audio candidates
and offsets, metadata preservation, missing-detection exclusion and continuation.
For the existing `1000r1_data.csv`, the pre-refactor and refactored calculations
were compared with motion preprocessing disabled, fit start 25 s and fixed windows
`[11.149896, 12.149896]` and `[21.149896, 22.149896]`. Fitted constants, every torque
sample, forward reconstruction and both signed impulses were identical. The
impulses were -247.57673961727286 and 276.5888171556742 dyn·cm·s. These fixed windows
are a numerical regression check, not confirmation of the audio event times.

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
