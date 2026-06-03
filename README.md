# Transient Dynamics Sprinkler — `main.py`

Written by Rachel Bertaud, Colorado School of Mines, 2026.

---

## What This Code Does

Extracts the **torque signal** of a sprinkler system from measured angular position data. It works by assuming the free-decay tail of the sprinkler's angular position fits a **damped harmonic oscillator ODE**, then uses that fit to identify the system's natural frequency (ω) and damping coefficient (γ). With those parameters known, a Fourier-based method reconstructs the full torque signal. A forward ODE solve is then run to verify the result. Refer to the end for the matematical breakdown.

---

## System Overview

```
Step 0 → Set user-defined inputs (switches, data selection)
Step 1 → Read CSV data, interactively pick fit window
Step 2 → Estimate ω (natural frequency) and γ (damping coefficient)
Step 3 → Clean noise from data (optional)
Step 4 → Extract torque signal via Fourier transform
Step 5 → Forward ODE solve to verify extracted torque signal
Step 6 → Save torque signal to CSV and plot results
```

---

## Dependencies

```
numpy
scipy
matplotlib
```

Install with:
```bash
pip install numpy scipy matplotlib
```

---

## Data Format

Input data is a `.csv` file with **two columns** (header row is skipped):

| Column 0 | Column 1 |
|----------|----------|
| time (t) | angular position (φ) |

**File naming convention:** `{direction}_{Re}_trial{N}.csv`

- `direction`: `forward` or `rev`
- `Re`: Reynolds number (e.g. `500`, `1000`)
- `N`: trial number (e.g. `1`, `2`)

Example: `forward_1000_trial1.csv`

---

## User-Defined Inputs (Step 0)

At the top of `main.py`, set these variables before running:

| Variable | Type | Description |
|---|---|---|
| `plot_switch` | `0` or `1` | `1` = show all plots at end, `0` = suppress all plots |
| `write_t` | `0` or `1` | `1` = save extracted torque signal to CSV, `0` = skip saving |
| `fit_switch` | `0` or `1` | `1` = refine ω/γ estimates with a curve fit, `0` = use raw peak estimates |
| `proc_data_switch` | `0` or `1` | `1` = clean noise from data, `0` = use raw data |
| `spin_dir` | string | Spin direction — `"forward"` or `"rev"` |
| `re` | string | Reynolds number matching the filename of interest — e.g. `"1000"` |
| `trial` | int | Trial number matching the filename of interest — e.g. `1` |
| `data_dir` | string | **Full local path** to the directory where your CSV data files live |

---

## How to Run

```bash
python main.py
```

The script is **interactive** — it will open one plot window with a slider. Drag the slider to select the **start time** for the free-decay fit (should be near a peak a bit after the sprinkler begins ringing down). Close the window when done. After the window is closed, the script runs automatically and prints estimated/fitted values to the terminal.

---

## Outputs

**Terminal output:**
- Estimated γ, ω, and ODE constants (c1, c2) from peak analysis
- Fitted γ, ω, c1, c2 after curve fitting (if `fit_switch = 1`)
- Error between forward solve for $\phi$ us9ing extracted torquer and experimental $\phi$ data.
- Path of saved signal file (if `write_t = 1`)

**Saved file (if `write_t = 1`):**

Output is written to `data_dir` with the name `{spin_dir}_{Re}_trial{N}_signal.csv` — e.g. `forward_1000_trial1_signal.csv`. If a file with that name already exists, it is deleted and replaced. The file contains two columns covering the second half of the stitched signal (the physically meaningful forcing window):

| Column | Description |
|--------|-------------|
| `t` | time |
| `torque_signal` | real part of extracted torque signal |

**Plots (if `plot_switch = 1`):**
1. **Analytical solution** — raw data (scatter) overlaid with the damped oscillator fit (line); zoomed to the fit region, with the fit start point highlighted
2. **Franken signal** — the stitched signal used for the FFT, with color-coded segments: prepended constant (blue) and original data (pink)
3. **Extracted torque signal** — real part of the torque signal over the forcing window
4. **ODE verification** — forward-solve for φ(t) (line) overlaid on measured data (scatter). forward solve uses extracted torque signal.

---

## File Structure

```
Transient_Dynamics_Sprinkler/
├── main.py                          # Main code — run this
└── Main_Functions/
    ├── dataread_funcs.py            # CSV reading, peak finding, interactive slider plots
    ├── estimate_funcs.py            # ω and γ estimation and curve fitting
    ├── process_funcs.py             # Noise removal and signal stitching ("franken" signal)
    ├── fft_funcs.py                 # Fourier-based torque extraction and forward ODE solve
    └── plot_funcs.py                # All plotting functions
```

---

## Method Summary

1. **ω estimation:** Period is measured from peak-to-peak spacing in the free-decay tail. `ω_d = 2π / T_mean`, then `ω = sqrt(ω_d² + γ²)`.

2. **γ estimation:** Peak amplitudes decay as `e^(-γt)`. An exponential is fit to normalized peak amplitudes to extract γ.

3. **ODE fit (optional):** Both ω and γ are further refined by fitting the full free-decay segment to the analytical damped oscillator solution using `scipy.optimize.curve_fit`.

4. **Signal stitching ("Franken" signal):** A constant pre-forcing segment is prepended to the data. If the data is being processes, this includes the analytical free-decay solution is appended after the forcing ends. This creates a periodic-compatible signal for the FFT.

5. **Torque extraction (FFT):** The torque signal `f(t)` is extracted by deconvolving the ODE's transfer function `G(k)` in Fourier space:
   `F(k) = Y(k) / G(k)`, where `G(k) = -1 / (k² + 2iγk - ω²)`

6. **Verification (forward ODE):** The extracted torque signal is used as forcing when solving ODE (RK45) to reconstruct φ(t).

---

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

