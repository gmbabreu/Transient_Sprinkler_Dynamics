"""Scientific stages of the transient-sprinkler analysis, in execution order.

Time coordinates are original video seconds. Angles are radians after loading;
physical quantities use cgs units. These two functions only calculate results:
main.py owns files, selections, validation, reports, plots and exports.
"""

import numpy as np
from scipy import integrate
from scipy.signal import cheby2, sosfiltfilt

from Main_Functions.dataread_funcs import fit_segments
from Main_Functions.estimate_funcs import est_a, est_b, get_constants, fit_phi
from Main_Functions.process_funcs import combine_data
from Main_Functions.fft_funcs import torque_solver, phi_from_torque
from Main_Functions.impulse_funcs import analyze_impulses


def reconstruct_torque(full_t, full_phi, t_peaks, phi_peaks, t_fit, settings, physics):
    """Fit the free decay, reconstruct torque, and verify the forward response.

    Inputs are the loaded/validated recording, its peaks, selected fit start [s],
    and existing analysis/physics settings. Return measured-grid arrays, fitted
    parameters and physical parameters for main.py's reporting and selections.
    Numerical stages retain the original order, including the forward solve
    before dimensional scaling. No measurement-noise correction is applied.
    """
    # SECTION TWO - ESTIMATE AND FIT THE FREE-DECAY OSCILLATOR
    # Assume the selected tail is unforced: phi'' + 2*b*phi' + c^2*phi = 0.
    # phi = exp(-b*(t-t0)) * (C_1*cos(a*(t-t0)) + C_2*sin(a*(t-t0))),
    # a = sqrt(c^2-b^2). Time is in s, phi/C_1/C_2 in rad, b/c in s^-1.
    # Keep nearest-peak segmentation, peak estimates and optional fit unchanged.
    fit_index, fit_peaks_index, t_fit, t_seg, phi_seg = fit_segments(
        full_t, full_phi, t_peaks, t_fit)

    # Estimate a from data.
    a_est = est_a(fit_peaks_index, t_fit, t_peaks)
    # Estimate b from data.
    b_est = est_b(fit_peaks_index, t_peaks, phi_peaks)
    # Estimate c from a and b.
    c_est = np.sqrt((a_est*a_est) + (b_est*b_est))
    # Get the ODE constants using the estimates for gamma and omega.
    C1_est, C2_est = get_constants(b_est, c_est, full_phi[fit_index])
    if settings["fit_enabled"]:
        # Fit the data.
        b, c, C_1, C_2 = fit_phi(t_seg, phi_seg, b_est, c_est, C1_est, C2_est,
                                full_t[fit_index])
    else:
        # Set end values as estimates if no fit is requested.
        b = b_est
        c = c_est
        C_1 = C1_est
        C_2 = C2_est

    # Define the analytical solution in terms of the newly found values.
    def phi_an(t):
        """Retain the analytical fit accepted by the existing padding API.

        combine_data currently prepends constant initial-value padding; it
        does not append this analytical curve or replace measured samples.
        """
        t0 = t_peaks[fit_peaks_index]
        a = np.sqrt(c**2 - b**2)
        return np.exp(-b * (t - t0)) * (C_1 * np.cos(a * (t - t0))
                                        + C_2 * np.sin(a * (t - t0)))

    # SECTION THREE - PREPARE THE SIGNAL AND PAD FOR THE FFT
    # Uniformly sampled phi(t) [rad] is retained, including initial noise:
    # preprocessing/remove_noise stays disabled. Prepend phi[0] at the same dt
    # until N_f = 2*max(2048, next_power_of_two(N)); do not append the fit.
    # Padded time remains in original video seconds (including negative padding).
    fourier_t, fourier_phi, N_added = combine_data(full_t, full_phi, phi_an)

    # SECTION FOUR - RECONSTRUCT FORCING AND FILTER
    # For phi'' + 2*b*phi' + c^2*phi = f, f has units rad/s^2.
    # The helper uses f = IFFT(FFT(phi)/G_hat(-k)),
    # G_hat(k) = -1/(k^2 + 2*i*b*k - c^2), k spaced by 2*pi/L.
    # Preserve its Nyquist convention, L = t[-1]-t[0], real projection, and
    # optional zero-phase Chebyshev-II filter with normalized Wn (not Hz).
    N_f = len(fourier_phi)
    # Keep the existing FFT period and normalized Chebyshev cutoff convention;
    # batching changes which trials run, not the numerical reconstruction.
    L = fourier_t[-1] - fourier_t[0]
    # Output the extracted torque signal, then retain its real part.
    extracted_torque = torque_solver(N_f, b, c, L, fourier_phi)
    extracted_torque = np.real(extracted_torque)
    if settings["lowpass_enabled"]:
        sos = cheby2(N=settings["cheby2_order"], rs=settings["cheby2_rs"],
                     Wn=settings["cheby2_wn"], btype="low", output="sos")
        extracted_torque = sosfiltfilt(sos, extracted_torque)

    # Remove exactly the prepended padding. These timestamps then match the
    # input recording, so exports and integration use only measured samples.
    fourier_t_seg = fourier_t[N_added:]
    fourier_phi_seg = fourier_phi[N_added:]
    extracted_torque_seg = extracted_torque[N_added:]

    # SECTION FIVE - SOLVE THE FORWARD OSCILLATOR AND CHECK RECONSTRUCTION
    # Solve phi'' + 2*b*phi' + c^2*phi = f with phi=phi'=0 initially.
    # The existing helper uses PCHIP forcing and RK45 with its default tolerances.
    # Relative error = sqrt(integral((phi_forward-phi)^2 dt)) /
    #                  sqrt(integral(phi^2 dt)); both curves are in rad.
    # Keep this solve before physical scaling, as in the original execution order.
    phi_forward = phi_from_torque(N_f, fourier_t_seg, extracted_torque_seg, b, c)
    # SciPy's trapezoidal rule also supports environments using NumPy 1.x.
    error_forward = (np.sqrt(integrate.trapezoid((phi_forward - fourier_phi_seg)**2,
                                        x=fourier_t_seg))
                     / np.sqrt(integrate.trapezoid(fourier_phi_seg**2, x=fourier_t_seg)))

    # SECTION SIX - PHYSICAL PARAMETERS AND DIMENSIONAL TORQUE
    # I = kappa/c^2 [g cm^2], torque = I*f [dyn cm], kappa [dyn cm/rad].
    # EXISTING CONVENTION: delta = b*I, although the fitted/forward ODE uses
    # 2*b*phi'. Retain this damping-factor inconsistency and all derived values.
    # xi = delta/(2*sqrt(kappa*I)); w0 [s^-1], T0 and Td [s].
    kappa = physics["kappa_dyn_cm_per_rad"]
    # The solver returns forcing per unit inertia. Scale by I to obtain
    # dimensional torque in dyn cm, using the same relation as before.
    I = kappa/c**2  # solve for inertia, in units of g cm^2
    delta = b*I  # solve for damping coefficient, in units of g cm^2 s^-1
    xi = delta/(2*np.sqrt(kappa*I))  # solve for damping ratio, dimensionless
    w0 = (delta/I)/(2*xi)  # solve for natural frequency, units of s^-1
    T0 = (2*np.pi)/w0  # solve for natural, undamped period, units of s
    Td = 1/(xi*w0)  # solve for decay timescale, units of s

    true_torque = extracted_torque_seg * I
    # This full-record cumulative curve is for plotting. The two measured
    # event impulses below are integrated separately within their own bounds.
    cummInt = integrate.cumulative_trapezoid(np.real(true_torque), fourier_t_seg, initial=0)

    return dict(time=fourier_t_seg, torque=true_torque, cumulative=cummInt,
                phi_forward=phi_forward, error_forward=error_forward,
                b=b, c=c, kappa=kappa, inertia=I, delta=delta, xi=xi,
                w0=w0, T0=T0, Td=Td)


def measure_event_impulses(time, torque, full_t, windows, physics, reynolds, direction):
    """Finish section six after main.py selects startup/shutdown windows [s].

    J = integral(torque dt) [dyn cm s] for each separate event, with linear
    interpolation at exact endpoints and trapezoidal quadrature. Use raw torque:
    no baseline subtraction and no integral across the interval between events.

    Existing theory: J_start = -s*rho*Q*G, J_stop = +s*rho*Q*G [dyn cm s].
    Q = Re*flow_factor [cm^3/s] is total flow through both arms, signed positive
    for f and negative for r; rho [g/cm^3], signed one-arm G [cm^2], s = +/-1.
    Keep these implemented signs/conventions without reconciling other documents.

    The same helper calculates both impulses, theory and window sensitivity in
    its existing order. Sensitivity contracts then expands each end by 0.25 s;
    invalid variants are excluded. No integration endpoints are clamped.
    """
    flow_factor = physics["flow_rate_per_re_cm3_s"]
    return analyze_impulses(
        time, torque, full_t,
        startup_window=windows[0], shutdown_window=windows[1],
        rho_g_cm3=physics["rho_g_cm3"],
        flow_rate_cm3_s=None if flow_factor is None else reynolds * flow_factor,
        geometry_factor_cm2=physics["geometry_factor_cm2"],
        sign_convention=physics["sign_convention"], spin_dir=direction,
    )
