import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

def torque_solver(N_f, b, c, L, fourier_phi):
    def G_hat(k):
        return -1.0 / (k**2 + 2.0*1j*b*k - c**2)
    
    kk = (2*np.pi / L) * np.concatenate([
        np.arange(0, N_f//2),
        [0],
        np.arange(-N_f//2 + 1, 0)
    ])
    
    extracted_torque = np.fft.ifft(np.fft.fft(fourier_phi) / G_hat(-kk))
    
    return extracted_torque

def phi_from_torque(N_f, fourier_t, extracted_torque, b, c):
    extracted_torque_seg = np.real(extracted_torque)  # signal from 0 to end of full_t

    t_seg = fourier_t
    
    # interpolant for the forcing
    para = PchipInterpolator(t_seg, extracted_torque_seg, extrapolate=False)

    def dydt(t, y):
        p = para(t) if t_seg[0] <= t <= t_seg[-1] else 0.0
        return [y[1], -(2*b)*y[1] - (c**2)*y[0] + p]

    y0 = [0.0, 0.0]
    sol = solve_ivp(dydt, [t_seg[0], t_seg[-1]], y0, t_eval=t_seg, method='RK45')

    phi_gen = sol.y[0]
    return phi_gen
