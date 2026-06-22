import numpy as np
from scipy.optimize import curve_fit

def est_omega_d(peak_index, t_insert, t_peaks):

    # determine where the nearest peak is to where we want to insert the data
    insert_near_peak = np.argmin(np.abs(t_peaks - t_insert))

    # pulls out segment of data between peak used for fit
    # and peak close to insert
    peak_seg = t_peaks[peak_index:insert_near_peak]
    

    # does difference of every other point in peak_seg
    period_est = peak_seg[2:] - peak_seg[:-2]

    #print("Mean: ", np.mean(period_est))
    
    #estimates period!
    omega_d = (2*np.pi)/np.mean(period_est)
    
    return omega_d

# def est_gamma(peak_index, t_peaks, y_peaks):
#     # picks every other peak from peak index to end
#     indices = np.arange(peak_index, len(t_peaks), 2)

#     # makes segements from those indices
#     t_gamma = t_peaks[indices]
#     y_gamma = y_peaks[indices]

#     # move t to zero
#     t_gamma = t_gamma - t_gamma[0]

#     # normalize y data to be between 1 and 0
#     y_gamma = y_gamma / y_gamma[0]

#     # analytical solution from least squares:
#     # gamma = -sum(x * ln(y)) / sum(x^2)
#     log_y = np.log(y_gamma)
#     gamma = -np.dot(t_gamma, log_y) / np.dot(t_gamma, t_gamma)

#     return gamma

def est_gamma(peak_index, x_peaks, y_peaks):
    # step by 2 starting from peak_index
    indices = np.arange(peak_index, len(x_peaks), 2)
    x_gamma = x_peaks[indices]
    y_gamma = y_peaks[indices]

    # start from 0, normalize
    x_gamma = x_gamma - x_gamma[0]
    y_gamma = y_gamma / y_gamma[0]

    # define fit and solve
    def exp_model(x, gamma):
        return np.exp(-gamma * x)

    popt, _ = curve_fit(exp_model, x_gamma, y_gamma, p0=0.5)
    gamma = popt[0]

    return gamma

def get_constants(gamma_est, omega_est, phi_t):
    # defining a constant that is used later
    gam = np.sqrt((omega_est*omega_est) - (gamma_est*gamma_est))
    print("gam: ", gam)

    # In shifted coordinates tau = t - t_peak:
    # phi(tau) = e^(-gamma*tau) * (c1*cos(gam*tau) + c2*sin(gam*tau))
    # At tau=0: phi(0) = c1 = phi_t
    # At tau=0: phi'(0) = -gamma*c1 + gam*c2 = 0 => c2 = (gamma/gam)*c1
    c_1 = phi_t
    c_2 = (gamma_est/gam) * phi_t

    return c_1, c_2

def fit_phi(t_seg, y_seg, gamma, w_0, c1, c2, t0):
    def ode_model(t, gamma, w0, c1, c2):
        wd = np.sqrt(w0**2 - gamma**2)
        print("wd :", wd)
        return np.exp(-gamma * (t - t0)) * (c1 * np.cos(wd * (t - t0)) + c2 * np.sin(wd * (t - t0)))

    p0 = [gamma, w_0, c1, c2]
    bounds = ([0, 0, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf])
    popt, _ = curve_fit(ode_model, t_seg, y_seg, p0=p0, bounds=bounds)

    gamma, w_0, c1, c2 = popt
    return gamma, w_0, c1, c2
