import numpy as np
from scipy.optimize import curve_fit

def est_a(fit_index, t_fit, t_peaks):

    # determine where the nearest peak is to where we want to estimate the data
    estimate_near_peak = np.argmin(np.abs(t_peaks - t_fit))

    # pulls out segment of data between peak used for fit
    # and peak close to estimate value chosen
    peak_seg = t_peaks[fit_index:]
    

    # difference of every other point in peak_seg
    # estimating period B
    B = peak_seg[2:] - peak_seg[:-2]

    #print("Mean: ", np.mean(period_est))
    
    # obtain estimate for a = (b^2 + c^2)
    a_est = (2*np.pi)/np.mean(B)
    
    return a_est


def est_b(peak_index, t_peaks, phi_peaks):
    # step by 2 starting from peak_index
    print(peak_index)
    print(len(t_peaks))
    indices = np.arange(peak_index, len(t_peaks), 2)
    t_peaks = t_peaks[indices]
    phi_peaks = phi_peaks[indices]

    # normalize time and position data in terms of minimum t in t_peaks
    t_peaks = t_peaks - t_peaks[0]
    phi_peaks = phi_peaks / phi_peaks[0]

    # define standard exponential fit 
    def exp_model(x, gamma):
        return np.exp(-gamma * x)

    # fits data to exponential model defined above
    popt, _ = curve_fit(exp_model, t_peaks, phi_peaks, p0=0.5)

    # grabs output for decay parameter b
    b = popt[0]

    return b

def get_constants(b_est, c_est,  phi_fit):
    # defining estimate for a
    a_est = np.sqrt((c_est*c_est) - (b_est*b_est))


    # In shifted coordinates ts = t - t_fit:
    # phi_s(ts) = e^(-b*ts) * (c1*cos(a*ts) + c2*sin(a*ts))
    # At ts=0: phi_s(0) = c1 = phi(t_fit) = 
    # At ts=0: phi_s'(0) = -b*c1 + a*c2 = 0 => c2 = (b/a)*c1
    C_1 = phi_fit
    C_2 = (b_est/a_est) * phi_fit

    # derivied from general solution to ode during time of zero forcing^

    return C_1, C_2

def fit_phi(t_seg, phi_seg, b_est, c_est, C_1, C_2, t0):
    def ode_model(t, b, c, C_1, C_2):
        a = np.sqrt(c**2 - b**2)
        # print("wd :", wd)
        return np.exp(-b * (t - t0)) * (C_1 * np.cos(a * (t - t0)) + C_2 * np.sin(a * (t - t0)))

    p0 = [b_est, c_est, C_1, C_2]
    bounds = ([0, 0, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf])
    popt, _ = curve_fit(ode_model, t_seg, phi_seg, p0=p0, bounds=bounds)

    b, c, C_1, C_2 = popt
    return b, c, C_1, C_2
