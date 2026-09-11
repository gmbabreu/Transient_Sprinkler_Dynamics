import numpy as np

def combine_data(full_t, full_phi, phi_an, proc_data_switch):

    # if(fit_index == (len(full_phi) - 1)):
    #     proc_data_switch = 0
    
    dt = np.mean(np.diff(full_t))

    t_end = None

    t_neg = np.arange(-2048 * dt, 0, dt)
    phi_neg = np.full(2048, full_phi[0])

    
    fourier_t = np.concatenate([t_neg, full_t[:2048]])
    fourier_phi = np.concatenate([phi_neg, full_phi[:2048]])

    return fourier_t, fourier_phi, t_end


def remove_noise(full_t, full_phi, threshold):
    
    dy = np.abs(np.diff(full_phi) / np.diff(full_t))  # derivative
    start_index = np.argmax(dy > threshold)

    start_mean = np.mean(full_phi[:start_index])
    if (start_mean < 0):
        full_phi[start_index:] += start_mean
    else:
        full_phi[start_index:]  -= start_mean
    
    for i in range(start_index):
        full_phi[i] = 0

    return full_phi

