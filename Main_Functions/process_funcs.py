import numpy as np

def combine_data(full_t, full_phi, phi_an):
    
    dt = np.mean(np.diff(full_t))


    N = len(full_phi)

    if N < 2048:
        N_pad = 2048 - N
    else:
        print("this code is NOT built to handle data of length  > 2048!!!")

    N_added = 2048 + N_pad
    t_neg = np.arange(-(N_added) * dt, 0, dt)
    phi_neg = np.full(N_added, full_phi[0])

    
    fourier_t = np.concatenate([t_neg, full_t[:2048]])
    fourier_phi = np.concatenate([phi_neg, full_phi[:2048]])

    return fourier_t, fourier_phi, N_added


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

