import numpy as np

from scipy.signal import butter, filtfilt, welch

def combine_data(full_t, full_y, index, phi_an, proc_data_switch):

    if(index == (len(full_y) - 1)):
        proc_data_switch = 0
    
    dt = np.mean(np.diff(full_t))

    t_end = None

    t_neg = np.arange(-2048 * dt, 0, dt)
    y_neg = np.full(2048, full_y[0])

    
    franken_t = np.concatenate([t_neg, full_t[:2048]])
    franken_y = np.concatenate([y_neg, full_y[:2048]])



    # OLD CODE TO INSERT ANALYITCAL SOLUTION INTO TAIL - NOT VERY GOOD
    # there is something in the tail torque that brings error
    # of forward down significantly


    spot = full_t[index]  # t value at peak

    if(proc_data_switch == 1):
     # define end of positive t
        n_end = 2048 - len(full_t[:index])
        t_end = np.linspace(spot + dt, spot + dt * n_end, n_end)
    else:
        t_end = None

    if(proc_data_switch == 1):
         # define negative t
        t_neg = np.arange(-2048 * dt, 0, dt)
        y_neg = np.zeros(len(t_neg))
    else:
        t_neg = np.arange(-2048 * dt, 0, dt)
        y_neg = np.full(2048, full_y[0])



    if(proc_data_switch == 1):
         # combine
        franken_t = np.concatenate([t_neg, full_t[:index], t_end])
        franken_y = np.concatenate([y_neg, full_y[:index], phi_an(t_end)])
    else:
        franken_t = np.concatenate([t_neg, full_t[:2048]])
        franken_y = np.concatenate([y_neg, full_y[:2048]])

    return franken_t, franken_y, t_end


def remove_noise(full_t, full_y, threshold):
    
    dy = np.abs(np.diff(full_y) / np.diff(full_t))  # derivative
    start_index = np.argmax(dy > threshold)

    start_mean = np.mean(full_y[:start_index])
    if (start_mean < 0):
        full_y[start_index:] += start_mean
    else:
        full_y[start_index:]  -= start_mean
    
    for i in range(start_index):
        full_y[i] = 0

    return full_y

