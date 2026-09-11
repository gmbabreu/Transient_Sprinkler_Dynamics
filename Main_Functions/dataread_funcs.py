import os
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


def read_data(fname):
    data = np.loadtxt(fname, delimiter=',',skiprows=1)

    #pull out data
    full_t = data[:,0]
    full_phi = data[:,1]*np.pi/180 # converting from degrees to radians
    N = len(full_t)

    #find peaks using function form scipy.signal
    loc, _ = find_peaks(np.abs(full_phi))
    t_peaks = full_t[loc]
    phi_peaks = full_phi[loc]

    return full_t, full_phi, N, t_peaks, phi_peaks

def plot_data(t_data, phi_data, flag):
    # defines figure
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)

    # plots data
    line, = ax.plot(t_data, phi_data, 'o-',color='pink')
    vline = ax.axvline(x=t_data[0], color='purple', linestyle='--')

    # add slider
    ax_slider = plt.axes([0.2, 0.08, 0.6, 0.03])
    slider = Slider(ax_slider, 'Time of Choice', t_data.min(), t_data.max(), valinit=t_data[0], color="purple")

    def update(val):
        vline.set_xdata([slider.val, slider.val])
        fig.canvas.draw_idle()

    slider.on_changed(update)
    if(flag == 1):
        plt.title("Please use slider to pick time to start data fit. Exit plot once done.")
    else:
        plt.title("Please use slider to pick time to end data fit. Exit plot once done.")
    plt.show()

    return slider.val  # returns the final value when window is closed


def fit_segments(full_t, full_phi, t_peaks, t_fit):
    fit_peaks_index = np.argmin(np.abs(t_peaks - t_fit))
    fit_index = np.where(full_t == t_peaks[fit_peaks_index])[0][0]
    
    t_seg = full_t[fit_index:]
    phi_seg = full_phi[fit_index:]
    
    return fit_index, fit_peaks_index, t_fit, t_seg, phi_seg    
