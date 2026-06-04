'''
    Written by Rachel Bertaud during PhD at Colroado School of Mines, 2026
    
    To analyze the transient behavior of the sprinkler, it is imperitive to know
    the natural frequency (omega) and the damping coefficent (gammma) of the system.
    By assuming that the tail of the measured angular position data of the sprinkler
    fits the ODE for a damped, hormoic oscillator, we have found a method to find
    those values, and use them to extract the torque signal of the system.
    
    This code:
        0. Read user defined inputs that control outputs
        1. Reads in angular position data in csv format and request fit time from user
        2. Estimates omega (natural frequency) and gamma (damping coefficent)
        3. Cleans noise from the data before and after forcing
        4. Uses Fourier tranforms to extract the underlying torque signal
        5. Does a forward ODE solve using the extracted torquue signal to verify the results
        6. Computes the integral of torque over time for the system
        7. Saves extracted torque signal to data directory in .csv format and plots results of code!

    All sections of this code are clearly labelled as above for
    process trasparency.
'''

# SECTION ZERO - USER DEFINED INPUTS
###################################################################################################

# switches plots on (1) or off (0)
plot_switch = 1

# saves extracted torque signal (1) or does not (0)
write_t = 0

# fit data after producing estimates of gamma and omega (1) or no fit (0)
fit_switch = 1

# processes the data before forcing and after forcing to remove noise (1)
# or uses raw data (0)
proc_data_switch = 1

# define the spin direction of data to use
# reads from data file name, i.e. for "forward_500_trail1" put "forward" here
# rev for reverse and forward for forward
spin_dir = "rev"

# define reynolds number of data to use
# reads from data file name, i.e. for "forward_500_trail1" put "500" here 
re = "1000"

# define trail number  of data to use
# reads from data file name, i.e. for "forward_500_trail1" put "1" here 
trial = 1

# define where data is stored on local machine
data_dir = "/Users/rachelbertaud/code/Sprinkler/Transient_Dynamics/Sprinkler_Data/"

lowpass_switch = 1

# DEFINE DEPENDENCIES AND UDFs
import os
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy import integrate
from scipy.signal import butter, sosfiltfilt, cheby2, filtfilt

# enters path where the function files are 
sys.path.append(os.path.join(os.path.dirname(__file__), 'Main_Functions'))

# DATA READ FUNCS
from dataread_funcs import read_data, plot_data, fit_segments

# ESTIMATE FUNCS
from estimate_funcs import est_omega_d, est_gamma, get_constants, fit_phi

# PROCESS FUNCS
from process_funcs import combine_data, remove_noise

# PLOT FUNCS
from plot_funcs import plot_analytical, plot_franken, plot_phi_gen, plot_torque, plot_torque_int

# FOURIER TRANSFORM FUNCS
from fft_funcs import torque_solver, phi_from_torque


# SECTION ONE - READ ANGULAR POSITION DATA AND REQUEST FIT LOCATION
###################################################################################################

os.chdir(data_dir)
data_name = spin_dir + "_" + re + "_trial" + str(trial)
fname = data_name + ".csv"

print("------------------------------")
print("Data coming from ", fname, " in directory ", data_dir)
if(fit_switch == 1):
    print("Data is being fit!")
if(proc_data_switch == 1):
    print("Data is being processed!")

if spin_dir == "forward":
    spin_switch = 1
elif spin_dir == "rev":
    spin_switch = 0
else:
    raise ValueError(f"Invalid direction '{parts[0]}' in file name - please use 'forward' or 'rev'")


# read data for x,y data, get size of data, and find peaks of data
full_t, full_y, N, t_peaks, y_peaks = read_data(fname)


# lets user look at plot and define fitting target
t_target = 49.07 #plot_data(full_t, full_y, 1)

# lets user look at plot and define where we will input the data
# this is only really used for the processing (which we dont do anymore)

t_insert = full_t.max() 
# t_insert = 54.08 # plot_data(full_t, full_y, 0)


# returns segment of t and y after UD time for fitting target
index, peak_index, insert_index, t_seg, y_seg = fit_segments(full_t, full_y, t_peaks, t_target, t_insert)

# SECTION TWO - ESTIMATE OMEGA AND GAMMA
###################################################################################################

# estimates period from data
omega_d = est_omega_d(peak_index, t_insert, t_peaks)

# estimates gamma from data
gamma_est = est_gamma(peak_index, t_peaks, y_peaks)


# estimate omega from period and gamma
omega_est = np.sqrt((omega_d*omega_d) + (gamma_est*gamma_est))

# gets constants for ODE using ours estimates for gamma and omega
# given fit t value
c1_est, c2_est = get_constants(gamma_est, omega_est, full_y[index])

print("---------------------------------------")
print("-----ESTIMATES FROM DATA (NO FIT)_-----")
print("---------------------------------------")
print("Estimate of gamma: ", gamma_est)
print("Estimate of omega: ", omega_est)
# print("Estimate of c1: ", c1_est)
# print("Estimate of c2: ", c2_est)

# if user wants a fit...
if(fit_switch == 1):
    # fit the data
    gamma, omega, c1, c2 = fit_phi(t_seg, y_seg, gamma_est, omega_est, c1_est, c2_est, full_t[index])
    print("------------------------------")
    print("-------VALUES AFTER FIT-------")
    print("------------------------------")
    print("Gamma: ", gamma)
    print("Omega: ", omega)
    # print("c1: ", c1)
    # print("c2: ", c2)
    print("------------------------------")
    

else:
    # set end values as estimates
    gamma = gamma_est
    omega = omega_est
    c1 = c1_est
    c2 = c2_est
    
# define the analytical solution in terms of new found values
def phi_an(t):
    t0 = t_peaks[peak_index]
    wd = np.sqrt(omega**2 - gamma**2)
    return np.exp(-gamma * (t - t0)) * (c1 * np.cos(wd * (t - t0)) + c2 * np.sin(wd * (t - t0)))

# SECTION THREE - CLEANS NOISE FROM DATA
###################################################################################################

idexComb = index
# threshold is how much you want to clean the data! bigger = smoother ( .01 is nice for rev , 0.5 for forward )
if(proc_data_switch == 1):
    full_y = remove_noise(full_t, full_y, threshold=.01)
    idexComb = insert_index
    


franken_t, franken_y, t_end = combine_data(full_t, full_y, idexComb, phi_an, proc_data_switch)

# SECTION FOUR - FOURIER TRANSFORM
###################################################################################################

N_f = 2*2048
L = franken_t[-1] - franken_t[0]
h = L/N_f

# outputs extracted torque signal
signal = torque_solver(N_f, gamma, omega, L, franken_y)

stopper = False
for i in range(len(franken_y)):
    if stopper:
        break
    else:
        if (franken_y[i] != 0):
            turnOnIdx = i - 2048
            stopper = True
        
print("turnOnIdx :", turnOnIdx)
print("index: ", index)

signal = np.real(signal)
plt.plot(franken_t, signal, label="Before smooth")

if(lowpass_switch == 1):
    sos = cheby2(N=4, rs=60, Wn=0.5, btype="low", output="sos")
    signal = sosfiltfilt(sos, signal)


    # sos = cheby2(N=7, rs=40, Wn=0.8, btype="low", output="sos")
    # signal[turnOnIdx:(index+2048)] = sosfiltfilt(sos, signal[turnOnIdx:(index+2048)])

    # sos = cheby2(N=7, rs=40, Wn=0.1, btype="low", output="sos")
    # signal[(index + 2048):] = sosfiltfilt(sos, signal[(2048 + index):])
    
    # b,a = butter(N=7, Wn=0.5, btype="low")
    # signal[turnOnIdx:(2048 + index)] = filtfilt(b, a, signal[turnOnIdx:(2048 + index)])

    # b,a = butter(N=7, Wn=0.1, btype="low")
    # signal[(2048 + index):] = filtfilt(b, a, signal[(2048 + index):])
    
plt.plot(franken_t, signal, label="after smooth")
plt.legend()
plt.show()




# calculates torque signal integral
torque_integral = np.trapezoid(signal, franken_t)

# going forwards
new_t   = franken_t[N_f//2:]
new_y   = franken_y[N_f//2:]
new_sig = signal[N_f//2:]

# SECTION FIVE - DO FORWARD PROBLEM WITH TORQUE SIGNAL
###################################################################################################

phi_gen = phi_from_torque(N_f, franken_t, signal, gamma, omega)


err = np.sqrt(np.trapezoid(((phi_gen - new_y)**2), x=new_t)) / np.sqrt(np.trapezoid((new_y**2), x=new_t))
print("Error forward: ", err)

# SECTION SIX - COMPUTE TORQUE INTEGRAL
###################################################################################################

cummInt = integrate.cumulative_trapezoid(np.real(new_sig), new_t, initial=0)

# SECTION SEVEN - SAVE SIGNAL AND PLOT RESULTS
###################################################################################################

if(write_t == 1):
    out_fname = data_name + "_signal.csv"
    out_path = os.path.join(data_dir, out_fname)
    if os.path.exists(out_path):
        os.remove(out_path)
        print("Existing file removed: ", out_path)
    np.savetxt(out_path, np.column_stack([franken_t[N_f//2:], np.real(signal[N_f//2:])]),
           delimiter=',', header='t,torque_signal', comments='')
    print("Signal saved to: ", out_path)


if(plot_switch == 1):
    plot_analytical(full_t, full_y, phi_an(full_t), index)
    plot_franken(franken_t, franken_y, full_t, index, t_end)
    plot_torque(franken_t, signal)
    plot_torque_int(new_t, cummInt)
    plot_phi_gen(new_t, new_y, phi_gen)
    plt.show()
