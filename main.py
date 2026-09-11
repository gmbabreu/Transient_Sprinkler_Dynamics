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

# experimental spring constant kappa
kappa = 10 # units of dyn * cm (g * cm^2/s^2)

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
spin_dir = "r"

# define reynolds number of data to use
# reads from data file name, i.e. for "forward_500_trail1" put "500" here 
re = "1000"

# define trail number  of data to use
# reads from data file name, i.e. for "forward_500_trail1" put "1" here 
trial = 1

# define where data is stored on local machine
data_dir = "C:\\Users\\gabreu\\Desktop\\Transient_Sprinkler_Dynamics\\Data_Generation\\Data"
lowpass_switch = 1

# Spring constant in dyn cm/rad
K = 32102.99

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
from estimate_funcs import est_a, est_b, get_constants, fit_phi

# PROCESS FUNCS
from process_funcs import combine_data, remove_noise

# PLOT FUNCS
from plot_funcs import plot_analytical, plot_franken, plot_phi_gen, plot_torque, plot_torque_int

# FOURIER TRANSFORM FUNCS
from fft_funcs import torque_solver, phi_from_torque


# SECTION ONE - READ ANGULAR POSITION DATA AND REQUEST FIT LOCATION
###################################################################################################

os.chdir(data_dir)
data_name = re + spin_dir + str(trial) + "_data"
fname = data_name + ".csv"

print("------------------------------")
print("Data coming from ", fname, " in directory ", data_dir)
if(fit_switch == 1):
    print("Data is being fit!")
if(proc_data_switch == 1):
    print("Data is being processed!")

print("Data incoming is in degrees and being converted to radians!!")

if spin_dir == "f":
    spin_switch = 1
elif spin_dir == "r":
    spin_switch = 0
else:
    raise ValueError(f"Invalid direction '{parts[0]}' in file name - please use 'forward' or 'rev'")


# read data for time t, angle phi  data, get size of data, and find peaks of data
full_t, full_phi, N, t_peaks, phi_peaks = read_data(fname) # converts from degrees to radians


# lets user look at plot and define fitting target point t_fit (t_f)
t_fit = 47.05 #plot_data(full_t, full_phi, 1)

# returns segment of t and phi after user defined time for fitting target
fit_index,fit_peaks_index, t_fit, t_seg, phi_seg = fit_segments(full_t, full_phi, t_peaks, t_fit)

# SECTION TWO - ESTIMATE OMEGA AND GAMMA
###################################################################################################

# estimates a from data
a_est = est_a(fit_peaks_index, t_fit, t_peaks)

# estimates b from data
b_est = est_b(fit_peaks_index, t_peaks, phi_peaks)


# estimate c from a and b
c_est = np.sqrt((a_est*a_est) + (b_est*b_est))

# gets constants for ODE using ours estimates for gamma and omega
# given fit t value
C1_est, C2_est = get_constants(b_est, c_est, full_phi[fit_index])

print("---------------------------------------")
print("-----ESTIMATES FROM DATA (NO FIT)_-----")
print("---------------------------------------")
print("Estimate of b: ", b_est)
print("Estimate of c: ", c_est)

# if user wants a fit...
if(fit_switch == 1):
    # fit the data
    b, c, C_1, C_2 = fit_phi(t_seg, phi_seg, b_est, c_est, C1_est, C2_est, full_t[fit_index])
    print("------------------------------")
    print("-------VALUES AFTER FIT-------")
    print("------------------------------")
    print("b: ", b)
    print("c: ", c)
    print("------------------------------")
    

else:
    # set end values as estimates if no fit
    b = b_est
    c = c_est
    C_1 = C1_est
    C_2 = C2_est
    
# define the analytical solution in terms of new found values
def phi_an(t):
    t0 = t_peaks[fit_peaks_index]
    a = np.sqrt(c**2 - b**2)
    return np.exp(-b * (t - t0)) * (C_1 * np.cos(a * (t - t0)) + C_2 * np.sin(a * (t - t0)))

# SECTION THREE - CLEANS NOISE FROM DATA
###################################################################################################

idexComb = fit_index
# threshold is how much you want to clean the data! bigger = smoother ( .01 is nice for rev , 0.5 for forward )
if(proc_data_switch == 1):
    full_phi = remove_noise(full_t, full_phi, threshold=.5)
    


fourier_t, fourier_phi, t_end = combine_data(full_t, full_phi, phi_an, proc_data_switch)

# SECTION FOUR - FOURIER TRANSFORM
###################################################################################################

N_f = 2*2048
L = fourier_t[-1] - fourier_t[0]
h = L/N_f

# outputs extracted torque signal
extracted_torque = torque_solver(N_f, b, c, L, fourier_phi)

stopper = False
for i in range(len(fourier_phi)):
    if stopper:
        break
    else:
        if (fourier_phi[i] != 0):
            turnOnIdx = i - 2048
            stopper = True
        
# print("turnOnIdx :", turnOnIdx)
# print("index: ", index)

extracted_torque = np.real(extracted_torque)
# plt.plot(fourier_t, extracted_torque, label="Before smooth")

if(lowpass_switch == 1):
    sos = cheby2(N=4, rs=60, Wn=0.5, btype="low", output="sos")
    extracted_torque = sosfiltfilt(sos, extracted_torque)


    # sos = cheby2(N=7, rs=40, Wn=0.8, btype="low", output="sos")
    # signal[turnOnIdx:(index+2048)] = sosfiltfilt(sos, signal[turnOnIdx:(index+2048)])

    # sos = cheby2(N=7, rs=40, Wn=0.1, btype="low", output="sos")
    # signal[(index + 2048):] = sosfiltfilt(sos, signal[(2048 + index):])
    
    # b,a = butter(N=7, Wn=0.5, btype="low")
    # signal[turnOnIdx:(2048 + index)] = filtfilt(b, a, signal[turnOnIdx:(2048 + index)])

    # b,a = butter(N=7, Wn=0.1, btype="low")
    # signal[(2048 + index):] = filtfilt(b, a, signal[(2048 + index):])
    
# plt.plot(fourier_t, extracted_torque, label="after smooth")
# plt.legend()
# plt.show()


# calculates torque signal integral
torque_integral = np.trapezoid(extracted_torque, fourier_t)

# going forwards
fourier_t_seg   = fourier_t[N_f//2:]
fourier_phi_seg   = fourier_phi[N_f//2:]
extracted_torque_seg = extracted_torque[N_f//2:]

# SECTION FIVE - DO FORWARD PROBLEM WITH TORQUE SIGNAL
###################################################################################################

phi_forward = phi_from_torque(N_f, fourier_t, extracted_torque, b, c)


error_forward = np.sqrt(np.trapezoid(((phi_forward - fourier_phi_seg)**2), x=fourier_t_seg)) / np.sqrt(np.trapezoid((fourier_phi_seg**2), x=fourier_t_seg))
print("Error forward: ", error_forward)

# SECTION SIX - COMPUTE TORQUE INTEGRAL
###################################################################################################



I = kappa/c**2 # solve for inertia, in units of g cm^2
delta = b*I # solve for damping coefficent, in units of g cm^2 s^-1
xi = delta/(2*np.sqrt(kappa*I)) # solve for damping ration, dimensionless
w0 = (delta/I)/(2*xi) # solve for natural frequency, units of s^-1
T0 = (2*pi)/w0 # solve for natural, undamped period, units of s
Td = 1/(xi*w0) # solve for decay timesacle, units of s


print("-------------------------------")
print("-------SYSTEM PARAMETERS-------")
print("-------------------------------")
print("Moment of Intertia, I: ", I, " g cm^2")
print(r"Damping coefficent, $\delta$: ", delta, " g cm^2 s^-1")
print("-------------------------------")

true_torque = extracted_torque_seg*I # units of dyn * cm

cummInt = integrate.cumulative_trapezoid(np.real(true_torque), fourier_t_seg, initial=0)

# SECTION SEVEN - SAVE SIGNAL AND PLOT RESULTS
###################################################################################################



if(write_t == 1):
    out_fname = data_name + "_signal.csv"
    out_path = os.path.join(data_dir, out_fname)
    if os.path.exists(out_path):
        os.remove(out_path)
        print("Existing file removed: ", out_path)
    np.savetxt(out_path, np.column_stack([fourier_t[N_f//2:], np.real(extracted_torque[N_f//2:])]),
           delimiter=',', header='t,extracted_torque', comments='')
    print("Extractes torque saved to: ", out_path)


if(plot_switch == 1):
    # plot_analytical(full_t, full_phi, phi_an(full_t), fit_index)
    # plot_franken(fourier_t, fourier_phi, full_t, fit_index, t_end)
    plot_torque(fourier_t, true_torque)
    plot_torque_int(fourier_t_seg, cummInt)
    # plot_phi_gen(fourier_t_seg, fourier_phi_seg, phi_forward)
    plt.show()
