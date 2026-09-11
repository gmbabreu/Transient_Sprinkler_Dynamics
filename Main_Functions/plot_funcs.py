import matplotlib.pyplot as plt
import numpy as np

def plot_analytical(full_t, full_y, y_fit, index):
    fig, ax = plt.subplots()

    # data as points
    ax.scatter(full_t, full_y, color='hotpink', s=10, label='Original Data', zorder=2)

    # analytical fit as line
    ax.plot(full_t, y_fit, color='mediumpurple', linewidth=2, label='Analytical Solution', zorder=3)

    # lime green circle around index
    ax.scatter(full_t[index], full_y[index], s=200, facecolors='none', 
               edgecolors='limegreen', linewidths=2.5, zorder=4, label='Start of Data Used to Find Gamma/Omega')

        # zoom around fit region
    ax.set_xlim(full_t[index] - 1, full_t[-1] + 1)
    y_region = full_y[index:]
    ax.set_ylim(y_region.min() - 0.1, y_region.max() + 0.1)

    ax.legend()
    ax.set_xlabel('time, seconds')
    ax.set_ylabel('Angular Displacement, unsure')
    ax.set_title('Plot of Analytical Solution from gamma/omega and Original Data')
    plt.tight_layout()
    # plt.show()


def plot_franken(franken_t, franken_y, full_t, index, t_end):
    has_tail = t_end is not None
    
    if has_tail:
        n_neg = len(franken_t) - len(full_t[:index]) - len(t_end)
    else:
        checker = True
        n_neg = 2048

        

    t_zeros = franken_t[:n_neg]
    
    if has_tail:
        t_orig  = franken_t[n_neg:n_neg + index]
        t_tail  = franken_t[n_neg + index:]
    else:
        t_orig = franken_t[2048:4096]

    y_zeros = franken_y[:n_neg]
    
    if has_tail:
        y_orig  = franken_y[n_neg:n_neg + index]
        y_tail  = franken_y[n_neg + index:]
    else:
        y_orig = franken_y[2048:4096]

    plt.figure()
    plt.plot(t_zeros, y_zeros, color='cornflowerblue', linewidth=2, label='added zeros')
    plt.plot(t_orig,  y_orig,  color='hotpink',        linewidth=2, label='original data')
    
    if has_tail:
        plt.plot(t_tail,  y_tail,  color='mediumpurple',   linewidth=2, label='analytical tail')
    plt.legend()
    plt.xlabel('time, seconds')
    plt.ylabel('Angular Displacement, unsure')
    plt.tight_layout()
    plt.title('Plot of Data Used to Find Torque Signal')
    # plt.show()

def plot_torque(time, signal):
    signal = np.real(signal)
    plt.figure()
    plt.plot(time, signal, color='hotpink')
    plt.xlabel('time, seconds')
    plt.ylabel(r'Torque, dyn $\cdot$ cm')
    plt.title('Extracted Torque Signal from Angular Data')
    
def plot_phi_gen(full_t, full_y, phi_gen):
    fig1, ax1 = plt.subplots()

    ax1.scatter(full_t, full_y, color='hotpink', s=10, label='Experimental Data', zorder=2)
    ax1.plot(full_t, phi_gen, color='mediumpurple', linewidth=2, label='Signal Produced by Extracted Torque', zorder=3)

    ax1.set_xlabel('time, seconds')
    ax1.set_ylabel('Angular Displacement, radians')
    ax1.legend()
    plt.title('Plot of Angular Data Produced by Extracted Torque Signal')
    plt.tight_layout()

def plot_torque_int(new_t, cummInt):
    plt.figure()
    plt.plot(new_t, cummInt, color='hotpink')
    plt.xlabel('time, seconds')
    plt.ylabel('Integral of Torque (dyn cm)')
    plt.title('Cummulative Integral of Torque Signal')
