# %%
import os, sys
import pandas as pd

exec(open(os.environ['PYTHONSTARTUP']).read())
exec(open(STARTUP_2022_particle_waves).read())


plot_path = mconfig['paths']['plot'] + '/postprocessing/unified_spectrum_test/'

# %%

import UnifiedWaveSpectrum.longwave as USpec

from matplotlib.ticker import (
    AutoLocator, AutoMinorLocator)

def wavenumber2frequency(k):
    k = np.array(k, float)
    near_zero = np.isclose(k, 0)
    k[near_zero] = np.inf
    k[~near_zero] = np.sqrt(9.81 * k[~near_zero]) / (2 * np.pi)
    return k

def frequency2wavenumber(f):
    return (2 * np.pi * f)**2 / 9.81

def wavenumber2period(k):
    k = np.array(k, float)
    near_zero = np.isclose(k, 0)
    k[near_zero] = np.inf
    k[~near_zero] = 2 * np.pi / np.sqrt(9.81 * k[~near_zero])
    return k
def period2wavenumber(T):
    if type(T) == np.ndarray:
        T[T == 0] = 1e-1
    else:
        if T == 0:
            T = 1e-1
    kk = (2 * np.pi / T )**2 / 9.81
    return kk


def add_second_xaxis(ax, unit = 'frequency'):

    if unit == 'frequency':
        functions = (wavenumber2frequency, frequency2wavenumber)
        xlabel = 'f [Hz]'
    elif unit == 'period':
        functions = (wavenumber2period, period2wavenumber)
        xlabel = 'T [s]'
    else:
        raise ValueError('unit must be either frequency or period')

    ax2 = ax.secondary_xaxis('top', functions=functions)
    #ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.set_xlabel(xlabel)

    return ax2

# %% Test peak wavenumber and Omega functions

u_10_range = np.arange(3, 23, 2)

def X_non_dim(X, u_10):
    return X * 9.81 / u_10**2

X_range = np.arange(10, 1e3) * 1e3 # din fetch in km 


for ui in u_10_range:
    X = X_non_dim(X_range, ui)
    Omega = USpec.Omega_1D(X)
    plt.plot(X, 1/Omega, label = f'u_10 = {ui} m/s')
    print(f'u10 = {ui} m/s , inverse wage age min = {min(Omega)}, wave age = {max(1/Omega)}')
plt.legend(ncol = 1)
plt.xlabel('non-dim X') 
plt.ylabel('wave age')

plt.title('inverse wave age for different wind speeds and fetches', loc='left')

M.save_anyfig(plt.gcf(), path = plot_path , name =f'inverse_waveage_test.png')

# %%
import imp 
imp.reload(USpec)

# set inverse wage age parameter to 1.2 (fully deloped sea)
Omega = 0.85

USpec.k_p(Omega, u_10_range)

plt.plot(u_10_range, USpec.k_p(Omega, u_10_range))
plt.title('predicted Peak wavenumber', loc='left')

plt.xlabel('u_10 [m/s]')
plt.ylabel('k_p [rad/m]')
plt.grid(True)


# %%

#define wavenumber vector in the range of 0.01 to 1.0 rad/m
k = np.linspace(0.002, 5, 600)

# define peak wavenumber with windspeed of 10 m/s
u_10 = 10 # m/w
k_p = USpec.k_p(Omega, u_10)

L_PM = USpec.L_PM(k, k_p)
J_p = USpec.J_p(k, k_p, Omega)
S_l, F_p = USpec.S_l(k, k_p, Omega, return_F_p = True)


F = M.figure_axis_xy(4.5, 6.5, container = True)
gs = GridSpec(2, 1)

# Create the first subplot on the grid
ax1 = F.fig.add_subplot(gs[0, 0])

# Plot everything but S_l on the first subplot
ax1.plot(k, L_PM, label = 'Pierson-Moskowitz')
ax1.plot(k, J_p, label = 'J_p JONSWAP')
ax1.plot(k, F_p, label = 'F_p', color = 'k')
ax1.legend(ncol = 1)
ax1.set_xlabel('k [rad/m]')
ax1.set_xlim(k[0], k[-1])
ax1.set_xlim(period2wavenumber(20), period2wavenumber(3))
ax11 = add_second_xaxis(ax1, unit = 'period')
ax11.set_xticks([12, 8, 6, 4])

ax1.set_ylabel('Current spectrum [rad^2]')

# Create the second subplot on the grid
ax2 = F.fig.add_subplot(gs[1, 0])

# Plot S_l on the second subplot
ax2.plot(k, S_l, label = 'Unified Spectrum S_l', color = 'k')
ax2.legend(ncol = 1)
ax2.set_xlabel('k [rad/m]')
ax2.set_xlim(k[0], k[-1])
ax2.set_xlim(period2wavenumber(20), period2wavenumber(3))
ax22 = add_second_xaxis(ax2, unit = 'period')
ax22.set_xticks([12, 8, 6, 4])


ax2.set_ylabel('Displacement spectrum [m^3/rad]')

plt.suptitle('Example of Unified Spectrum components\n for u10=' + str(u_10) + 'm/s' )

plt.grid(True)
plt.tight_layout() 

F.save_light(path = plot_path , name =f'S_l_single_test.png')

# %%
k = period2wavenumber(np.linspace(0.1, 40, 1000))
#k = np.linspace(0.002, 1e2, 10000)

F = M.figure_axis_xy(5, 4, container = False)

for ui in u_10_range:
    k_p = USpec.k_p(Omega, ui)
    S_l = USpec.S_l(k, k_p, Omega, return_F_p = False)
    plt.loglog(k, S_l, label = 'b_l', color = 'k', zorder=4)
    F.ax.axvline(k_p, color = 'b', linestyle = '-', lw= 0.5, zorder=2)

plt.title('Unified Spectrum S_l for different wind speeds', loc='left')

plt.xlabel('wavenumber k [rad/m]')
plt.ylabel('Elevation Spectrum S_l [m^2 s]')
plt.grid()
plt.ylim(1e-15, 1e3)

F.save_light(path = plot_path , name =f'S_l_test.png')


# %% Testing Omega_2D function: 1st, with angle definition.
import imp 
imp.reload(USpec)


u_10 = 10
theta = np.linspace( -np.pi, np.pi, 100)
c_p = np.arange(1, 20, 0.3)
TT, CC = np.meshgrid(theta, c_p)

Omega_picles = USpec.Omega_2D( u_10 , CC , method = 'directional', theta = TT)
#Omega_picles[Omega_picles<0] = 0


F =M.figure_axis_xy(7.7, 4.5, container = True)

plt.subplot(1,2,1)

clev  = np.arange(0 , 8, 0.25)
plt.contourf(TT, CC, Omega_picles, levels = clev)
#plt.colorbar()

plt.xlabel('Wind-Wave Angle theta [rad]')
plt.ylabel('Phase speed c_p [m/s]')
plt.title('Inverse wave age parameter Omega\namplitude-angle definition', loc='left')

plt.plot(0, u_10, 'ok', markersize = 10, markerfacecolor = 'w', markeredgewidth = 2, label = 'wind direction')
plt.legend()


# % 2nd, with vector definition 
dx = 0.2
c_p_x = np.arange(-15, 15 + dx , dx)
c_p_y = np.arange(-15, 15 + dx, dx)
CC_x, CC_y = np.meshgrid(c_p_x, c_p_y)

#TT, CC = np.meshgrid(theta, c_p)
u_10_vec = np.array([5, 5])

Omega_vector = CC_x * 0 
for cxi,i in zip(c_p_x, np.arange(len(c_p_x))):
    for cxy,j in zip(c_p_y, np.arange(len(c_p_y))):
        Omega_vector[i, -j] = USpec.Omega_2D(  u_10_vec, np.array([ cxi ,cxy ]))

plt.subplot(1,2,2)
#M.figure_axis_xy(6.5, 5.5, container = False)

plt.contourf(CC_x, CC_y, Omega_vector.T, levels = clev)
plt.colorbar()

plt.xlabel('cp_x [m/s]')
plt.ylabel('cp_y [m/s]')
plt.title('Inverse wave age parameter Omega\nvector definition', loc='left')

ulim = np.linalg.norm(u_10_vec) *1.25
plt.xlim(-ulim, ulim)
plt.ylim(-ulim, ulim)


plt.plot(u_10_vec[0],u_10_vec[1], 'ok', markersize = 10, markerfacecolor = 'w', markeredgewidth = 2)
plt.legend()

F.save_light(path = plot_path , name =f'Omega_2D_test.png')

# %%
