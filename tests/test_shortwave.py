# %%
import os, sys
import pandas as pd

exec(open(os.environ['PYTHONSTARTUP']).read())
exec(open(STARTUP_2022_particle_waves).read())


plot_path = mconfig['paths']['plot'] + '/postprocessing/unified_spectrum_test/'

import UnifiedWaveSpectrum as USpec
import UnifiedWaveSpectrum.longwave
import UnifiedWaveSpectrum.shortwave
import imp

# %% test a_m
u_range = np.arange(0.1, 20, 0.2)

a_m = np.zeros_like(u_range)
for u_10, i in zip(u_range, range(len(u_range))) :
    a_m[i] = USpec.shortwave.a_m(u_10)

F = M.figure_axis_xy( 5.4, 3.5)
plt.plot(u_range, a_m, label='a_m', color='k')

plt.xlabel('$u_{10}$ [m/s]')
plt.ylabel('$a_m$ [-]')
plt.title('Generalized Phillips-Kitaigorodskii\nequilibrium range parameter for short waves a_m', loc='left')  
F.save(path = plot_path, name='test_HW_a_m')

# %%

k =np.logspace(-1,4, 500)

c_m , k_m = USpec.shortwave.c_m(return_wavenumber=True)
Omega = 0.85
u_10 = 15 # m/s
k_p = USpec.longwave.k_p(u_10, Omega)
F_m = USpec.shortwave.F_m(k, k_m) * USpec.longwave.L_PM(k, k_p)
#F_m =USpec.longwave.L_PM(k, k_p/2)

F = M.figure_axis_xy( 5.4, 3.5)

plt.plot(k,  F_m , label='F_m', color='k')

#add labels and title
plt.xlabel('$k$ [rad/m]')
plt.ylabel('$F_m$ [-]')
plt.title('Short-wave side effect function F_m', loc='left')

F.save(path = plot_path, name='test_HW_F_m')


# %% extend dispersion relation

# def c(k, k_m = 370.88):
#     """
#     phase speed for a given wave number (k).
    
#     Parameters:
#     k (numpy.ndarray): Wave number vector
#     k_m (float): Secondary (gravity-capillary) peak wave number in the curvature spectrum. Default is 370.88 rad/m.

#     Returns:
#     numpy.ndarray: The phase speed at the given wave number
#     """
#     g= 9.81
#     return np.sqrt( (g / k) * (1  + k/k_m))


imp.reload(USpec)

k =np.logspace(-1,3, 500)

font_for_pres()
F = M.figure_axis_xy( 4.5, 4.5, view_scale = 0.6, container = True)

plt.loglog(k, USpec.shortwave.c(k) , c='k', label='c(k) Jaehne and Riemer (1990)')
plt.loglog(k, np.sqrt(9.81/k) , label='c(k) deep water')

plt.legend()

# %%
k =np.logspace(-2,4, 800)
Omega = 0.85

imp.reload(USpec)

font_for_pres()
F = M.figure_axis_xy( 7.4, 4.5, view_scale = 0.6, container = True)

plt.suptitle('Shortwave spectrum', fontsize = 14)

plt.subplot(1,2,1)

for u_10 in np.arange(3, 21, 4):
    a_m = USpec.shortwave.a_m(u_10)
    k_p = USpec.longwave.k_p(u_10, Omega)
    S_h = USpec.shortwave.S_h(k, a_m, k_p)
    plt.loglog( k, S_h , label=f'u_10 = {u_10} m/s' )

plt.legend()

#add labels and title
plt.xlabel('$k$ [rad/m]')
plt.ylabel('$S_h$ [m$^3$/ rad]')
plt.ylim(1e-15, 1e0)

plt.title('Elevation spectrum S_h', loc='left')

plt.subplot(1,2,2)

for u_10 in np.arange(3, 21, 4):
    a_m = USpec.shortwave.a_m(u_10)
    k_p = USpec.longwave.k_p(u_10, Omega)
    S_h = USpec.shortwave.S_h(k, a_m, k_p)
    plt.loglog( k, S_h * k**3, label=f'u_10 = {u_10} m/s' )

plt.legend()

#add labels and title
plt.xlabel('$k$ [rad/m]')
plt.ylabel('$B_h$ [rad$^2$]')
plt.title('Curvature spectrum B_h', loc='left')
plt.ylim(1e-5, 1e0)

F.save(path = plot_path, name='test_HW_S_h')
# %%