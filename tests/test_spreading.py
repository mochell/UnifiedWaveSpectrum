# %%
import os, sys
import pandas as pd
import imp

exec(open(os.environ['PYTHONSTARTUP']).read())
exec(open(STARTUP_2022_particle_waves).read())


plot_path = mconfig['paths']['plot'] + '/postprocessing/unified_spectrum_test/'
# %%
import UnifiedWaveSpectrum.spreading as USpec
import UnifiedWaveSpectrum.unified as unified


# %%
imp.reload(USpec)

Omega = 0.85
u_10 = 10 # m/s
k =np.logspace(-2,3, 800)

font_for_pres()
F = M.figure_axis_xy( 5.4, 3.5)

for u_10 in np.arange(3, 23, 2):
    k_p = USpec.k_p(Omega, u_10)
    plt.semilogx(k,  USpec.Delta(k, k_p, u_10) , 'k', lw = 0.8)

    # testing  sub-function
    # u_star = USpec.u_star(u_10)
    # plt.semilogx(k,  USpec.Delta_ustar(k, k_p, u_star) )

plt.grid()
plt.title('Unified Spreading Function Delta', loc='left')

F.save(path = plot_path, name='test_spreading_Delta')

# %%

# from
#https://github.com/pakodekker/oceansar/blob/ffc6f13be34fef8ab9e7570a005b130e6462688e/oceansar/spread/elfouhaily.py

# rho = 1000.   # Density of water in kg/m^3
# S = 0.072     # Surface tension of water in N/m
# X_0 = 22e3    # Dimensionless fetch\
# g = 9.81

# #@njit(parallel=True)
# def elfouhaily(k, theta, U_10, Omega_c):
#     # Eq. 3 (below)
#     k_0 = g/U_10**2
#     # Eq. 4 (below)
#     #X = k_0*fetch
#     # Eq. 37
#     #Omega_c = 0.84*np.tanh((X/X_0)**(0.4))**(-0.75)
    
#     cK = np.sqrt(g/k + S/rho*k)

#     # Eq. 3 (below)
#     k_p = k_0*Omega_c**2
#     cK_p = np.sqrt(g/k_p + S/rho*k_p)

#     # Eq. 24
#     k_m = np.sqrt(rho*g/S)
#     cK_m = np.sqrt(g/k_m + S/rho*k_m)

#     # (McDaniel, 2001, above Equation 3.9)
#     C_10 = (0.8 + 0.065*U_10) * 1e-3
#     # Eq. 61
#     ustar = np.sqrt(C_10)*U_10

#     # Eq. 59
#     a_0 = np.log(2.)/4.
#     a_p = 4.
#     a_m = 0.13*ustar/cK_m

#     # Eq. 57
#     Delta = np.tanh(a_0 + a_p*(cK/cK_p)**2.5 + a_m*(cK_m/cK)**2.5)
#     # Eq. 49
#     G = np.where((theta > -np.pi/2.) & (theta < np.pi/2.), (1. + Delta*np.cos(2.*theta))/(np.pi), 0)

#     print(ustar, cK_m, a_m)

#     return G, Delta


# G, Delta1 = elfouhaily(k, 0, u_10, Omega  )
# plt.semilogx(k,  Delta1 , 'r-')

# %%
imp.reload(USpec)
k =np.logspace(-2,3, 800)
varphi  = np.linspace(-np.pi, np.pi, 100)
kk, vv = np.meshgrid(k, varphi)


Phi = USpec.Phi(kk, vv, k_p, u_10) 
font_for_pres()
F =M.figure_axis_xy(8.7, 6.5, container = True, view_scale= 0.7)
plt.suptitle('Directional distribution function Phi')

clev = np.arange(0, 0.4, 0.05)

ax3 = plt.subplot(2,2,1)

plt.contourf(vv, kk,  Phi, levels = clev )

ax3.set_yscale('log')
ax3.set_ylim(5, 1e3)
ax3.set_xlabel('directional angle phi [rad]')
ax3.set_ylabel('wavenumber k [rad/m]')
ax3.set_title('Wind Sea', loc='left')

ax4 = plt.subplot(2,2,2, projection='polar')
plt.contourf(vv, kk,  Phi, levels = clev )
plt.colorbar()
ax4.set_yscale('log')
ax4.set_ylim(5,  1e3)


ax1 =  plt.subplot(2,2,3)
#USpec.Phi(ki, vi, k_p, u_10)
plt.contourf(vv, kk,  Phi, levels = clev )
#plt.colorbar()
ax1.set_yscale('log')
ax1.set_xlabel('directional angle phi [rad]')
ax1.set_ylabel('wavenumber k [rad/m]')
ax1.set_ylim(8 * 1e-3, 2)
ax1.set_title('Swell', loc='left')

# make polar plot
ax2 = plt.subplot(2,2,4, projection='polar')
plt.contourf(vv, kk,  Phi, levels = clev )
plt.colorbar()
ax2.set_yscale('log')
ax2.set_ylim(8 * 1e-3, 2)

F.save(path = plot_path, name='test_spreading_Phi')

# %%

(Phi < 0).sum()

#np.trapz(Phi, k, axis = 1)
plt.plot(k, Phi.sum(axis= 0) / (2 *np.pi) )
plt.ylim(0, 4)

Phi.sum(0)

np.trapz(Phi, varphi +np.pi, axis = 0)

plt.plot(k, np.trapz(Phi, dx = np.diff(varphi)[0], axis = 0) )

plt.plot( k, np.trapz(Phi, varphi, axis = 0) )


np.diff(varphi)

# %%

varphi  = np.linspace(-np.pi, np.pi, 100)
D= 1
spread = (1 + D * np.cos(2 * varphi)) #/ (2 * np.pi)

plt.plot(varphi/np.pi, spread)

#     G = 