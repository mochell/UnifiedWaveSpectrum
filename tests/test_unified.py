# %%
import os, sys

exec(open(os.environ['PYTHONSTARTUP']).read())
exec(open(STARTUP_2022_particle_waves).read())


plot_path = mconfig['paths']['plot'] + '/postprocessing/unified_spectrum_test/'

import imp
import UnifiedWaveSpectrum as USpec
import UnifiedWaveSpectrum.longwave
import UnifiedWaveSpectrum.shortwave
import UnifiedWaveSpectrum.unified
import UnifiedWaveSpectrum.spreading

# %%

#k = period2wavenumber(np.linspace(0.001, 20, 1000))
#k = np.linspace(0.001, 1e4, 1000)
k =np.logspace(-3,4, 800)

Omega = 0.85

F = M.figure_axis_xy(7.5, 4, container = True, view_scale = 0.7)

u_10_range = np.arange(3, 23, 2)
#u_10_range = np.arange(3, 23, 5)

ax1 = plt.subplot(121)
ax2 = plt.subplot(122)

imp.reload(USpec.unified)

for ui in u_10_range:
    
    k_p = USpec.longwave.k_p(Omega, ui)
    S_l = USpec.longwave.S_l(k, k_p, Omega, return_F_p = False)

    a_m = USpec.shortwave.a_m(ui)
    S_h = USpec.shortwave.S_h(k, a_m, k_p, return_F_m = False)

    Psi = USpec.unified.Psi(k, 0, S_h, S_l, 0).squeeze()

    ax1.loglog(k, S_l , label = 'S_l', color = 'b', lw= 0.5, zorder=4)
    ax1.loglog(k, S_h , label = 'S_h', color = 'r', lw= 0.5, zorder=4)
    ax1.loglog(k, (S_l + S_h) , label = 'Psi', color = 'k', zorder=4)
    #F.ax.axvline(k_p, color = 'b', linestyle = '-', lw= 0.5, zorder=2)

    ax2.loglog(k, S_l * k**3, label = 'S_l', color = 'b', lw= 0.5, zorder=4)
    ax2.loglog(k, S_h * k**3, label = 'S_h', color = 'r', lw= 0.5, zorder=4)
    ax2.loglog(k, (S_l + S_h) * k**3 , label = 'Psi', color = 'k', zorder=4)

ax1.set_title('Unified Spectrum S', loc='left')
ax2.set_title('Unified Curvature Spectrum B', loc='left')

ax1.set_xlabel('wavenumber k [rad/m]')
ax1.set_ylabel('Elevation Spectrum S [m^3 rad]')

ax2.set_xlabel('wavenumber k [rad/m]')
ax2.set_ylabel('Curvature Spectrum B [rad^2]')

ax1.grid()
ax2.grid()

ax1.set_ylim(1e-15, 1e3)
ax2.set_ylim(1e-4, 1)
#plt.xlim(1e-2, 1e4)

F.save(path = plot_path , name =f'S_unified_test')

# %% Test 2d Spectrum 
imp.reload(USpec.spreading)

def Psi_generator(u_10, c_p, return_components = False):

    #Psi = USpec.unified.Psi(k, 0, S_h, S_l, 0).squeeze()
    Omega = USpec.longwave.Omega_1D_simple(u_10, c_p)
    Omega = np.where(Omega >4 , 4, Omega)
    def Psi(k, theta):
        k_p = USpec.longwave.k_p(Omega, u_10)
        S_l = USpec.longwave.S_l(k, k_p, Omega, return_F_p = False)

        a_m = USpec.shortwave.a_m(u_10)
        S_h = USpec.shortwave.S_h(k, a_m, k_p, return_F_m = False)

        Phi = USpec.spreading.Phi(k, theta, k_p, u_10)

        if return_components:
            return S_l, S_h, Phi, Omega, k_p
        else:
            return (S_l + S_h) * Phi    

    return Psi




k =np.logspace(-3,4, 600)
theta = np.linspace( -np.pi/2, np.pi/2, 100)
#Psi1= Psi_generator(10, 1, return_components=True)
Psi2= Psi_generator(10, 9, return_components=False)

#S_l, S_h, Phi, Omega, k_p = 

k2 =np.logspace(-3,4, 3)
Psi1(k2, 0)



# %%

kk, tt = np.meshgrid(k, theta)
Psi_i = Psi2(kk, tt)
B_i = Psi_i * k**3



F= M.figure_axis_xy(4, 7.7, container = True, view_scale = 0.7)

# Create a gridspec instance
gs = GridSpec(8, 1)#, height_ratios=[1, 1, 1])

# 1st panel: 2D spectrum
ax1 = plt.subplot(gs[0:4])
Blevs = np.linspace(0, 1, 11)[1:-1]* B_i.max()

c1 = ax1.contourf(np.log(kk), tt, Psi_i , cmap = 'PuBuGn') 
c2 = ax1.contour(np.log(kk), tt, B_i, colors='k', levels = Blevs , linewidths = 0.8 )

plt.colorbar(c1, ax=ax1, orientation='horizontal', label ='Elevation Spectrum S [m^3 rad]', aspect=30, pad=0.15)
ax1.set_title('2D Elevation & Curvature Spectrum', loc ='left')
ax1.set_ylabel('angle theta [rad]')
ax1.set_xlabel('log wavenumber k')

# 2nd panel: 1D spectrum
ax2 = plt.subplot(gs[4:6])
#ax2.plot(np.log(k), Psi_i.T[:, ::10] , 'k-')
ax2.plot(k, Psi_i.max(axis=0) , 'k-')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylim(1e-15, 1e3)
ax2.set_title('Omnidirectional Elevation Spectrum', loc='left')
ax2.set_ylabel('Elevation S [m^3 rad]')


# 3rd panel: Curvature Spectrum
ax3 = plt.subplot(gs[6:8])
ax3.plot(k, (Psi_i * kk**3).max(axis=0) , 'k-')
ax3.set_xscale('log')
#ax3.set_yscale('log')
ax3.set_title('Omnidirectional Curvature Spectrum', loc='left')
ax3.set_ylabel('Curvature B [rad^2]')
ax3.set_xlabel('wavenumber k [rad/m]')
