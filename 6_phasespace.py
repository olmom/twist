#%% Poincare model with twist and without twist
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema
from scipy.integrate import odeint
import os
from utils.poincare import PoincareOscillator
import matplotlib as mpl
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = '13'
plt.rcParams['legend.fontsize'] = '12'
plt.rcParams['xtick.labelsize'] = '12'
plt.rcParams['ytick.labelsize'] = '12'

############### CHOOSE PARAMETERS HERE ##############
###############========================##############
twist = np.array([0, 0.10, -0.10])
amplitude_rr = 0.05

F_zg = 0.05
pulse = .7
#################################n####################

# time array
total_days, dt = 8, 0.01
n_days = 4 #to show in time series
t = np.arange(0, total_days*24, dt)

# for Arnold tongue
Ts = np.arange(19., 29.1, 0.1)
Fs = np.arange(0., 0.10, 0.001)

# plot results
fig7a = plt.figure(figsize = (11, 7)) #interaction with pulse-like perturbations
fig7b = plt.figure(figsize = (11, 7)) #interaction with periodic inputs

for i in range(len(twist)):
    # construct Poincare Oscillator
    poincare_obj = PoincareOscillator(
        amp = 1,
        lam = amplitude_rr,
        eps = twist[i],
        tau = 24
    )

    # plot settings
    ax1 = fig7a.add_subplot(4,3,i+1)
    ax2 = fig7a.add_subplot(4,3,i+4)
    ax3 = fig7a.add_subplot(4,3,i+7)
    ax4 = fig7a.add_subplot(4,3,i+10)
    ax11 = fig7b.add_subplot(3,3,i+1)
    ax5 = fig7b.add_subplot(3,3,i+7)
    ax6 = fig7b.add_subplot(3,3,i+4)

    ax1.set_xlim([-1.7,1.7]); ax1.set_ylim([-1.7, 1.7])
    ax1.set_aspect(1.0/ax1.get_data_ratio(), adjustable='box')
    ax1.set_title('twist $\epsilon={}$'.format(format(twist[i], '.2f')))
    ax11.set_xlim([-1.7,1.7]); ax11.set_ylim([-1.7, 1.7])
    ax11.set_aspect(1.0/ax11.get_data_ratio(), adjustable='box')
    ax11.set_title('twist $\epsilon={}$'.format(format(twist[i], '.2f')))

    # define and plot isochrones
    r_isochrones = np.arange(.001, 3, .001)
    phi_isochrones = poincare_obj.isochrones(r_isochrones)
    for isoc in range(len(phi_isochrones)):
        r_isoc, phi_isoc = r_isochrones, phi_isochrones[isoc]
        x_isoc = r_isoc * np.cos(phi_isoc)
        y_isoc = r_isoc * np.sin(phi_isoc)
        ax1.plot(x_isoc, y_isoc, c='silver')
        ax11.plot(x_isoc, y_isoc, c='silver')
    ax1.plot(0,0,color='k', marker='o',markersize=5)
    ax1.set_xlabel('$x$'); ax1.set_ylabel('$y$')
    ax11.plot(0,0,color='k', marker='o',markersize=5)
    ax11.set_xlabel('$x$'); ax1.set_ylabel('$y$')

    # solve a poincare object with no perturbation (control)
    y0 = [1,0]
    sol_ctrl = odeint(poincare_obj.dynamics_cartesian, y0, t)
    x_ctrl = sol_ctrl[:,0]
    y_ctrl = sol_ctrl[:,1]

    ax2.plot(t[0:int(n_days*24/dt)]/24, x_ctrl[0:int(n_days*24/dt)], 
             color='k', linestyle='--', lw=1.50)

    # solve a poincare object with perturbation
    phase_perturbation = np.arange(0.5, 24.5, .5) #apply at all CTs, in h
    y0 = np.tile([1,0], len(phase_perturbation))

    poincare_obj.F_perturbation = pulse 
    poincare_obj.t_perturbation = phase_perturbation

    sol = odeint(poincare_obj.dynamics_cartesian, y0, t)
    x = sol[:,0::2] 
    y = sol[:,1::2] 

    ax1.plot(x[:,6], y[:,6], color='crimson', lw=1.5)
    ax1.plot(x_ctrl, y_ctrl, color='k', lw=1)
    ax11.plot(x[:,6], y[:,6], color='crimson', lw=1.5)
    ax11.plot(x_ctrl, y_ctrl, color='k', lw=1)    
    ax2.plot(t[0:int(n_days*24/dt)]/24, x[0:int(n_days*24/dt),6], 
             color='crimson', lw=1.50)
    ax2.set_xlabel('time (days)'); ax2.set_ylabel('$x$')
    ax2.set_aspect(0.5/ax2.get_data_ratio(), adjustable='box')

    # compute phase shift (perturbed - control) and PRC
    idx_max_ctrl = argrelextrema(x_ctrl, np.greater)[0]
    x_max_ctrl = x_ctrl[idx_max_ctrl]
    t_max_ctrl = t[idx_max_ctrl]   

    phase_shifts = []
    for j in range(x.shape[1]):
        idx_max_pert = argrelextrema(x[:,j], np.greater)[0] 
        idx_min_pert = argrelextrema(x[:,j], np.less)[0] 
        x_max_pert = x[idx_max_pert, j]
        x_min_pert = x[idx_min_pert, j]
        t_max_pert = t[idx_max_pert]
        phase_shift = t_max_pert[-1] - t_max_ctrl[-1]
        phase_shifts.append(phase_shift)
        # only for CT=3 plot period-amplitude correlations:
        # amplitude defined as peak-to-trough distance,
        # period defined as peak-to-peak distance
        if j == 6:
            peak_trough = x_max_pert[0:6] - x_min_pert[0:6] 
            peak_peak = np.diff(t_max_pert)    
    phase_shifts = np.asarray(phase_shifts)
    phase_shifts[phase_shifts < -12] += 24
    phase_shifts[phase_shifts > 12] -= 24

    ax3.plot(peak_peak[1:6], peak_trough[1:6], '--', color='crimson')
    ax3.plot(peak_peak[1:6], peak_trough[1:6], 'o', color='crimson')
    ax3.set_xlabel('peak-to-peak distance (h)')
    ax3.set_ylabel('peak-to-trough\ndistance (a.u.)')
    ax3.set_xlim([23.0, 25.0])
    ax3.set_aspect(0.5/ax3.get_data_ratio(), adjustable='box')

    ax4.axhline(y=0, linestyle='--', color='grey', lw=0.5)
    ax4.plot(phase_perturbation, phase_shifts, 
             'k-', markersize=2)
    ax4.set_xlabel('circadian time CT'); ax4.set_ylabel('phase shift')
    ax4.set_xticks([0,6,12,18,24]); ax4.set_yticks([-12,-6,0,6,12]) 
    ax4.set_aspect(0.5/ax4.get_data_ratio(), adjustable='box')

    # solve a poincare object with zeitgeber input 
    zeitgeber_periods = np.arange(20,29,.1)
    zeitgeber_periods = np.arange(28,19,-0.1)
    y0 = np.tile([1,0], len(zeitgeber_periods))

    poincare_obj.F_perturbation = False 
    poincare_obj.T = zeitgeber_periods
    poincare_obj.F = F_zg  

    t = np.arange(0,50*24,dt)

    sol_zg = odeint(poincare_obj.dynamics_cartesian, y0, t)
    x = sol_zg[:,0::2] 
    y = sol_zg[:,1::2] 

    # plot arnold tongue + resonance curve
    eps = twist[i]
    periods1 = np.ones((len(Fs),len(Ts)))
    phases1  = np.ones((len(Fs),len(Ts)))
    amplitudes1  = np.ones((len(Fs),len(Ts)))

    for j, v in enumerate(Fs):
        print(j,v)
        # load data
        namepath = 'results/arnold_eps={}'.format(format(eps, '.2f'))
        name = namepath + '/F={}.npy'.format(format(v, '.3f'))
        res        = np.load(name)[::-1,1].T #y is displayed

        periods1[j] = res[0]
        phases1[j]  = res[1]
        amplitudes1[j]  = res[11]
        df = pd.DataFrame(res[1:10])
        df2 = pd.DataFrame(res[11:20])
        means = df.mean(axis = 0)
        sds   = df.std(axis = 0)
        means2 = df2.mean(axis = 0)
        sds2   = df2.std(axis = 0)
        means[sds>0.1] = np.nan #remove results without 1:1 entrainment
        means[means>12] -=24
        phases1[j] = means
        amplitudes1[j] = means2
        periods1[j][np.argwhere(np.isnan(phases1[j]))] = np.nan
        amplitudes1[j][np.argwhere(np.isnan(phases1[j]))] = np.nan

    periods1 = np.flip(periods1, axis=1)
    amplitudes1 = np.flip(amplitudes1, axis=1)
    phases1 = np.flip(phases1, axis=1)

    norm = mpl.colors.Normalize(vmin=1.,vmax=4.)
    c = ax6.imshow(amplitudes1, origin='lower', aspect='auto',
            cmap='plasma', norm=norm, zorder=3)    
    idxs_xlabels = [10,50,90]
    ax6.set_xticks(idxs_xlabels)
    xtick_labs = Ts[np.asarray(idxs_xlabels)]
    ax6.set_xticklabels([str(round(float(l), 2)) for l in xtick_labs])
    ax6.set_xlabel('$T$ (h)')
    idxs_ylabels = [0, 49, 99]
    ax6.set_yticks(idxs_ylabels)
    ytick_labs = Fs[np.asarray(idxs_ylabels)]
    ax6.set_yticklabels([str(round(float(l), 2)) for l in ytick_labs])    
    ax6.set_ylabel('zeitg. strength $F$')

    divider = make_axes_locatable(ax6)
    cax = divider.append_axes('bottom', size='10%', pad=0.55)
    cbar = fig7b.colorbar(c, cax=cax, orientation='horizontal')
    cbar.set_label('rel. amplitude (a.u.)', labelpad=8.)       

    # resonance curve values
    F_idx = np.where(Fs == F_zg)[0][0]
    ax5.plot(Ts, amplitudes1[F_idx, :],'k-', markersize=2)
    ax5.set_xlabel('zeitgeber period (h)'); 
    ax5.set_ylabel('rel. amplitude')
    ax5.set_ylim([0.0,3.9])
    ax5.set_aspect(0.5/ax5.get_data_ratio(), adjustable='box')
    ax5.set_xlim([Ts[0], Ts[-1]])
    ax5.set_xticks([20, 24, 28])
    ax5.text(.05, 0.9, '$F={}$'.format(format(F_zg, '.2f')), 
        ha='left', va='top', transform=ax5.transAxes)
    ax5.set_aspect(0.5/ax5.get_data_ratio(), adjustable='box')

fig7a.subplots_adjust(top=0.95,
        bottom=0.07,
        left=0.10,
        right=0.9,
        hspace=0.185,
        wspace=0.8)
fig7b.subplots_adjust(top=0.95,
        bottom=0.07,
        left=0.10,
        right=0.9,
        hspace=0.275,
        wspace=0.8)        

#####################################################
#####################################################
#%% Rigid oscillators are more resistant to twist effects
############### CHOOSE PARAMETERS HERE ##############
###############========================##############
twist = 0.05
amplitude_rr = np.array([0.05, 0.1, 1.0])
pulse = .7
#####################################################

t = np.arange(0, total_days*24, dt)
fig8 = plt.figure(figsize=(11,7.5))

for i in range(len(amplitude_rr)):
     # construct Poincare Oscillator
    poincare_obj = PoincareOscillator(
        amp = 1,
        lam = amplitude_rr[i],
        eps = twist,
        tau = 24
    )   

    ax1 = fig8.add_subplot(3,3,i+1) #phase spaces
    ax2 = fig8.add_subplot(3,3,i+4) 
    ax3 = fig8.add_subplot(3,3,i+7) 
    ax1.set_xlim([-1.6,1.6]); ax1.set_ylim([-1.6, 1.6])
    ax1.set_aspect(1.0/ax1.get_data_ratio(), adjustable='box')

    # isochrones
    r_isochrones = np.arange(.001, 3, .001)
    phi_isochrones = poincare_obj.isochrones(r_isochrones)
    for isoc in range(len(phi_isochrones)):
        r_isoc, phi_isoc = r_isochrones, phi_isochrones[isoc]
        x_isoc = r_isoc * np.cos(phi_isoc)
        y_isoc = r_isoc * np.sin(phi_isoc)
        ax1.plot(x_isoc, y_isoc, c='silver')
    ax1.plot(0,0,color='k', marker='o',markersize=5)
    ax1.set_xlabel('$x$'); ax1.set_ylabel('$y$')

    # Poincare object unperturbed and perturbed
    sol_ctrl = odeint(poincare_obj.dynamics_cartesian, [1,0], t)
    poincare_obj.F_perturbation = pulse
    sol_pert = odeint(poincare_obj.dynamics_cartesian, [1,0], t)

    ax1.plot(sol_pert[:,0], sol_pert[:,1],c='crimson',lw=1.5)
    ax1.plot(sol_ctrl[:,0], sol_ctrl[:,1],c='k',lw=1)
    ax1.set_xlabel('$x$'); ax1.set_ylabel('$y$')
    ax1.set_title('$\lambda={}$,\n$\epsilon={} $'.format(
        format(amplitude_rr[i], '.2f'), format(twist, '.2f')))

    # PRC also affected by lambda
    y0 = np.tile([1,0], len(phase_perturbation))

    poincare_obj.F_perturbation = pulse 
    poincare_obj.t_perturbation = phase_perturbation

    sol = odeint(poincare_obj.dynamics_cartesian, y0, t)
    x = sol[:,0::2]; x_ctrl = sol_ctrl[:,0] 

    idx_max_ctrl = argrelextrema(x_ctrl, np.greater)[0]
    x_max_ctrl = x_ctrl[idx_max_ctrl]
    t_max_ctrl = t[idx_max_ctrl]   

    phase_shifts = []
    for j in range(x.shape[1]):
        idx_max_pert = argrelextrema(x[:,j], np.greater)[0] 
        x_max_pert = x[idx_max_pert, j]
        t_max_pert = t[idx_max_pert]
        phase_shift = t_max_pert[-1] - t_max_ctrl[-1]
        phase_shifts.append(phase_shift)
        
        if j == 6:
            peak_trough = x_max_pert[0:6] - x_min_pert[0:6] 
            peak_peak = np.diff(t_max_pert)  

    phase_shifts = np.asarray(phase_shifts)
    phase_shifts[phase_shifts < -12] += 24
    phase_shifts[phase_shifts > 12] -= 24

    ax2.plot(peak_peak[1:6], peak_trough[1:6], '--', color='crimson')
    ax2.plot(peak_peak[1:6], peak_trough[1:6], 'o', color='crimson')
    ax2.set_xlabel('peak-to-peak distance (h)')
    ax2.set_ylabel('peak-to-trough\ndistance (a.u.)')
    ax2.set_xlim([23.5, 24.5])
    ax2.set_aspect(0.66/ax2.get_data_ratio(), adjustable='box')

    ax3.axhline(y=0, linestyle='--', color='grey', lw=0.5)
    ax3.plot(phase_perturbation, phase_shifts, 
             'k-', markersize=2)
    ax3.set_xlabel('circadian time CT'); ax3.set_ylabel('phase shift (h)')
    ax3.set_xticks([0,6,12,18,24]); ax3.set_yticks([-12,-6,0,6,12]) 
    ax3.set_aspect(0.66/ax3.get_data_ratio(), adjustable='box')


fig8.subplots_adjust(top=0.88,
        bottom=0.11,
        left=0.11,
        right=0.9,
        hspace=0.205,
        wspace=0.35)       


#####################################################
#####################################################
#%% Driving a system with large twist -- chaotic dynamics 
# (shear-induced chaos)
# choose parameters here
twist = np.array([0.0, 2.0]) #loss of sync for large twist (~1)
period, amplitude_rr = 24., .05

# turn on zeitgeber
F_zg, T_zg = 0.25, 24.5

# time array
total_days = 10
n_days = 5 #to show in time series
t = np.arange(0, total_days*24, dt)

# turn on Zg input and see response of oscillators with diff. twist
n_oscs = len(twist)
y0 = np.tile([1,0], n_oscs)

colors = ['darkgray', 'royalblue']

suppfig1 = plt.figure(figsize=(11,3.5))
suppfig1.suptitle('period=${}$h, $\lambda={}$h-1'.format(period,amplitude_rr) + \
             ', turning on ZG: T={}h, F={}'.format(T_zg, F_zg))
ax = suppfig1.add_subplot(111)
ax.axvline(x=t[int(len(t)/5)]/24, color='black', lw=.5, linestyle='dashed')

# object for timeseries without and with zeitgeber
obj_no_zg = PoincareOscillator(amp=1, tau=period, lam=amplitude_rr, 
                               eps=twist)
obj = PoincareOscillator(amp=1, tau=period, lam=amplitude_rr, 
                         eps=twist, F=F_zg, T=T_zg)

# solve ODEs
sol_nozg = odeint(obj_no_zg.dynamics_cartesian, y0, t[0:int(len(t)/5)])
sol_zg = odeint(obj.dynamics_cartesian, sol_nozg[-1,:], t[int(len(t)/5):])
zg = np.concatenate((np.repeat(0, int(len(t)/5)),
                     obj.zeitgeber(t[int(len(t)/5):])))
sol = np.concatenate((sol_nozg, sol_zg))

for o in range(len(twist)):
    x_osc = sol[:, 0::2][:,o]
    y_osc = sol[:, 1::2][:,o]
    ax.plot(t/24, x_osc, c=colors[o], 
            label='twist={}'.format(twist[o]), alpha=.9, lw=2.0)            

ax.plot(t/24, zg, c='k', 
        linestyle='dashed', label='zeitgeber')            
ax.legend(framealpha=0, loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_xlabel('time (days)')
ax.set_ylabel('$x$ concentration')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

suppfig1.subplots_adjust(
    top=0.88,
    bottom=0.135,
    left=0.075,
    right=0.86,
    hspace=0.2,
    wspace=0.2
)

#########################################

# save figures
isExist = os.path.exists('./figures/')
if not isExist:  
    os.makedirs('./figures/')
#fig7a.savefig('./figures/fig7a.pdf', format='pdf')
#fig7b.savefig('./figures/fig7b.pdf', format='pdf')
#suppfig1.savefig('./figures/suppfig1.pdf', format='pdf') 
#fig8.savefig('./figures/fig8.pdf', format='pdf')