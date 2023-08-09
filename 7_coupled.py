#%% Poincare model with different twist values: study of response to coupling
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import os
from utils.poincare import PoincareOscillator
from utils.rhythmic_parameters import RhythmicParameters
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = '13'
plt.rcParams['legend.fontsize'] = '12'
plt.rcParams['xtick.labelsize'] = '12'
plt.rcParams['ytick.labelsize'] = '12'
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

######################################################
######################################################

#%% 1. How do period and amplitude of a single cell oscillator
# change in response to coupling for different twist values?
twist_values = np.arange(-0.5, 0.2, .1)
K_values     = np.arange(0, 0.11, .01)
combinations = np.array([(tw, k) for tw in twist_values for k in K_values])
period, amplitude_rr, amplitude_rr2 = 24., .05, 1.0

# couple two identical oscillators (same period, amplitude, lambda, epsilon)
n_oscs = 2
y0 = np.tile([1,0], n_oscs)   
total_days, dt = 80,.01
t = np.arange(0, total_days*24, dt)

periods, periods2, amplitudes, amplitudes2 = [], [], [], []
for c in combinations:
    tw, k = c 
    obj = PoincareOscillator(amp=1, tau=period, lam=amplitude_rr, 
                             eps=tw, K=k)
    obj2 = PoincareOscillator(amp=1, tau=period, lam=amplitude_rr2, 
                             eps=tw, K=k)
    sol = odeint(obj.dynamics_cartesian, y0, t)
    sol2 = odeint(obj2.dynamics_cartesian, y0, t)
    x, x2 = sol[:, 0::2], sol2[:, 0::2]

    per = RhythmicParameters().periods(t[-int(20*24/dt):], 
                                          x[-int(20*24/dt):])
    amp = RhythmicParameters().amplitudes(t[-int(20*24/dt):], 
                                                x[-int(20*24/dt):])                                               
    print('lam={}, twist={}, K={} -> per={}, amp={}'.format(
        format(amplitude_rr, '.2f'), format(tw, '.2f'), 
        format(k, '.2f'), per, amp))
    per2 = RhythmicParameters().periods(t[-int(20*24/dt):], 
                                          x2[-int(20*24/dt):])
    amp2 = RhythmicParameters().amplitudes(t[-int(20*24/dt):], 
                                                x2[-int(20*24/dt):])
    per, per2, amp, amp2 = per[1], per2[1], amp[1], amp2[1]
    periods.append(per); amplitudes.append(amp)
    periods2.append(per2); amplitudes2.append(amp2)
    
fig9 = plt.figure(figsize=(11,4))
suppfig2 = plt.figure(figsize=(11,7))
suppfig3 = plt.figure(figsize=(11,4))
colors = plt.cm.Dark2_r(np.linspace(0,1,len(twist_values)))
alphas = [1, .25, .25, .25, .25, .25, 1]

ax1 = fig9.add_subplot(121)
axs4_4 = suppfig2.add_subplot(212)
axs5_1 = suppfig3.add_subplot(121)
axs5_2 = suppfig3.add_subplot(122)

for i in range(len(twist_values)):
    tw = round(twist_values[i], 2)
    ax1.plot(K_values, periods[(11*i):((i+1)*11)], 'o-', color=colors[i],
             label='$\epsilon={}$'.format(tw), markersize=4) 
    axs5_1.plot(K_values, periods[(11*i):((i+1)*11)], 'o-', color='dodgerblue',
             alpha=alphas[i], markersize=4)       
    axs5_1.plot(K_values, periods2[(11*i):((i+1)*11)], 'o-', color='crimson',
             alpha=alphas[i], markersize=4)              
    axs4_4.plot(K_values, amplitudes[(11*i):((i+1)*11)], 'o-', color=colors[i], 
             label='$\epsilon={}$'.format(tw), markersize=4) 
    axs5_2.plot(K_values, amplitudes[(11*i):((i+1)*11)], 'o-', color='dodgerblue', 
             alpha=alphas[i], label='$\lambda={}$'.format(amplitude_rr),
             markersize=4) 
    axs5_2.plot(K_values, amplitudes2[(11*i):((i+1)*11)], 'o-', color='crimson', 
             alpha=alphas[i], label='$\lambda={}$'.format(amplitude_rr2), 
             markersize=4) 

ax1.legend(framealpha=0,loc='center left', bbox_to_anchor=(1, 0.5))
axs5_2.legend(framealpha=0,loc='center left', bbox_to_anchor=(1, 0.5))
ax1.set_yticks([6, 12, 18, 24, 30, 36, 42])
axs5_1.set_yticks([6, 12, 18, 24, 30, 36, 42])
ax1.set_ylim([6.5,44.5]); axs5_1.set_ylim([6.5,44.5])
axs4_4.set_ylim([1.8, 4.35]); axs5_2.set_ylim([1.8, 4.3])
ax1.set_aspect(1.0/ax1.get_data_ratio(), adjustable='box')
axs5_1.set_aspect(1.0/axs5_1.get_data_ratio(), adjustable='box')
axs4_4.set_aspect(1.0/axs4_4.get_data_ratio(), adjustable='box')
axs5_2.set_aspect(1.0/axs5_2.get_data_ratio(), adjustable='box')
ax1.set_xlabel('coupling strength $K$'); 
axs5_1.set_xlabel('coupling strength $K$'); 
axs4_4.set_xlabel('coupling strength $K$')
axs5_2.set_xlabel('coupling strength $K$')
ax1.set_ylabel('period (h)'); axs4_4.set_ylabel('rel. amplitude')
axs5_1.set_ylabel('period (h)'); axs5_2.set_ylabel('rel. amplitude')

######################################################################
######################################################################
#%% 2. Now couple three cells with different twist values (otherwise identical)
# choose parameters here
twist = np.array([0.10, 0.0, -0.10]) #loss of sync for large twist (~1)
K_coup = 0.1 #turn on coupling

# time array
total_days = 10
n_days = 5 #to show in time series
t = np.arange(0, total_days*24, dt)

n_oscs = len(twist)
y0 = np.tile([1,0], n_oscs)   

axs4_1 = suppfig2.add_subplot(231)
axs4_2 = suppfig2.add_subplot(232)
axs4_3 = suppfig2.add_subplot(233)
ax2 = fig9.add_subplot(122)
ax2.set_title('turning on mean field coupling at K=${}$'.format(K_coup))
ax2.axvline(x=t[int(len(t)/5)]/24, color='black', lw=.5, linestyle='dashed')

# object for timeseries without and with zeitgeber
obj_no_coup = PoincareOscillator(amp=1, tau=period, lam=amplitude_rr, 
                                 eps=twist)
obj = PoincareOscillator(amp=1, tau=period, lam=amplitude_rr, 
                         eps=twist, K=K_coup)
obj1 = PoincareOscillator(amp=1, tau=period, lam=amplitude_rr, 
                         eps=np.tile(twist[0], n_oscs), K=K_coup)
obj2 = PoincareOscillator(amp=1, tau=period, lam=amplitude_rr, 
                         eps=np.tile(twist[1], n_oscs), K=K_coup)                 
obj3 = PoincareOscillator(amp=1, tau=period, lam=amplitude_rr, 
                         eps=np.tile(twist[2], n_oscs), K=K_coup)

# solve ODEs
sol_nocoup = odeint(obj_no_coup.dynamics_cartesian, y0, t[0:int(len(t)/5)])
sol_coup = odeint(obj.dynamics_cartesian, sol_nocoup[-1,:], t[int(len(t)/5):])
sol_coup_1 = odeint(obj1.dynamics_cartesian, sol_nocoup[-1,:], t[int(len(t)/5):])
sol_coup_2 = odeint(obj2.dynamics_cartesian, sol_nocoup[-1,:], t[int(len(t)/5):])
sol_coup_3 = odeint(obj3.dynamics_cartesian, sol_nocoup[-1,:], t[int(len(t)/5):])
mean_field = [obj.mean_field(sol_coup[i,:]) \
              for i in range(int(len(sol_coup)))]
mean_field1 = [obj.mean_field(sol_coup_1[i,:]) \
              for i in range(int(len(sol_coup_1)))]
mean_field2 = [obj.mean_field(sol_coup_2[i,:]) \
              for i in range(int(len(sol_coup_2)))]
mean_field3 = [obj.mean_field(sol_coup_3[i,:]) \
              for i in range(int(len(sol_coup_3)))]
mean_field = np.asarray(mean_field)[:,0]
mean_field1 = np.asarray(mean_field1)[:,0]
mean_field2 = np.asarray(mean_field2)[:,0]
mean_field3 = np.asarray(mean_field3)[:,0]
MF = np.concatenate((np.repeat(0, int(len(t)/5)), mean_field))
MF1 = np.concatenate((np.repeat(0, int(len(t)/5)), mean_field1))
MF2 = np.concatenate((np.repeat(0, int(len(t)/5)), mean_field2))
MF3 = np.concatenate((np.repeat(0, int(len(t)/5)), mean_field3))

sol = np.concatenate((sol_nocoup, sol_coup))
sol1 = np.concatenate((sol_nocoup, sol_coup_1))
sol2 = np.concatenate((sol_nocoup, sol_coup_2))
sol3 = np.concatenate((sol_nocoup, sol_coup_3))

for o in range(len(twist)):
    x_osc = sol[:, 0::2][:,o]
    x_osc_1 = sol1[:, 0::2][:,o]
    x_osc_2 = sol2[:, 0::2][:,o]
    x_osc_3 = sol3[:, 0::2][:,o]
    axs4_1.plot(t/24, x_osc_2, c=colors[-2],
                label='twist={}'.format(twist[o]), alpha=.5) 
    axs4_2.plot(t/24, x_osc_1, c=colors[-1], 
                label='twist={}'.format(twist[o]), alpha=.5) 
    axs4_3.plot(t/24, x_osc_3, c=colors[-3], 
                label='twist={}'.format(twist[o]), alpha=.5)         
    ax2.plot(t/24, x_osc, c=colors[-(o+1)],
                label='$\epsilon={}$'.format(twist[o]), alpha=.5)            
axs4_1.plot(t[int(len(t)/5):]/24, MF2[int(len(t)/5):], c='k', 
            linestyle='dashed', label='mean field')            
axs4_2.plot(t[int(len(t)/5):]/24, MF1[int(len(t)/5):], c='k', 
            linestyle='dashed', label='mean field')            
axs4_3.plot(t[int(len(t)/5):]/24, MF3[int(len(t)/5):], c='k', 
            linestyle='dashed', label='mean field')            
ax2.plot(t[int(len(t)/5):]/24, MF[int(len(t)/5):], c='k', 
         linestyle='dashed', label='mean field')            
ax2.legend(framealpha=0, loc='center left', bbox_to_anchor=(1, 0.5))                         
ax2.set_xlabel('time (days)');axs4_3.set_xlabel('time (days)')
axs4_1.set_xlabel('time (days)');axs4_2.set_xlabel('time (days)')
ax2.set_ylabel('$x$'); axs4_3.set_ylabel('$x$')
axs4_2.set_ylabel('$x$'); axs4_1.set_ylabel('$x$')
ax2.set_aspect(0.33/ax2.get_data_ratio(), adjustable='box')
axs4_1.set_aspect(0.5/axs4_1.get_data_ratio(), adjustable='box')
axs4_2.set_aspect(0.5/axs4_2.get_data_ratio(), adjustable='box')
axs4_3.set_aspect(0.5/axs4_3.get_data_ratio(), adjustable='box')
axs4_1.set_ylim([-2.2,2.2]); axs4_2.set_ylim([-2.2,2.2])
axs4_3.set_ylim([-2.2,2.2]); ax2.set_ylim([-2.2,2.2])

fig9.subplots_adjust(
    top=0.945,
    bottom=0.150,
    left=0.055,
    right=0.90,
    hspace=0.465,
    wspace=0.255
)

suppfig2.subplots_adjust(
    top=0.88,
    bottom=0.11,
    left=0.11,
    right=0.9,
    hspace=0.2,
    wspace=0.35
)
suppfig3.subplots_adjust(
    top=0.8,
    bottom=0.345,
    left=0.11,
    right=0.9,
    hspace=0.2,
    wspace=0.35
)

######################################################################
######################################################################
#%% 3. Twist induced chaos upon coupling for high values of epsilon
# choose parameters here
twist = np.array([-0.1, 0.0, 2.0]) #loss of sync for large twist (~1)

# time array
total_days = 15
n_days = 5 #to show in time series
t = np.arange(0, total_days*24, dt)

n_oscs = len(twist)
y0 = np.tile([1,0], n_oscs)   

suppfig4 = plt.figure(figsize=(11,7))
suppfig4.suptitle('period=${}$h, $\lambda={}$h-1'.format(period,amplitude_rr) + \
                  ', turning on K at Kc={}'.format(K_coup))
ax = suppfig4.add_subplot(211)
ax.axvline(x=t[int(len(t)/5)]/24, color='black', lw=.5, linestyle='dashed')
colors2 = ['#736EB3', '#D95E00', 'limegreen'] 
lws = [1.5, 1.5, 2.0]

# object for timeseries without and with coupling
obj_no_coup = PoincareOscillator(amp=1, tau=period, lam=amplitude_rr, 
                                 eps=twist)
obj = PoincareOscillator(amp=1, tau=period, lam=amplitude_rr, 
                         eps=twist, K=K_coup)

# solve ODEs
sol_nocoup = odeint(obj_no_coup.dynamics_cartesian, y0, t[0:int(len(t)/5)])
sol_coup = odeint(obj.dynamics_cartesian, sol_nocoup[-1,:], t[int(len(t)/5):])

mean_field = [obj.mean_field(sol_coup[i,:]) \
              for i in range(int(len(sol_coup)))]
mean_field = np.asarray(mean_field)[:,0]
MF = np.concatenate((np.repeat(0, int(len(t)/5)), mean_field))
sol = np.concatenate((sol_nocoup, sol_coup))

for o in range(len(twist)):
    x_osc = sol[:, 0::2][:,o]
    y_osc = sol[:, 1::2][:,o]
    ax.plot(t[0:int(len(t)/5)]/24, x_osc[0:int(len(t)/5)], c=colors2[o], 
            alpha=.4, lw=lws[o])            
    ax.plot(t[int(len(t)/5):]/24, x_osc[int(len(t)/5):], c=colors2[o], 
            label='twist={}'.format(twist[o]), alpha=.9, lw=lws[o] )            
    
    ax2 = suppfig4.add_subplot(2,3,o+4)
    ax2.plot(x_osc, y_osc, c=colors2[o])
    ax2.plot(x_osc[0:int(len(t)/5)], y_osc[0:int(len(t)/5)], c='k', lw=.75)
    ax2.set_xlim([-1.8, 1.8]); ax2.set_ylim([-1.8, 1.8])
    ax2.set_xlabel('$x$'); ax2.set_ylabel('$y$')
    ax2.set_aspect(1/ax2.get_data_ratio(), adjustable='box')

ax.plot(t/24, MF, c='k', 
        linestyle='dashed', label='mean field')            
ax.legend(framealpha=0, loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_xlabel('time (days)')
ax.set_ylabel('$x$ concentration')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
suppfig4.subplots_adjust(
    top=0.88,
    bottom=0.135,
    left=0.075,
    right=0.86,
    hspace=0.44,
    wspace=0.245
)



# save figures
isExist = os.path.exists('./figures/')
if not isExist:  
    os.makedirs('./figures/')
#fig9.savefig('./figures/fig9.pdf', format='pdf')
#suppfig2.savefig('./figures/suppfig2.pdf', format='pdf') 
#suppfig3.savefig('./figures/suppfig3.pdf', format='pdf') 
#suppfig4.savefig('./figures/suppfig4.pdf', format='pdf')