import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema
from scipy.integrate import odeint
import os
import matplotlib.pyplot as plt
from utils.almeida import AlmeidaOscillator
from utils.rhythmic_parameters import RhythmicParameters 
from utils.almeida import y0_Gp, y0_VE
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = '11'
plt.rcParams['legend.fontsize'] = '10'
plt.rcParams['xtick.labelsize'] = '10'
plt.rcParams['ytick.labelsize'] = '10'
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
np.random.seed(230307)

############### CHOOSE PARAMETERS HERE ##############
###############========================##############
n_oscs = 40 
frac_variation = 0.10 #fraction of the parameter to change around
                      #its default value
amplitude_threshold = 90 #percentage of last peak-to-trough amplitude
                         #compared to the first, to remove non LC oscs
param_to_change = ['V_E', 'G_p'] #parameters whose twist effects are shown
                                 #in figure5
#####################################################

# 1. Find and change activation of E box (V_E)
import inspect
par_names      = inspect.getfullargspec(AlmeidaOscillator().__init__)[0][1:]
default_values = list(AlmeidaOscillator().__init__.__defaults__)
idx_par        = np.where(np.asarray(par_names) == param_to_change[0])[0][0]
changing_par_values = \
    np.linspace(default_values[idx_par] * (1-frac_variation), 
                default_values[idx_par] * (1+frac_variation), 
                n_oscs) 
par_values = default_values.copy()
par_values[idx_par] = changing_par_values
par_values = np.asarray(par_values)

dv = np.asarray(default_values)
obj_default = AlmeidaOscillator(dv[0], dv[1], dv[2], dv[3], dv[4], dv[5],
                                dv[6], dv[7], dv[8], dv[9], dv[10], dv[11],
                                dv[12], dv[13], dv[14], dv[15], dv[16],
                                dv[17])    
obj = AlmeidaOscillator(par_values[0], par_values[1], par_values[2],
                        par_values[3], par_values[4], par_values[5],
                        par_values[6], par_values[7], par_values[8],
                        par_values[9], par_values[10], par_values[11],
                        par_values[12], par_values[13], par_values[14],
                        par_values[15], par_values[16], par_values[17])

# time array and initial conditions
total_days, dt = 100, 0.01
t = np.arange(0, total_days*24, dt)
y0 = y0_VE # initial conditions after removing transients in prior simulation

# solution of the control Almeida + computation of amplitudes & period
sol_def = odeint(obj_default.dynamics, y0[0:8], t)
sol_def = sol_def[-int(40*24/dt):,:]/sol_def[-int(40*24/dt):,:].mean(axis=0) - 1
period_default = RhythmicParameters().determine_period_singleosc(
        t[0:int(40*24/dt)], sol_def[-int(40*24/dt):, 0])
amplitudes_default = []
for v in range(sol_def.shape[1]):
    amp_def = RhythmicParameters().determine_amplitude_singleosc(
        t[0:int(40*24/dt)], sol_def[-int(40*24/dt):, v])
    amplitudes_default.append(amp_def)
amplitudes_default = np.asarray([amplitudes_default[v][1] \
                                 for v in range(sol_def.shape[1])])

# solution of Almeida object after changing V_E -- parametric twist study
sol = odeint(obj.dynamics, y0, t)
sol_LC = sol[-int(40*24/dt):, :]
t_LC = t[-int(40*24/dt):]

x = sol_LC[:,0::8] #BMAL
y = sol_LC[:,5::8] #CRY

# normalize absolute values to their means & center around 0
x_norm = x/x.mean(axis=0) - 1
y_norm = y/y.mean(axis=0) - 1

# compute periods, amplitudes, magnitudes
amps, periods = [], []
amps_y, periods_y = [], []
summ, summy = [], [] #summaries containing parameter, per, amp
for osc in range(n_oscs):
    print(osc)
    amp = RhythmicParameters().determine_amplitude_singleosc(
        t_LC, x_norm[:,osc])
    ampy = RhythmicParameters().determine_amplitude_singleosc(
        t_LC, y_norm[:,osc])
    per = RhythmicParameters().determine_period_singleosc(
        t_LC, x_norm[:,osc])           
    pery = RhythmicParameters().determine_period_singleosc(
        t_LC, y_norm[:,osc])            
    amps.append(amp); periods.append(per)
    amps_y.append(ampy); periods_y.append(pery)          
    summ.append(np.array([obj.V_E[osc],per,amp]))
    summy.append(np.array([obj.V_E[osc],pery,ampy]))

amps, periods = np.asarray(amps), np.asarray(periods)
amps_y, periods_y = np.asarray(amps_y), np.asarray(periods_y)

# remove period doubling: if sd of amp is >5% of the avg amp => PD
PD = [True if amps[o][2]>.05*amps[o][1] else False for o in range(n_oscs)]
PD = np.invert(np.asarray(PD))

amps_y = np.asarray([amps_y[o][1] for o in range(n_oscs)] )
amps = np.asarray([amps[o][1] for o in range(n_oscs)] )
amps, amps_y = amps[PD], amps_y[PD]
periods, periods_y = periods[PD], periods_y[PD]

fig5 = plt.figure(figsize=(6,6))
ax1 = fig5.add_subplot(221)
ax2 = fig5.add_subplot(422)
ax3 = fig5.add_subplot(424)

colors = ['forestgreen', 'salmon']
ax1.plot(np.sort(periods[amps>.1]), 
         np.sort(amps[amps>.1])[::-1], 'o', label='BMAL1', 
    c=colors[0], alpha=.5, markersize=5)
ax1.plot(np.sort(periods_y[amps_y>.1]), 
         np.sort(amps_y[amps_y>.1]), 'o', label='CRY', 
    c =colors[1], alpha=.5, markersize=5)
ax2.plot(t[0:int(120/dt)]/24, x_norm[-int(120/dt):,0] + 1, alpha=.75, 
    c=colors[0], linestyle='dashed')
ax2.plot(t[0:int(120/dt)]/24, x_norm[-int(120/dt):,-1] + 1, alpha=.75, 
    c=colors[0])
ax3.plot(t[0:int(120/dt)]/24, y_norm[-int(120/dt):,0] + 1, alpha=.75, 
    c=colors[1])
ax3.plot(t[0:int(120/dt)]/24, y_norm[-int(120/dt):,-1] + 1, alpha=.75, 
    c=colors[1],
    linestyle='dashed')

ax1.legend(framealpha=0, loc='center left')
ax1.set_title('twist for variations in ${}$'.format(param_to_change[0]))
ax1.set_xlim([23.8,25.8]); ax1.set_ylim([1.7,3.8])
ax1.set_xlabel('period (h)'); ax1.set_ylabel('rel. amplitude')
ax2.set_ylabel('rel. BMAL1 conc.')
ax2.set_yticks([0, 2, 4])
ax3.set_xlabel('time (days)'); ax3.set_ylabel('rel. CRY conc. ')

ax1.set_aspect(1.0/ax1.get_data_ratio(), adjustable='box')
ax2.set_aspect(0.5/ax2.get_data_ratio(), adjustable='box')
ax3.set_aspect(0.5/ax3.get_data_ratio(), adjustable='box')

#########################################
#########################################

# 2. Find and change degradation of PER (G_p)
par_names      = inspect.getfullargspec(AlmeidaOscillator().__init__)[0][1:]
default_values = list(AlmeidaOscillator().__init__.__defaults__)
idx_par        = np.where(np.asarray(par_names) == param_to_change[1])[0][0]
changing_par_values = \
    np.linspace(default_values[idx_par] * (1-frac_variation), 
                default_values[idx_par] * (1+frac_variation), 
                n_oscs) 
par_values = default_values.copy()
par_values[idx_par] = changing_par_values
par_values = np.asarray(par_values)

dv = np.asarray(default_values)
obj_default = AlmeidaOscillator(dv[0], dv[1], dv[2], dv[3], dv[4], dv[5],
                                dv[6], dv[7], dv[8], dv[9], dv[10], dv[11],
                                dv[12], dv[13], dv[14], dv[15], dv[16],
                                dv[17])    
obj = AlmeidaOscillator(par_values[0], par_values[1], par_values[2],
                        par_values[3], par_values[4], par_values[5],
                        par_values[6], par_values[7], par_values[8],
                        par_values[9], par_values[10], par_values[11],
                        par_values[12], par_values[13], par_values[14],
                        par_values[15], par_values[16], par_values[17])

# time array and initial conditions
total_days, dt = 100, 0.01
t = np.arange(0, total_days*24, dt)
y0 = y0_Gp # initial conditions after removing transients in prior simulation

# solution of the control Almeida + computation of amplitudes & period
sol_def = odeint(obj_default.dynamics, y0[0:8], t)
sol_def = sol_def[-int(40*24/dt):,:]/sol_def[-int(40*24/dt):,:].mean(axis=0) - 1
period_default = RhythmicParameters().determine_period_singleosc(
        t[0:int(40*24/dt)], sol_def[-int(40*24/dt):, 0])
amplitudes_default = []
for v in range(sol_def.shape[1]):
    amp_def = RhythmicParameters().determine_amplitude_singleosc(
        t[0:int(40*24/dt)], sol_def[-int(40*24/dt):, v])
    amplitudes_default.append(amp_def)
amplitudes_default = np.asarray([amplitudes_default[v][1] \
                                 for v in range(sol_def.shape[1])])

# solution of Almeida object after changing G_p -- parametric twist study
sol = odeint(obj.dynamics, y0, t)
sol_LC = sol[-int(40*24/dt):, :]
t_LC = t[-int(40*24/dt):]

x = sol_LC[:,0::8] #BMAL
y = sol_LC[:,5::8] #CRY

# normalize absolute values to their means & center around 0
x_norm = x/x.mean(axis=0) - 1
y_norm = y/y.mean(axis=0) - 1

# compute periods, amplitudes, magnitudes
amps, periods = [], []
amps_y, periods_y = [], []
summ, summy = [], []
for osc in range(n_oscs):
    print(osc)
    amp = RhythmicParameters().determine_amplitude_singleosc(
        t_LC, x_norm[:,osc])
    ampy = RhythmicParameters().determine_amplitude_singleosc(
        t_LC, y_norm[:,osc])
    per = RhythmicParameters().determine_period_singleosc(
        t_LC, x_norm[:,osc])           
    pery = RhythmicParameters().determine_period_singleosc(
        t_LC, y_norm[:,osc])            
    amps.append(amp); periods.append(per)
    amps_y.append(ampy); periods_y.append(pery)          
    summ.append(np.array([obj.G_p[osc],per,amp]))
    summy.append(np.array([obj.G_p[osc],pery,ampy]))

amps, periods = np.asarray(amps), np.asarray(periods)
amps_y, periods_y = np.asarray(amps_y), np.asarray(periods_y)

# remove period doubling: if sd of amp is >5% of the avg amp => PD
PD = [True if amps[o][2]>.05*amps[o][1] else False for o in range(n_oscs)]
PD = np.invert(np.asarray(PD))

amps_y = np.asarray([amps_y[o][1] for o in range(n_oscs)] )
amps = np.asarray([amps[o][1] for o in range(n_oscs)] )
amps, amps_y = amps[PD], amps_y[PD]
periods, periods_y = periods[PD], periods_y[PD]

ax4 = fig5.add_subplot(223)
ax5 = fig5.add_subplot(426)
ax6 = fig5.add_subplot(428)

ax4.plot(periods[amps>.1], 
         np.sort(amps[amps>.1])[::-1], 'o', label='BMAL1', 
    c=colors[0], alpha=.5, markersize=5)
ax4.plot(periods_y[amps_y>.1], 
         np.sort(amps_y[amps_y>.1])[::-1], 'o', label='CRY', 
    c=colors[1], alpha=.5, markersize=5)
ax5.plot(t[0:int(120/dt)]/24, x_norm[-int(120/dt):,0] + 1, alpha=.75, 
    c=colors[0])
ax5.plot(t[0:int(120/dt)]/24, x_norm[-int(120/dt):,-1] + 1, alpha=.75,
    linestyle='dashed', c=colors[0])
ax6.plot(t[0:int(120/dt)]/24, y_norm[-int(120/dt):,0] + 1, alpha=.75, 
    c=colors[1])
ax6.plot(t[0:int(120/dt)]/24, y_norm[-int(120/dt):,-1] + 1, alpha=.75, 
    linestyle='dashed', c=colors[1])

ax4.legend(framealpha=0, loc='center left')
ax4.set_xlim([23.8,25.8]); ax4.set_ylim([1.7,3.8])
ax4.set_title('twist for variations in $\gamma_p$')
ax4.set_xlabel('period (h)'); ax4.set_ylabel('rel. amplitude')
ax5.set_ylabel('rel. BMAL1 conc.')
ax5.set_yticks([0, 2, 4])
ax6.set_xlabel('time (days)'); ax6.set_ylabel('rel. CRY conc. ')

ax4.set_aspect(1.0/ax4.get_data_ratio(), adjustable='box')
ax5.set_aspect(0.5/ax5.get_data_ratio(), adjustable='box')
ax6.set_aspect(0.5/ax6.get_data_ratio(), adjustable='box')

fig5.subplots_adjust(
top=0.905,
bottom=0.11,
left=0.11,
right=0.945,
hspace=0.615,
wspace=0.325
)
isExist = os.path.exists('./figures/')
if not isExist:  
    os.makedirs('./figures/')
#fig5.savefig('./figures/fig5.pdf', format='pdf')    


