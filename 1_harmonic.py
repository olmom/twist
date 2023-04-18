#%% Time series of a harmonic/duffing oscillator starting at different ics
# default parameters: m=1, k=0.2, beta=-1,0,1, gamma=0
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from utils.harmonic import HarmonicOscillator
from utils.rhythmic_parameters import RhythmicParameters
import os
import matplotlib.cm as cm
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = '13'
plt.rcParams['legend.fontsize'] = '12'
plt.rcParams['xtick.labelsize'] = '12'
plt.rcParams['ytick.labelsize'] = '12'
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
from utils.colors import BlueBlackRed

############### CHOOSE PARAMETERS HERE ##############
###############========================##############
n_oscs = 10
k = 0.2 # spring coefficient
damping = 0 #friction
beta = np.array([0, -1, +1]) #spring coeff. of nonlinear elasticity
#####################################################
#####################################################

# time array and initial conditions
dt = 0.01
t = np.arange(0,300,dt)
y0_2 = np.linspace(0.01, .1, n_oscs)
y0 = [[.1, y0_2[i]] for i in range(n_oscs)]
y0 = np.hstack(y0)

fig2 = plt.figure(figsize=(11,8.5)) 
n_cols = 4
colors = cm.rainbow(np.linspace(0, 1, n_oscs))

for i in range(len(beta)):
    # build object with the desired harmonic nonlinearity
    harmonic_obj = HarmonicOscillator(
        k = k,
        m = 1, 
        gamma = damping, 
        beta = beta[i]
    )
    # solve ODEs of object
    sol = odeint(harmonic_obj.dynamics, y0, t)
    x, v = sol[:,0::2], sol[:,1::2]

    # plot panels
    ax0 = fig2.add_subplot(len(beta), n_cols, n_cols*i+1) #time series
    ax1 = fig2.add_subplot(len(beta), n_cols, n_cols*i+2) #time series
    ax2 = fig2.add_subplot(len(beta), n_cols, n_cols*i+3) #phase space
    ax3 = fig2.add_subplot(len(beta), n_cols, n_cols*i+4) #twist

    # compute amplitude and period
    periods, amplitudes = [], []
    for o in range(n_oscs):
        per = RhythmicParameters().determine_period_singleosc(t, x[:,o])
        amp = RhythmicParameters().determine_amplitude_singleosc_peaktrough(
            t, x[:,o])[1]
        periods.append(per); amplitudes.append(amp)
    periods, amplitudes = np.asarray(periods), np.asarray(amplitudes)
    ax3.axvline(x=periods[0]/periods[0]*24, 
                linestyle='dashed', color='silver', lw=0.75)
    ax3.plot(periods/periods[0]*24, amplitudes, 
             marker='o', c='k', markersize=3)

    # time series and phase spaces, normalized to a 24h period 
    # and color-coded with respect to their periods
    min_period, max_period = 21, 27
    norm = plt.Normalize(min_period, max_period)
    cm = BlueBlackRed 
    z = 24*periods/periods[0]
    for o in range(n_oscs):
        if (o==0) or (o==3) or (o==6) or (o==9):
            ax1.plot(t[0:int(100/dt)]/periods[0], x[0:int(100/dt),o], 
                     c=cm(norm(z[o])), alpha=.7, lw=1)
            ax2.plot(x[:,o], v[:,o], c=cm(norm(z[o])), alpha=.7, lw=1)
        ax3.plot(z[o], amplitudes[o], 'o', c=cm(norm(z[o])))

    # plot panels: timeseries, phase space, twist
    ax1.set_xlabel('time (days)'); ax1.set_ylabel('displacement $x$')
    ax2.set_xlabel('displacement $x$'); ax2.set_ylabel('velocity $v$')
    ax3.set_xlabel('period (h)'); ax3.set_ylabel('amplitude')

    soft_or_hard = 'soft' if beta[i]<0 else 'hard'
    title = r'Harmonic oscillator ($\beta=0$)' if beta[i]==0 else \
        r'Duffing {} oscillator, $\beta={}$'.format(soft_or_hard, beta[i])
    ax0.set_title(title); ax0.axis('off')
    ax1.set_xticks([0,2,4,6])
    ax1.set_aspect(0.75/ax1.get_data_ratio(), adjustable='box')
    ax2.set_ylim([-0.225, 0.225]); ax2.set_xlim([-0.3, 0.3])
    ax2.set_aspect(0.75/ax2.get_data_ratio(), adjustable='box')
    ax3.set_xlim([21.5, 28]); ax3.set_ylim([0.1,0.6])
    ax3.set_aspect(0.75/ax3.get_data_ratio(), adjustable='box')

sm = plt.cm.ScalarMappable(cmap=cm, norm=norm)
sm.set_array([])
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
axins = inset_axes(ax2,
                   width="200%",  
                   height="7.5%",
                   loc='lower left',
                   borderpad=-5.5
                   )
fig2.colorbar(sm, ticks=np.linspace(min_period, max_period, 
                                    max_period-min_period+1), 
              boundaries=np.arange(min_period, max_period, 0.05),
              cax=axins, orientation='horizontal', label='period (h)')   

fig2.subplots_adjust(
    top=0.955,
    bottom=0.135,
    left=0.05,
    right=0.98,
    hspace=0.20,
    wspace=0.49)    

isExist = os.path.exists('./figures/')
if not isExist:  
    os.makedirs('./figures/')
#fig2.savefig('./figures/fig2.pdf', format='pdf')
