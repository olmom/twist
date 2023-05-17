import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema
from scipy.integrate import odeint
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd 
from utils.goodwin import GonzeOscillator
from utils.rhythmic_parameters import RhythmicParameters 
from sklearn.linear_model import LinearRegression
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = '13'
plt.rcParams['legend.fontsize'] = '12'
plt.rcParams['xtick.labelsize'] = '12'
plt.rcParams['ytick.labelsize'] = '12'
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
np.random.seed(230305)

############### CHOOSE PARAMETERS HERE ##############
###############========================##############
n_oscs = 100
frac_variation = 0.10 #fraction of the parameter to change around
                      #its default value
amplitude_threshold = 90 #percentage of last peak-to-trough amplitude
                         #compared to the first, to remove non LC oscs
#####################################################

# time array and initial conditions
total_days, dt = 200, 0.01
t = np.arange(0, total_days*24, dt)
y0 = np.tile([1,0.5,0.2], n_oscs)

# checking effects of changing multiple pars at the same time
# 1. changing k4, k6
obj = GonzeOscillator()
obj.k4 = np.random.uniform(obj.k4*(1-frac_variation), 
                           obj.k4*(1+frac_variation), n_oscs) 
obj.k6 = np.random.uniform(obj.k6*(1-frac_variation), 
                           obj.k6*(1+frac_variation), n_oscs)                            
sol = odeint(obj.dynamics, y0, t)
t_LC = t[0:int(40/dt)]
x = sol[-int(40/dt):, 0::3]/sol[-int(40/dt):, 0::3].mean(axis=0) -1

amps,periods= [], []
for osc in range(n_oscs):
    amp = RhythmicParameters().determine_amplitude_singleosc_peaktrough(
        t_LC, x[:,osc])
    if (amp[0].shape[0] == 0) or (abs(amp[1]) < 1e-8):
        amp = (np.array([0]),0,0)
        per = 0    
    else:
        amp = amp 
        per = RhythmicParameters().determine_period_singleosc(
            t_LC, x[:,osc]) 
    max_magn = argrelextrema(sol[:,0::3][:,osc], np.greater)[0][-15:]
    min_magn = argrelextrema(sol[:,0::3][:,osc], np.less)[0][-15:]
    amps.append(amp); periods.append(per)
amps,periods = np.asarray(amps),np.asarray(periods)
first_peak_trough = [amps[o][0][0] for o in range(n_oscs)] 
last_peak_trough = [amps[o][0][-1] for o in range(n_oscs)] 
ratio = np.asarray(last_peak_trough)/np.asarray(first_peak_trough)
amps = np.asarray([amps[o][1] for o in range(n_oscs)] )
amps = amps[ratio > amplitude_threshold/100]
periods = periods[ratio > amplitude_threshold/100]

n_bins=10
df1 = pd.DataFrame({"periods": periods[amps>0.1],
                   "amplitudes": amps[amps>0.1]})
df1['change'] = np.repeat('k4, k6', df1.periods.shape[0])


###############################

# 2. changing k2, k4
obj = GonzeOscillator()
obj.k2 = np.random.uniform(obj.k2*(1-frac_variation), 
                           obj.k2*(1+frac_variation), n_oscs) 
obj.k4 = np.random.uniform(obj.k4*(1-frac_variation), 
                           obj.k4*(1+frac_variation), n_oscs)                            
sol = odeint(obj.dynamics, y0, t)
t_LC = t[0:int(40/dt)]
x = sol[-int(40/dt):, 0::3]/sol[-int(40/dt):, 0::3].mean(axis=0) -1

amps,periods= [], []
for osc in range(n_oscs):
    amp = RhythmicParameters().determine_amplitude_singleosc_peaktrough(
        t_LC, x[:,osc])
    if (amp[0].shape[0] == 0) or (abs(amp[1]) < 1e-8):
        amp = (np.array([0]),0,0)
        per = 0    
    else:
        amp = amp 
        per = RhythmicParameters().determine_period_singleosc(
            t_LC, x[:,osc]) 
    max_magn = argrelextrema(sol[:,0::3][:,osc], np.greater)[0][-15:]
    min_magn = argrelextrema(sol[:,0::3][:,osc], np.less)[0][-15:]
    amps.append(amp); periods.append(per)
amps,periods = np.asarray(amps),np.asarray(periods)
first_peak_trough = [amps[o][0][0] for o in range(n_oscs)] 
last_peak_trough = [amps[o][0][-1] for o in range(n_oscs)] 
ratio = np.asarray(last_peak_trough)/np.asarray(first_peak_trough)
amps = np.asarray([amps[o][1] for o in range(n_oscs)] )
amps = amps[ratio > amplitude_threshold/100]
periods = periods[ratio > amplitude_threshold/100]

n_bins=10
df2 = pd.DataFrame({"periods": periods[amps>0.1],
                   "amplitudes": amps[amps>0.1]})
df2['change'] = np.repeat('k2, k4', df2.periods.shape[0])

# df summarizing periods&amplitudes of co-changes of k2&k4 or k4&k6
df = df2.merge(df1, how='outer') 

########################

# seaborn's lmplot with Pearson's R
# https://stackoverflow.com/questions/25579227/seaborn-lmplot-with-equation-and-r2-text
g = sns.lmplot(x='periods', y='amplitudes', data=df, 
               col='change', col_wrap=3, 
               line_kws={'color': 'crimson',
                         'lw': 1.5,
                         'alpha': .7}, 
               scatter_kws={'color': 'black',
                            's': 10,
                         'alpha': .5},
               height=3, aspect=1,
               sharex=False, sharey=False)
def annotate(data, **kws):
    r, p = stats.spearmanr(data['periods'], data['amplitudes'])
    signif = 'n.s.' if p>.05 else ('*' if (p<.05) and (p>.01) else\
                                   '**' if (p<.01) and (p>.001) else '***')
    n = data['periods'].shape[0]
    ax = plt.gca()
    ax.text(.05, .9, '$r={:.2f}$, {}\n$n={}$'.format(r, signif, n),
            transform=ax.transAxes, size=10)               

g.map_dataframe(annotate)
titles = ['$k_2$, $k_4$', '$k_4$, $k_6$']
axes = g.axes.flatten()
axes[0].set_ylabel("rel. amplitude")
for a in range(axes.shape[0]):
    axes[a].set_xlabel("period (h)")
    axes[a].set_title(titles[a])
g.fig.subplots_adjust(
    top=0.866,
    bottom=0.227,
    left=0.087,
    right=0.984,
    hspace=0.2,
    wspace=0.329
)

################################

# 3. now change one parameter at a time: k2, k4, k6 individually
import copy
np.random.seed(230224)

# construct object to change k2, k4, k6
myoscillator = GonzeOscillator()
obj = copy.copy(myoscillator)
parameters = ['k2', 'k4', 'k6']
par_value = [obj.k2, obj.k4, obj.k6]

ax3 = g.fig.add_subplot(133)
colors = ['#1B9E77', '#D95F02', '#7570B3']

for p in range(len(parameters)):
    # load data from bifurcation analyses, 
    # since it contains period and amp info
    data = pd.read_csv("./results/bifurcations/gonze_{}_allinfo.dat".format(
        parameters[p]), sep=" ", header=None)
    
    # dataframe with max, min, period values
    subset = data.loc[data[0]==3]
    k, max_x, min_x = np.asarray(subset[3]), np.asarray(subset[6]), \
        np.asarray(subset[9])
    periods = np.asarray(subset[5])

    # filter out deg values (k) that are within +-10% of the default value
    k_filt = k[(k >= par_value[p]*(1-frac_variation))  &  \
               (k <= par_value[p]*(1+frac_variation))]
    max_x = max_x[(k >= par_value[p]*(1-frac_variation)) & \
                  (k <= par_value[p]*(1+frac_variation))]
    min_x = min_x[(k >= par_value[p]*(1-frac_variation)) & \
                  (k <= par_value[p]*(1+frac_variation))]
    periods = periods[(k >= par_value[p]*(1-frac_variation)) & \
                      (k <= par_value[p]*(1+frac_variation))]

    # dataframe with unstable FP (~ mean of oscillation)
    # but the parameter values in this df != parameter values 
    # from previous dataframe -> to calculate amp we have to do regression
    subset2=data.loc[data[0]==2]
    k2, unstable_FP = subset2[3].values[1:], subset2[6].values[1:]
    k2_filt = k2[(k2 >= par_value[p]*(1-frac_variation))  &  \
                 (k2 <= par_value[p]*(1+frac_variation))]
    unstable_FP= unstable_FP[(k2 >= par_value[p]*(1-frac_variation))  &  \
                             (k2 <= par_value[p]*(1+frac_variation))]

    # do linear regression on the curve unstable_FP = f(k2_filt)
    x, y = k2_filt.reshape((-1,1)), unstable_FP
    model = LinearRegression().fit(x, y)
    slope, intercept = model.coef_, model.intercept_
    r_sq = model.score(x, y)
    UFP_pred = model.predict(k_filt.reshape((-1,1)))
    
    # study amplitude-period correlation in an ensemble of 100 oscillators
    # pick 100 random values from the arrays
    idxs = np.random.randint(0, len(k_filt), n_oscs)
    k_filt = k_filt[idxs]
    max_x = max_x[idxs]
    min_x = min_x[idxs]
    periods = periods[idxs]
    mean = UFP_pred[idxs]
    amplitudes = (max_x - min_x)/mean

    n = amplitudes[amplitudes>.1].shape[0]
    label = '$k_2$\n$(n={})$'.format(n) if (p == 0) else\
          ('$k_4$\n$(n={})$'.format(n) if (p == 1) else \
           '$k_6$\n$(n={})$'.format(n))
    x_limits = [22, 26]
    y_limits = [0, 1.9]
    ax3.set_xlabel('period (h)') 
    ax3.plot(periods[amplitudes>.1], amplitudes[amplitudes>.1], 'o', 
             label=label, c=colors[p], markersize=4, alpha=0.6)
    ax3.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25), 
              framealpha=0, ncol=3, columnspacing=.4, handletextpad=.03) 
    ax3.set_ylim(y_limits); ax3.set_xlim(x_limits)
    ax3.set_aspect(1.0/ax3.get_data_ratio(), adjustable='box')


plt.show()
isExist = os.path.exists('./figures/')
if not isExist:  
    os.makedirs('./figures/')
#g.fig.savefig('./figures/fig3b.pdf', format='pdf')    