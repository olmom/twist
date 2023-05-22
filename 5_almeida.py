from pylab import*
import pandas as pd
import seaborn as sns
from utils.almeida import AlmeidaOscillator 
from sklearn.linear_model import LinearRegression   
import copy
import os
from scipy.integrate import odeint
from utils.almeida import y0_Gp
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = '13'
plt.rcParams['legend.fontsize'] = '12'
plt.rcParams['xtick.labelsize'] = '12'
plt.rcParams['ytick.labelsize'] = '12'
current_palette = sns.color_palette()
colors = ['#1B9E77', '#D95F02', '#7570B3']
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

##########################################

# Parametric twist plots: ensembles of heterogeneous oscillators taken
# from bifurcation analyses from XPP AUTO
myoscillator = AlmeidaOscillator() 
obj = copy.copy(myoscillator)        

par_value = [obj.V_D, obj.G_p]
param_to_change = ['V_D', 'G_p']
frac_variation = 0.20 
n_oscs = 40

colors = ['forestgreen', 'salmon', '#1C7D86']
#colors = ['forestgreen', '#802E7B', '#EE42CA', '#F2B342', '#317EC2' ,'salmon', '#8A4F21', '#1C7D86']

fig5 = plt.figure(figsize=(6,6))

for p in range(len(param_to_change)):
    # load data
    data = pd.read_csv("./results/bifurcations/almeida_{}_allinfo.dat".format(
        param_to_change[p]), sep=" ", header=None)
    
    subset = data.loc[data[0]==3]
    subset2=data.loc[data[0]==2]
    
    idx_param = 3
    idx_period = 5
    
    idx_maxs = [6, 12, 13] #BMAL,PER,PERCRY
    labels = ['BMAL', 'PER', 'PER:CRY']
    
    for i in range(len(idx_maxs)):
        idx_max = idx_maxs[i]
        idx_min = idx_max + 8
    
        k = np.asarray(subset[idx_param])
        max_vble, min_vble = np.asarray(subset[idx_max]), np.asarray(subset[idx_min])
        periods = np.asarray(subset[idx_period])   
    
        # keep parameter values (k) that are within +-10% of the default value
        k_filt = k[(k >= par_value[p]*(1-frac_variation))  &  \
                   (k <= par_value[p]*(1+frac_variation))]
        max_vble = max_vble[(k >= par_value[p]*(1-frac_variation)) & \
                            (k <= par_value[p]*(1+frac_variation))]
        min_vble = min_vble[(k >= par_value[p]*(1-frac_variation)) & \
                            (k <= par_value[p]*(1+frac_variation))]
        periods = periods[(k >= par_value[p]*(1-frac_variation)) & \
                          (k <= par_value[p]*(1+frac_variation))]
    
        # dataframe with unstable FP (~ mean of oscillation)
        # but the parameter values in this df != parameter values 
        # from previous dataframe -> to calculate amp we have to do regression
        k2, unstable_FP = subset2[idx_param].values[1:], subset2[idx_max].values[1:]
        k2_filt = k2[(k2 >= par_value[p]*(1-frac_variation))  &  \
                     (k2 <= par_value[p]*(1+frac_variation))]
        unstable_FP = unstable_FP[(k2 >= par_value[p]*(1-frac_variation))  &  \
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
        max_vble, min_vble = max_vble[idxs], min_vble[idxs]
        periods = periods[idxs]
        mean = UFP_pred[idxs]
        amps = (max_vble - mean)/mean
    
        # default amplitude:
        k_filt_round = k_filt.flat[np.abs(k_filt - par_value[p]).argmin()]
        amp_default = amps[np.abs(k_filt - par_value[p]).argmin()]
        print(k_filt_round, amp_default)
        amp_ratio = amps/amp_default

        title_plot = '$V_D$' if param_to_change[p] == "V_D" else "$\gamma_P$"
    
        ax = fig5.add_subplot(2, 2, (2*p)+1)    
        ax.plot(periods,amp_ratio ,'o', c=colors[i], alpha=.5, markersize=5, label=labels[i])
        ax.set_xlabel('period (h)'); ax.set_ylabel('ratio of rel. amplitude\nto default amplitude')
        ax.set_title('twist for variations in {}'.format(title_plot))
        ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
        if p == 0:
            ax.set_xlim([22.2,28.2]); 
            ax.set_xticks([22.5, 25, 27.5])
        if p == 1:
            ax.set_xlim([23.2,26.8]);             
            ax.set_xticks([23.5, 25, 26.5])
        ax.set_ylim([0.65, 1.45])
        ax.legend(framealpha=0, loc='upper left', handletextpad=.1)

###########################
###########################

# Representative oscillations calculated through numerical integration
# For V_D
ax1 = fig5.add_subplot(422)
ax2 = fig5.add_subplot(424)

obj1 = copy.copy(myoscillator)   
obj1.V_D = np.array([obj1.V_D*(1-frac_variation), obj1.V_D*(1+frac_variation)])
total_days, dt = 100, 0.01
t = np.arange(0, total_days*24, dt)

sol1 = odeint(obj1.dynamics, y0_Gp[0:16], t)
sol1_LC = sol1[-int(40*24/dt):, :]
t_LC = t[-int(40*24/dt):]

x = sol1_LC[:,0::8] #BMAL
y = sol1_LC[:,6::8] #PER
z = sol1_LC[:,7::8] #PERCRY

x_norm = x/x.mean(axis=0) - 1
y_norm = y/y.mean(axis=0) - 1
z_norm = z/z.mean(axis=0) - 1

ax1.plot(t[0:int(120/dt)]/24, x_norm[-int(120/dt):,0] + 1, alpha=.75, 
    c=colors[0], linestyle='--')
ax1.plot(t[0:int(120/dt)]/24, x_norm[-int(120/dt):,1] + 1, alpha=.75, 
    c=colors[0])
ax2.plot(t[0:int(120/dt)]/24, y_norm[-int(120/dt):,0] + 1, alpha=.75, 
    c=colors[1])
ax2.plot(t[0:int(120/dt)]/24, y_norm[-int(120/dt):,1] + 1, alpha=.75, 
    c=colors[1], linestyle='--')
ax1.set_ylabel('BMAL1 conc.')
ax2.set_yticks([0, 1, 2]); ax1.set_yticks([0, 2, 4])
ax2.set_xlabel('time (days)'); ax2.set_ylabel('PER conc.')
ax1.set_aspect(0.5/ax1.get_data_ratio(), adjustable='box')
ax2.set_aspect(0.5/ax2.get_data_ratio(), adjustable='box')

# For G_p
obj2 = copy.copy(myoscillator)   
obj2.G_p = np.array([obj2.G_p*(1-frac_variation), obj2.G_p*(1+frac_variation)])

ax3 = fig5.add_subplot(426)
ax4 = fig5.add_subplot(428)

sol2 = odeint(obj2.dynamics, y0_Gp[0:16], t)
sol2_LC = sol2[-int(40*24/dt):, :]
t_LC = t[-int(40*24/dt):]

x = sol2_LC[:,0::8] #BMAL
y = sol2_LC[:,6::8] #PER
z = sol2_LC[:,7::8] #PERCRY

x_norm = x/x.mean(axis=0) - 1
y_norm = y/y.mean(axis=0) - 1
z_norm = z/z.mean(axis=0) - 1

ax3.plot(t[0:int(120/dt)]/24, x_norm[-int(120/dt):,0] + 1, alpha=.75, 
    c=colors[0])
ax3.plot(t[0:int(120/dt)]/24, x_norm[-int(120/dt):,1] + 1, alpha=.75, 
    c=colors[0], linestyle='--')
ax4.plot(t[0:int(120/dt)]/24, y_norm[-int(120/dt):,0] + 1, alpha=.75, 
    c=colors[1])
ax4.plot(t[0:int(120/dt)]/24, y_norm[-int(120/dt):,1] + 1, alpha=.75, 
    c=colors[1], linestyle='--')
ax3.set_ylabel('BMAL1 conc.')
ax4.set_yticks([0, 1, 2]); ax3.set_yticks([0, 2, 4])
ax4.set_xlabel('time (days)'); ax4.set_ylabel('PER conc.')
ax3.set_aspect(0.5/ax3.get_data_ratio(), adjustable='box')
ax4.set_aspect(0.5/ax4.get_data_ratio(), adjustable='box')

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