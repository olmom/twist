import os
import numpy as np
from scipy.integrate import odeint
from astropy.timeseries import LombScargle
import matplotlib.pyplot as plt
import multiprocessing as mp
from functools import partial
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = '13'
plt.rcParams['legend.fontsize'] = '12'
plt.rcParams['xtick.labelsize'] = '12'
plt.rcParams['ytick.labelsize'] = '12'
#%%----------------------------------------------------------------------------
# Function for oscillatory properties
#------------------------------------------------------------------------------
def Properties(res, T):
    """
    returns period, phase and amplitude of all variables of a clock model
    phase and amplitude are determined for 10 zeitgeber cycles(= 10 'days')
    """
    resu = np.zeros((len(res)-1, 21))
    lightonsets = np.arange(2*T, res[-1][-1]-T, T)
    for j in range(len(res)-1):
        #period via maximum peak in LombScargle-periodogram
        frequency, power = LombScargle(res[-1], res[j]).autopower()
        frequency        = 1./frequency
        power            = power[frequency > 0.1]
        frequency        = frequency[frequency > 0.1]   
        period = frequency[power==max(power)]
        resu[j][0] = period
        
        phi       = np.zeros(10)
        amplitude = np.zeros(10)
        for i in range(1,11):
            a = res[j][res[-1]>lightonsets[-i]]
            b = res[-1][res[-1]>lightonsets[-i]]
            a = a[b<lightonsets[-i]+T+\
                  (res[-1][-i*int(T)]-res[-1][-1-i*int(T)])]
            b = b[b<lightonsets[-i]+T+\
                  (res[-1][-i*int(T)]-res[-1][-1-i*int(T)])]
            ma = max(a)
            #phi as zero crossing of average of rising slope
            a0    = a - np.mean(a)
            zero  = a0[:-1][np.diff(np.sign(a0)) > 0.]
            if len(zero) >0.:
                izero = np.where(a0 == zero[0])
                phi[i-1] = b[izero[0][0]] - lightonsets[-i]
            else: phi[i-1] = np.nan
            mia = min(a)
            #amplitude as difference between minimum value and 
                #maximum value within the zeitgeber cycle
            amplitude[i-1] = ma - mia
        resu[j][1:11]  = phi
        resu[j][11:21] = amplitude
        
    return resu

#%%-------------------------------------------------------------------
# Amplitude-Phase-Model
#---------------------------------------------------------------------
def apm_cartesian (X, t, T=24.5, F=0.25, lam=0.05, amp=1., tau=24., eps=0):
    """
    amplitude phase model
    lam : amplitude relaxation rate 
    amp : amplitude 
    tau : period
    eps : twist
    F   : zeitgeber amplitude
    T   : zeitgeber period
    """

    x, y = X[0], X[1]
    
    Z = F * np.cos(2.*np.pi*t/T + np.pi/2) #zeitgeber

    r = np.sqrt(x**2. + y**2.)
    dphidt = 2.*np.pi/tau + eps*(amp-r)

    # ODEs, already transformed to cartesian coordinates
    dxdt = lam*x*(amp-r) - y*dphidt + Z
    dydt = lam*y*(amp-r) + x*dphidt

    return [dxdt, dydt]

#%%-----------------------------------------------------------------
# multiprocessing
#-------------------------------------------------------------------
def Worker(T, pars):
    print(mp.current_process(), T, pars[1])
    r = Entrain(pars, T)                                                    
    return r

def Entrainment(Ts, pars):    
    with mp.Pool(mp.cpu_count()) as pool:
        q = pool.map(partial(Worker, pars=pars), Ts)
    pool.close()
    pool.join()
    
    q = np.asanyarray(q)  
    return(q)

def Entrain(pars, T):
    pars[0] = T
    
    dt   = 0.01
    t    = np.arange(0, 80.*int(T), dt)
    S0   = [0.1, 0.9]
    
    S = odeint(apm_cartesian , S0, t, args=tuple(pars))
    S = np.transpose(S)
    res      = np.zeros([3, len(t)])
    res[0:2] = S
    res[2]   = t
    
    resu = Properties(res, T)
    return resu

#%%--------------------------------------------------------------------
# Simulations for the Arnold tongue -- Save results
#----------------------------------------------------------------------
twist = -0.10 #change this parameter and compare to any other twist value
pars  = [24.5, 0.05, 0.05, 1., 24., twist]

Ts = np.arange(19., 29.1, 0.1)
Fs = np.arange(0.000, 0.1, 0.001)

namepath = 'results/arnold_eps={}'.format(format(twist, '.2f'))
if not os.path.exists(namepath): os.makedirs(namepath)

for i, v in enumerate(Fs):
    print('---'+str(v)+'---')
    name = namepath + '/F={}.npy'.format(format(v, '.3f'))
    pars[1] = v
    r = Entrainment(Ts, pars)
    np.save(name, r)    