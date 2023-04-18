import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import newton 
from scipy.signal import argrelextrema

class RhythmicParameters:

    def __init__ (self, interp_length=15):
        self.interp_length = interp_length
        return 
        
    def determine_zeros(self, t, singleosc):                    
        """
        Compute zeros of a single oscillation through Newton method
        NOTE that the try/except is done because sometimes the line that
        joins 3 points together is almost parallel to x axis. In this case
        it is not able to calculate the root, and a ValueError arises (a
        value in x_new is above/below the interpolation range
        """ 
      
        singleosc = np.asarray(singleosc)
        results = []
        l = len(t)                         
        idx = np.arange(l, dtype=int)       
    
        for m in np.where(np.diff(singleosc >= 0 ) == True)[0]:
            res = []
            if idx[m] < (self.interp_length) or \
                    idx[m] > l - (self.interp_length):
                continue             
            min_i = max(idx[m] - self.interp_length, 0)   
            max_i = min(idx[m] + self.interp_length, l)
    
            time =         t[min_i:max_i]
            data = singleosc[min_i:max_i]
    
            t0 = time[0] #shift interval to avoid error by large numbers
            time = time - t0
    
            interpol = interp1d(time,data,kind='cubic',axis=0)
            try:
                root_t = newton(interpol, t[m]-t0)    
                singleosc_root = interpol(root_t).tolist()
                res.append([root_t+t0,singleosc_root])
                res = np.asarray(res[0])
    
            except ValueError:
                pass
            except RuntimeError:
                pass
            results.append(res)
        results = np.asarray(results)
        return results

    def determine_period_singleosc(self, t, singleosc):        
        """
        Determination of the period of a single oscillator:
        Two times the distance between two consecutive zeros
        """
        period = []
        times = np.hstack(self.determine_zeros(t, singleosc))[0::2]
        halfperiod = np.diff(times)
        halfperiod = halfperiod[halfperiod > 2] #filter close 0 crossings
        halfperiod_avg = halfperiod.mean()
        period = 2*halfperiod_avg
        return period

    def periods(self, t, x):
        """
        Determine periods of N oscillators 
        and return the mean and std dev
        """
        results = []
        for so in range(np.shape(x)[1]):
            period = self.determine_period_singleosc(t,x[:,so]).tolist()
            results.append(period)
        mean = np.mean(results)
        std  = np.std(results)
        return np.array(results), mean, std

    def determine_amplitude_singleosc_peaktrough(self, t, singleosc):        
        """
        Determination of the amplitude of a single oscillator:
        peak-to-trough distance
        """
        idx_max = argrelextrema(singleosc, np.greater)[0]
        idx_min = argrelextrema(singleosc, np.less)[0]
        osc_max = singleosc[idx_max]
        osc_min = singleosc[idx_min]
        n_cycles = 9
        amplitude = osc_max[-n_cycles:] - osc_min[-n_cycles:]
        mean = amplitude.mean()
        std = amplitude.std()

        return np.array(amplitude), mean, std

    def determine_amplitude_singleosc(self, t, singleosc):        
        """
        Determination of the amplitude of a single oscillator:
        peak-to-mean distance (for oscillations whose minimum
        is not evident, like some variables in the Almeida model that
        have spike-like dynamics)
        """
        idx_max = argrelextrema(singleosc, np.greater)[0]
        osc_max = singleosc[idx_max]
        n_cycles = 9

        # peak-to-mean
        dt = np.diff(t).mean()    
        mean_magnitude = singleosc[idx_max[-2]:idx_max[-1]].mean()
        amplitude = osc_max[-n_cycles:] - mean_magnitude
        mean = amplitude.mean()
        std = amplitude.std()

        return np.array(amplitude), mean, std

    def amplitudes(self, t, x):
        """
        Determine amplitudes of N oscillators 
        and return mean and std 
        """
        results = []
        for so in range(np.shape(x)[1]):
            amplitude = self.determine_amplitude_singleosc_peaktrough(
                t,x[:,so]
                )[1].tolist()
            results.append(amplitude)
        mean = np.mean(results)
        std  = np.std(results)
        return np.array(results), mean, std            