import numpy as np

class GonzeOscillator:

    def __init__ (self, k1=0.7, k2=0.35, k3=0.7, k4=0.35, k5=0.7, k6=0.35, 
                  K1=1, K2=1, K4=1, K6=1, n=4, 
                  F=False, T=False,
                  F_perturbation = False, t_perturbation = False):
        """
        Parameters of a modified Goodwin oscillator (Gonze oscillator)
        with Michaelian degradation kinetics 
        Reference: Gonze D, Bernard S, Waltermann C, Kramer A, Herzel H. 
        Spontaneous synchronization of coupled circadian oscillators. 
        Biophysical journal. 2005;89(1):120-9.
        """
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.k4 = k4
        self.k5 = k5
        self.k6 = k6
        self.K1 = K1
        self.K2 = K2
        self.K4 = K4
        self.K6 = K6
        self.n = n
        self.F = F #strength of zeitgeber input
        self.T = T #period of zeitgeber input
        self.F_perturbation = F_perturbation
        self.t_perturbation = t_perturbation
        self.dimension = 3

    def zeitgeber (self, t):
        """
        zeitgeber sinusoidal input
        """
        zg = self.F * np.cos (2*np.pi/self.T * t + np.pi/2)
        return zg

    def square_perturbation (self, t):
        """
        Square like perturbation
        """
        amp = self.F_perturbation
        t_0 = self.t_perturbation 
        width = .5 #semi-width (duration will be 2*width)
        return amp * (abs(t - t_0) < width)

    def dynamics (self, y0, t): 
        """
        modified Goodwin oscillator (Gonze oscillator) with Michaelian 
        degradation kinetics 
        Reference: Gonze D, Bernard S, Waltermann C, Kramer A, Herzel H. 
        Spontaneous synchronization of coupled circadian oscillators. 
        Biophysical journal. 2005;89(1):120-9.
        """ 
        f = np.zeros_like(y0)
        x = y0[0::3]
        y = y0[1::3]
        z = y0[2::3]

        fx = f[0::3]
        fy = f[1::3]
        fz = f[2::3]

        # perturbation
        x_pert = 0
        if self.F_perturbation != False:
            x_pert = self.square_perturbation(t)
        
        # zeitgeber
        zg = 0
        if self.F != False:
            zg = self.zeitgeber(t)           

        # odes
        fx[:] = self.k1*(self.K1**self.n/(self.K1**self.n + z**self.n)) - \
                self.k2*x/(self.K2 + x) + zg + x_pert
        fy[:] = self.k3*x - self.k4*y/(self.K4 + y)
        fz[:] = self.k5*y - self.k6*z/(self.K6 + z)
        
        return f       


#################################################
#################################################


class GoodwinOscillator:

    def __init__ (self, k1=1.0, k2=0.20, k3=1.0, k4=0.15, k5=1.0, k6=0.10, 
                  K1=1, n=9.5, 
                  F=False, T=False,
                  F_perturbation = False, t_perturbation = False):
        """
        Parameters of a Goodwin oscillator with linear degradation
        Reference: Goodwin BC. Oscillatory behavior in enzymatic control 
        processes
        Advances in enzyme regulation. 1965;3:425-37.
        Parameters from:
        https://journals.sagepub.com/doi/pdf/10.1177/074873099129001037
        https://www.sciencedirect.com/science/article/pii/S0022519396900673?via%3Dihub
        but modified because default parameters produced damped rhythms
        """
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.k4 = k4
        self.k5 = k5
        self.k6 = k6
        self.K1 = K1
        self.n = n
        self.F = F #strength of zeitgeber input
        self.T = T #period of zeitgeber input
        self.F_perturbation = F_perturbation
        self.t_perturbation = t_perturbation
        self.dimension = 3

    def zeitgeber (self, t):
        """
        zeitgeber sinusoidal input
        """
        zg = self.F * np.cos (2*np.pi/self.T * t + np.pi/2)
        return zg

    def square_perturbation (self, t):
        """
        Square like perturbation
        """
        amp = self.F_perturbation
        t_0 = self.t_perturbation 
        width = .5 #semi-width (duration will be 2*width)
        return amp * (abs(t - t_0) < width)

    def dynamics (self, y0, t): 
        """
        Goodwin oscillator with linear degradation
        Ref: Goodwin BC. Oscillatory behavior in enzymatic control processes. 
        Advances in enzyme regulation. 1965;3:425-37.
        """ 
        f = np.zeros_like(y0)
        x = y0[0::3]
        y = y0[1::3]
        z = y0[2::3]

        fx = f[0::3]
        fy = f[1::3]
        fz = f[2::3]

        # perturbation
        x_pert = 0
        if self.F_perturbation != False:
            x_pert = self.square_perturbation(t)
        
        # zeitgeber
        zg = 0
        if self.F != False:
            zg = self.zeitgeber(t)           

        # odes
        fx[:] = self.k1*((self.K1**self.n)/(self.K1**self.n + z**self.n)) - \
                self.k2*x + zg + x_pert
        fy[:] = self.k3*x - self.k4*y
        fz[:] = self.k5*y - self.k6*z
        
        return f       
