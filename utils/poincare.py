import numpy as np

class PoincareOscillator:

    def __init__ (self, amp, lam, tau, eps, 
                  F=False, T=False,
                  F_perturbation = False, t_perturbation = False,
                  K = False):
        """
        Parameters of a modified Poincare oscillator with twist
        """
        self.amp = amp #amplitude of Poincare oscillator
        self.lam = lam #amplitude relaxation rate of Poincare osc.
        self.tau = tau #intrinsic period of Poincare oscillator
        self.eps = eps #twist
        self.F = F #strength of zeitgeber input
        self.T = T #period of zeitgeber input
        self.F_perturbation = F_perturbation #strength of perturbation
        self.t_perturbation = t_perturbation #time of perturbation start
        self.K = K #coupling strength
        self.dimension = 2

    def zeitgeber (self, t):
        """
        zeitgeber sinusoidal input
        """
        zg = self.F * np.cos (2*np.pi/self.T * t + np.pi/2)
        return zg

    def mean_field (self, y0):
        """
        coupling through mean field: the effective concentration of x
        can be approximated with the average level of all x_i signals 
        or mean field MF, which acts on individual oscillators at a 
        coupling strength K
        """
        mf = np.zeros_like(y0)
        xi, yi = y0[0::2], y0[1::2]
        mf_x, mf_y = mf[0::2], mf[1::2]

        n_oscs = xi.shape[0]
        MF     = np.sum(xi, axis=0) / n_oscs 

        meanfx = MF 
        mf_x[:] = meanfx
        return mf    

    def square_perturbation (self, t):
        amp = self.F_perturbation
        t_0 = self.t_perturbation 
        width = .5 #semi-width (duration will be 2*width)
        return amp * (abs(t - t_0) <= width)

    def dynamics_polar (self, y0, t): 
        """
        ODEs for deterministic Poincare system with twist,
        in polar coordinates. 
        """ 
        f = np.zeros_like(y0)
        r_i = y0[0::2]
        phi_i = y0[1::2]

        fr   = f[0::2]
        fphi = f[1::2]
            
        omega = 2*np.pi/self.tau
    
        fr[:] = self.lam * r_i * (self.amp - r_i) 
        fphi[:] = omega + self.eps * (self.amp - r_i)
        
        return f

    def dynamics_cartesian (self, y0, t): 
        """
        ODEs for deterministic Poincare system with twist,
        in Cartesian coordinates. 
        """ 
        f = np.zeros_like(y0)
        xi = y0[0::2]
        yi = y0[1::2]

        fx = f[0::2]
        fy = f[1::2]

        r       = np.sqrt(xi**2 + yi**2)
        phi_dot = 2*np.pi/self.tau + self.eps*(self.amp-r)
        
        # perturbation
        x_pert = 0
        if self.F_perturbation != False:
            x_pert = self.square_perturbation(t)
        
        # zeitgeber
        zg = 0
        if self.F != False:
            zg = self.zeitgeber(t)

        # coupling 
        mean_field = 0
        if self.K != False:
            mean_field = self.K * self.mean_field(y0)[0::2]

        # odes
        fx[:] = self.lam*xi*(self.amp-r) - yi*phi_dot + \
                zg + x_pert + mean_field
        fy[:] = self.lam*yi*(self.amp-r) + xi*phi_dot              
        
        return f

    def isochrones (self, r_latent):
        """
        define isochrones as in Winfree's book:
        Winfree, AT: The Geometry of Biological Time (1980, Springer)
        (see Materials & Methods in paper for analytical calculation)
        """
        phi_isoch = np.arange(0, 2*np.pi + np.pi/4, np.pi/4) #every 45deg

        phi_latent = []
        for i in phi_isoch: #iterate over the different r
            phi_lat = i - self.eps/self.lam*np.log(r_latent)
            phi_latent.append(phi_lat)
        phi_latent = np.asarray(phi_latent)

        return phi_latent