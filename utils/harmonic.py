import numpy as np

class HarmonicOscillator:

    def __init__ (self, k, m, gamma, beta=False, F=False, T=False):
        """
        Parameters of a Harmonic oscillator
        If beta is given, then the harmonic oscillator ODE
        becomes nonlinear and this particular case of the 
        harmonic oscillator is called the Duffing oscillator
        """
        self.k = k #constant form Hook's Law, controls linear stiffness
        self.m = m #mass
        self.gamma = gamma #friction
        self.beta = beta #if beta = True we have the Duffing oscillator 
                         #controls the amout of nonlinearity in restoring F
        self.F = F #amplitude of driving force  
        self.T = T #period of driving force
        self.dimension = 2

    def driving_force(self, t):
        omega = 2*np.pi/self.T
        driv_force = self.F * np.cos(omega * t)
        return driv_force
        

    def dynamics (self, y0, t): 
        """
        Second-order linear ODE of a harmonic/Duffing oscillator:
         m*(x_dotdot) + gamma*(x_dot) + kx (+ beta*x**3) = 0
        can be converted into a linear system of 2 first order ODEs
        by doing the change of variable v = x_dot
        (_dot represents the derivative with respect to time)   
        """ 
        f = np.zeros_like(y0)
        x = y0[0::2]
        v = y0[1::2]

        fx = f[0::2]
        fv = f[1::2]
        
        driving_force = 0
        if (self.F != 0.0) or (self.F != False):
            driving_force = self.driving_force(t)

        fx[:] = v
        fv[:] = -(self.k/self.m)*x - (self.gamma/self.m)*v - \
                (self.beta/self.m)*x**3 - driving_force
        
        return f

    def potential(self, y):
        """
        Potential of a harmonic/Duffing oscillator:
        See calculation here
        https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/chaos01.htm
        """
        V = 0.5 * self.k * y**2 + 0.25 * self.beta * y**4
        return V