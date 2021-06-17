from OBE import *

class STIRAP(OpticalBlochEquation):
    def __init__(self,Omega1,Omega2,Delta,delta,t1,t2,tau1,tau2,tmax,init_state,Gamma1=0,Gamma2=0):
        OpticalBlochEquation.__init__(self,Omega1,Omega2,Delta,delta,Gamma1,Gamma2,tmax,init_state)
        self.t1 = t1
        self.t2 = t2
        self.tau1 = tau1
        self.tau2 = tau2
        self.omega1 = np.abs(np.real(self.Omega1*np.exp((-1/(2*self.tau1**2))*(self.time-self.t1)**2)))
        self.omega2 = np.abs(np.real(self.Omega2*np.exp((-1/(2*self.tau2**2))*(self.time-self.t2)**2)))

    def LBsuperoperator(self,omega1,omega2):
        H_a = self.Delta*self.g1.T*self.g1 + (self.Delta-self.delta)*self.g2.T*self.g2
        H_af = 0.5*omega1*(self.sigma1 + np.conj(self.sigma1.T)) + 0.5*omega2*(self.sigma2 + np.conj(self.sigma2.T))
        H = H_a + H_af
        #define superoperator
        H_eff = -1j*np.mat(kron(self.I3,H) - kron(np.conj(H.T),self.I3))
        L_eff = self.Gamma1*Dissipator(self.sigma1) + self.Gamma2*Dissipator(self.sigma2)
        S = H_eff+L_eff
        return S

    def LBDiagonalise(self,omega1,omega2):
        #define Hamiltonian
        evals, evecs = eig(self.LBsuperoperator(omega1,omega2))
        evecs = np.mat(evecs)
        return evals, evecs

    def getNextState(self,omega1,omega2):
        self.state = self.LBDiagonalise(omega1,omega2)[1]*np.mat(np.diag(np.exp(self.LBDiagonalise(omega1,omega2)[0]*self.dt)))*np.linalg.inv(self.LBDiagonalise(omega1,omega2)[1])*self.state

    def Trajectory(self):
        #define time array
        self.Initialise()
        for i in range(self.numt):
            self.saveTrajectory()
            self.getNextState(self.omega1[i],self.omega2[i])
        self.bloch1 = self.radius*self.bloch1
        self.bloch2 = self.radius*self.bloch2
        self.bloch1 = np.delete(self.bloch1,0,0)
        self.bloch2 = np.delete(self.bloch2,0,0)
        self.probability = np.delete(self.probability,0,0)
