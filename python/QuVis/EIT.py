from OBE import *

class EIT(OpticalBlochEquation):
    def __init__(self,Omega1,Omega2,Delta,delta,Gamma1,Gamma2,tmax,init_state,DeltaRange):
        OpticalBlochEquation.__init__(self,Omega1,Omega2,Delta,delta,Gamma1,Gamma2,tmax,init_state)
        # EIT state
        self.state_EIT = np.mat(init_state.reshape(9,1)) #matrix array
        #define detuning step
        self.numd = 1000
        self.Dmax = DeltaRange
        self.detuning = np.linspace(-1*self.Dmax,self.Dmax,self.numd)


    def LBsuperoperator_EIT(self,Delta=None):
        if Delta is None:
            Delta = self.Delta
        H_a = Delta*self.g1.T*self.g1 + (Delta-self.delta)*self.g2.T*self.g2
        H_af = 0.5*self.Omega1*(self.sigma1 + np.conj(self.sigma1.T)) + 0.5*self.Omega2*(self.sigma2 + np.conj(self.sigma2.T))
        H = H_a + H_af
        #H = 0.5*np.array([[0,0,self.Omega1],[0,-2*Delta,self.Omega2],[np.conj(self.Omega1),np.conj(self.Omega2),-2*(Delta-self.delta)]])
        #define superoperator
        H_eff = -1j*np.mat(kron(self.I3,H) - kron(np.conj(H.T),self.I3))
        L_eff = self.Gamma1*Dissipator(self.sigma1) + self.Gamma2*Dissipator(self.sigma2)
        S = H_eff+L_eff
        return S

    def LBDiagonalise_EIT(self,Delta=None):
        if Delta is None:
            Delta = self.Delta
        #diagonalize Hamiltonian
        evals, evecs = eig(self.LBsuperoperator_EIT(Delta))
        evecs = np.mat(evecs)
        return evals, evecs

    def getNextState_EIT(self,Delta=None):
        if Delta is None:
            Delta = self.Delta
        self.state_EIT = self.LBDiagonalise_EIT(Delta)[1]*np.mat(np.diag(np.exp(self.LBDiagonalise_EIT(Delta)[0]*self.tmax)))*np.linalg.inv(self.LBDiagonalise_EIT(Delta)[1])*self.init_state

    def Initialise_EIT(self):
        self.rhoba_r = np.zeros(1)
        self.rhoba_i = np.zeros(1)

    def saveSusceptibility(self):
        rhobar = np.real(self.state_EIT[2])
        rhobai = -1*np.imag(self.state_EIT[2])
        self.rhoba_r = np.append(self.rhoba_r,rhobar)
        self.rhoba_i = np.append(self.rhoba_i,rhobai)

    def Susceptibility(self):
        self.Initialise_EIT()
        for i in range(self.numd):
            self.getNextState_EIT(self.detuning[i])
            self.saveSusceptibility()
        self.rhoba_r = np.delete(self.rhoba_r,0)
        self.rhoba_i = np.delete(self.rhoba_i,0)
