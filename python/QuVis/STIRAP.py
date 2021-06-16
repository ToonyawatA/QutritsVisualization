from OBE import *

class Stirap(OpticalBlochEquation):
    def __init__(self,Omega1,Omega2,Delta,delta,Gamma1,Gamm2,t1,t2,tau1,tau2,tmax,init_state):
        OpticalBlochEquation.__init__(self,Omega1,Omega2,Delta,delta,tmax,init_state,Gamma1=0,Gamma2=0)
        self.t1 = t1
        self.t2 = t2
        self.tau1 = tau1
        self.tau2 = tau2
