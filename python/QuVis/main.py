from OBE import *

#initial init_state
psi0  = np.zeros(9)
psi0[0] = 1.0

#initiate project
project = OpticalBlochEquation(Omega1=6.76,Omega2=2.36,Delta=0.01,delta=19.21,Gamma1=0.25,Gamma2=0.25,tmax=1.81,init_state=psi0)

#Calculating it's trajectory
project.Trajectory()

#graph potting
project.makePlot(project.probability,population=True)
