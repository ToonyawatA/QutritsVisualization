from OBE import *
from STIRAP import *
from EIT import *

#initial init_state
psi0  = np.zeros(9)
psi0[0] = 1.0

#initiate project
#project = OpticalBlochEquation(Omega1=6.76,Omega2=2.36,Delta=0.01,delta=19.21,Gamma1=0.2,Gamma2=0.25,tmax=1.81,init_state=psi0)

#Calculating it's trajectory
#project.Trajectory()

#graph potting
#project.makePlot(population=True)

#pyvista plotting
#project.makePyvista()


###


#initial project2 STIRAP
#project2 = STIRAP(Omega1=2.51,Omega2=2.51,Delta=0.00,delta=0.00,t1=32.25,t2=20.4,tau1=5.91,tau2=5.91,tmax=53.95,init_state=psi0)

#Calculating it's trajectory
#project2.Trajectory()

#graph potting
#project2.makePlot(rabi=True)

#pyvista plotting
#project2.makePyvista()


###


#initiate project3 EIT
project3 = EIT(Omega1=0.2,Omega2=0.5,Gamma1=1.0,Gamma2=0.5,delta=0.0,Delta=0.0,tmax=40.0,init_state=psi0,DeltaRange=5.0)

#Calculating it's trajectory
project3.Trajectory()

#Calculating it's Susceptibility
project3.Susceptibility()

#graph potting
project3.makePlot(population=True)

#print(project3.probability)
#graph potting
project3.makePlot(susceptibility=True)
