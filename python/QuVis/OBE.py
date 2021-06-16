import numpy as np
from numpy import pi
import matplotlib
import matplotlib.pyplot as plt
import pyvista as pv
from pyvista import examples
from scipy.linalg import kron, eig
import time
import matplotlib.colors as mcolors
matplotlib.rcParams['text.usetex'] = True
from string import ascii_lowercase

def Dissipator(sigma):
    d = np.shape(sigma)[0]
    sigma = np.mat(sigma)
    I = np.eye(d)
    I = np.mat(I)
    return np.mat(kron(np.conj(sigma),sigma) - 0.5*kron(I,np.conj(sigma.T)*sigma) - 0.5*kron((np.conj(sigma.T)*sigma).T,I))


class OpticalBlochEquation:
    ### init ###
    def __init__(self,Omega1,Omega2,Delta,delta,Gamma1,Gamma2,tmax,init_state):
        self.Omega1 = Omega1
        self.Omega2 = Omega2
        self.Delta = Delta
        self.delta = delta
        self.Gamma1 = Gamma1
        self.Gamma2 = Gamma2
        self.tmax = tmax
        self.init_state = np.mat(init_state.reshape(9,1)) #matrix array
        self.state = np.mat(init_state.reshape(9,1)) #matrix array
        #define radius
        self.radius = 3.0
        #define quantum state
        self.g1 = np.mat(np.array([1,0,0]))
        self.g2 = np.mat(np.array([0,1,0]))
        self.e = np.mat(np.array([0,0,1]))
        #define matrices
        self.sigma1 = self.g1.T*self.e
        self.sigma2 = self.g2.T*self.e
        self.sigmax = np.mat(np.array([[0,1,0],[1,0,0],[0,0,0]]))
        self.sigmay = np.mat(np.array([[0,-1j,0],[1j,0,0],[0,0,0]]))
        self.sigmaz = np.mat(np.array([[1,0,0],[0,-1,0],[0,0,0]]))
        self.sigma0 = np.mat(np.array([[0,0,0],[0,0,0],[0,0,1]]))
        self.I3 = np.mat(np.eye(3))
        #define time step and
        self.numt = 1000
        self.time = np.linspace(0,self.tmax,self.numt)
        self.dt = self.time[1]-self.time[0]
        #define detuning step
        self.numd = 1000

    def Initialise(self):
        pg1 = np.real(self.init_state[0,0])
        pg2 = np.real(self.init_state[4,0])
        pe = np.real(self.init_state[8,0])
        #density matrix
        DM = self.init_state.reshape(3,3)
        #bloch1
        u1 = np.real(np.trace(self.sigmax*DM))
        v1 = np.real(np.trace(self.sigmay*DM))
        w1 = np.real(np.trace(self.sigmaz*DM))
        #bloch2
        r2 = np.real(np.trace(self.sigma0*DM))
        if(r2==1):
            u2 = r2*np.cos(0.5*pi*r2)
            v2 = r2*np.cos(0.5*pi*r2)
            w2 = r2*np.sin(0.5*pi*r2)
        else:
            u2 = r2*np.cos(0.5*pi*r2)*(np.real(self.init_state[2,0])/np.abs(self.init_state[2,0])) if np.abs(self.init_state[2,0])!=0 else 0
            v2 = r2*np.cos(0.5*pi*r2)*(np.imag(self.init_state[2,0])/np.abs(self.init_state[2,0])) if np.abs(self.init_state[2,0])!=0 else 0
            w2 = r2*np.sin(0.5*pi*r2) if np.abs(self.init_state[2,0])!=0 else 0
        self.probability = np.array([pg1,pg2,pe])
        self.bloch1 = np.array([u1,v1,w1])
        self.bloch2 = np.array([u2,v2,w2])

    def saveTrajectory(self):
        #state probability
        pg1 = np.real(self.state[0,0])
        pg2 = np.real(self.state[4,0])
        pe = np.real(self.state[8,0])
        #density matrix
        DM = self.state.reshape(3,3)
        #bloch1
        u1 = np.real(np.trace(self.sigmax*DM))
        v1 = np.real(np.trace(self.sigmay*DM))
        w1 = np.real(np.trace(self.sigmaz*DM))
        #bloch2
        r2 = np.real(np.trace(self.sigma0*DM))
        if(r2==1):
            u2 = r2*np.cos(0.5*pi*r2)
            v2 = r2*np.cos(0.5*pi*r2)
            w2 = r2*np.sin(0.5*pi*r2)
        else:
            u2 = r2*np.cos(0.5*pi*r2)*(np.real(self.state[2,0])/np.abs(self.state[2,0])) if np.abs(self.state[2,0])!=0 else 0
            v2 = r2*np.cos(0.5*pi*r2)*(np.imag(self.state[2,0])/np.abs(self.state[2,0])) if np.abs(self.state[2,0])!=0 else 0
            w2 = r2*np.sin(0.5*pi*r2) if np.abs(self.state[2,0])!=0 else 0
        parray = np.array([pg1,pg2,pe])
        b1array = np.array([u1,v1,w1])
        b2array = np.array([u2,v2,w2])
        self.probability = np.vstack((self.probability,parray))
        self.bloch1 = np.vstack((self.bloch1,b1array))
        self.bloch2 = np.vstack((self.bloch2,b2array))

    def LBsuperoperator(self):
        H_a = self.Delta*self.g1.T*self.g1 + (self.Delta-self.delta)*self.g2.T*self.g2
        H_af = 0.5*self.Omega1*(self.sigma1 + np.conj(self.sigma1.T)) + 0.5*self.Omega2*(self.sigma2 + np.conj(self.sigma2.T))
        H = H_a + H_af
        #define superoperator
        H_eff = -1j*np.mat(kron(self.I3,H) - kron(np.conj(H.T),self.I3))
        L_eff = self.Gamma1*Dissipator(self.sigma1) + self.Gamma2*Dissipator(self.sigma2)
        S = H_eff+L_eff
        return S

    def LBDiagonalise(self):
        #define Hamiltonian
        evals, evecs = eig(self.LBsuperoperator())
        evecs = np.mat(evecs)
        return evals, evecs

    def getNextState(self):
        self.state = self.LBDiagonalise()[1]*np.mat(np.diag(np.exp(self.LBDiagonalise()[0]*self.dt)))*np.linalg.inv(self.LBDiagonalise()[1])*self.state

    def Trajectory(self):
        #define time array
        self.Initialise()
        for i in range(self.numt-1):
            self.getNextState()
            self.saveTrajectory()
        self.bloch1 = self.radius*self.bloch1
        self.bloch2 = self.radius*self.bloch2

    def makePlot(self,arr,population=False,susceptibily=False,rabi=False):
        row, column = np.shape(arr)
        poplabel = [r'$|g_1\rangle$',r'$|g_2\rangle$',r'$|e\rangle$']
        plt.figure()
        if(population):
            for i in range(column):
                plt.plot(self.time,arr[:,i],label=poplabel[i])
            plt.xlabel(r'Time ($t$)')
            plt.ylabel(r'Population')
            plt.legend()
        elif(susceptibily):
            for i in range(column):
                plt.plot(self.time,arr[:,i],label=poplabel[i])
            plt.xlabel(r'Detuning ($\Delta$)')
            plt.ylabel(r'Susceptibility ($\chi$)')
            plt.legend()
        elif(rabi):
            for i in range(column):
                plt.plot(self.time,arr[:,i],label=poplabel[i])
            plt.xlabel(r'Time ($t$)')
            plt.ylabel(r'Rabi Frequency')
            plt.legend()
        plt.show()
