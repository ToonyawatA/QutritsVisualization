#!/usr/bin/env python
# coding: utf-8

# In[83]:


import numpy as np
import matplotlib.pyplot as plt
import time
import imageio
import os
import glob
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
from Ramsey_Sequence import Ramsey_Sequence


# In[84]:


class optical_bloch:
    def obe(bloch,t,omega,delta):
        gamma=0.0
        return np.array([-0.5*gamma*bloch[0]-delta*bloch[1],-0.5*gamma*bloch[1]+delta*bloch[0]-omega*bloch[2],-1.0*gamma*(bloch[2]-1)+omega*bloch[1]]) 

class halfpipulse:
    rabi = 74.0*np.pi
    detuning = 0.0
    time_interval = 0.5*np.pi/np.sqrt(rabi**2 + detuning**2)
    
class pipulse:
    rabi = 74.0*np.pi
    detuning = 0.0
    time_interval = np.pi/np.sqrt(rabi**2 + detuning**2)
    
class freepulse:
    rabi = 0.0
    detuning = 74.0*np.pi
    time_interval = 0.1
    
class Hadamard:
    rabi = 74.0*np.pi
    detuning = 74.0*np.pi
    time_interval = np.pi/np.sqrt(rabi**2 + detuning**2)
    


# In[85]:


pulse = Ramsey_Sequence(obe)
pulse.evolve(Hadamard)


# In[71]:


pulse.evolution


# In[ ]:




