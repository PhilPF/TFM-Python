# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 20:36:49 2022

@author: ppf
"""

import numpy as np
import matplotlib.pyplot as plt
from Jet import Jet

def EulerE(h,f,y):
    return y+f(y)*h

def f(y):
    return -y

t_0=0
jet_0=Jet([1.0,0.5,0.1])

t_F=10

t=[]
t.append(t_0)

jet=[]
jet.append(jet_0)

h=0.1
t.append(t[-1]+h)
while t[-1]<t_F:
    jet.append(EulerE(h,f,jet[-1]))
    t.append(t[-1]+h)

t[-1]=t_F
jet.append(EulerE(t_F-t[-2],f,jet[-1]))

sol=np.array([d.getJet() for d in jet]).transpose()

for s in sol:
    plt.plot(t,s)
    
print('steps: '+str(len(t)))
