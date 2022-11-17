# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 20:36:49 2022

@author: ppf
"""

import numpy as np
import matplotlib.pyplot as plt
from Jet import Jet
from Method import Explicit

def f(y):
    return Jet.exp(-y**2)

steps, t, sols =  Explicit(h=0.5,t_0=0, t_F=10, jet_0=Jet([1.5,0.5])).Euler(f)

for sol in sols:
     plt.plot(t,sol)
     
print('steps: '+ str(steps))