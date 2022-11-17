# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 20:36:49 2022

@author: ppf
"""

import numpy as np
import matplotlib.pyplot as plt

class Jet:
    def __init__(self,jet):
        self.jet=jet
        self.order=len(jet)
        
    def __str__(self):
        s = ""
        for i in range(self.order):
            if i==0: s+=str(self.jet[i])
            elif i==1: s+='+'+str(self.jet[i])+'x'
            else: s+='+'+str(self.jet[i])+'x^'+str(i)
        return str(s)
    
    def __add__(self, other):
        if isinstance(other, Jet) and self.order==other.getOrder():
            result = np.zeros(self.order)
            for i in range(self.order):
                result[i]=self.jet[i]+other.getJet()[i]
            return Jet(result)
        # if isinstance(other, float):
        #     return Jet(other+self.jet)
        
    def __mul__(self, other):
        # if isinstance(other, Jet) and self.order==other.getOrder():
        #     result = np.zeros(self.order)
        #     for i in range(self.order):
        #         result[i]=self.jet[i]+other.getJet()[i]
        #     return Jet(result)
        if isinstance(other, float):
             return Jet([other*d for d in self.jet])
         
    def getJet(self): 
        return self.jet
    
    def getOrder(self): 
        return self.order

def EulerE(h,f,y):
    return y+f(y)*h

def f(y):
    return y*(-1.0)

t_0=0
jet_0=Jet([1.0,0.5])

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

sol=*np.array([d.getJet() for d in jet]).transpose()

for s in sol:
    plt.plot(t,s)
