# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 19:23:50 2022

@author: ppf
"""

import numpy as np
from Jet import Jet

class Explicit:
    
    t=[]; jets=[]
    
    def __init__(self, *, h, t_0, t_F, jet_0):
        self.h=h
        self.t_0=t_0
        self.t_F=t_F
        self.jet_0=jet_0
        
    def __iterate(self, f, method):
        self.t.append(self.t_0)
        self.jets.append(self.jet_0)
        
        self.t.append(self.t[-1]+self.h)
        while self.t[-1]<self.t_F:
            self.jets.append(method(self.h, f, self.t[-1], self.jets[-1]))
            self.t.append(self.t[-1]+self.h)
            
        self.t[-1]=self.t_F
        self.jets.append(method(self.t_F-self.t[-2], f, self.t[-1], self.jets[-1]))
        
        return len(self.t)-2,self.t, np.array([d.getJet() for d in self.jets]).transpose()
        
    def Euler(self, f):
        self.s=1
        self.a=np.zeros((self.s,self.s))
        # self.b=np.zeros(self.s)
        self.c=np.zeros(self.s)
        
        self.b=[1]
        
        return self.__iterate(f,self.__RK)
    
    def RK4(self, f):
        self.s=4
        self.a=np.zeros((self.s,self.s))
        # self.b=np.zeros(self.s)
        # self.c=np.zeros(self.s)
        
        self.a[1,0]=0.5
        self.a[2,1]=0.5
        self.a[3,2]=1
        self.b=[1/6,1/3,1/3,1/6]
        self.c=[0,1/2,1/2,1]
        
        return self.__iterate(f,self.__RK)
    
    a=[];b=[];c=[];s=0
    def __RK(self, h, f, t, y):
        k=[]
        k.append(f(t,y))
        for i in range(1,self.s):
            k.append(f(t+h*sum(self.c[1:i+1]),y+h*Jet.dot(self.a[i,0:i],k[0:i])))
            
        return y+h*Jet.dot(self.b,k)  
    
    def AB2(self, f):
        extra_steps, self.t_0, self.jet_0 = Explicit(h=h,t_0=t_0,t_F=h+t_0,jet_0=jet_0).Euler(f)
        self.u=[0,1]
        self.v=[-1/2,3/2]
        return self.__iterate(f, self.__multistep)
    
    u=[];v=[];r=0
    def __multistep(self, h, f, t, y):
        k=[]
        for i in range(self.r-1):
            k.append(f(t[i],y[i]))
        return Jet.dot(self.u,y)-h*Jet.dot(self.v, k)
