# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 19:23:50 2022

@author: ppf
"""

import numpy as np

class Implicit:
    
    t=[]
    jets=[]
    
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
            self.jets.append(method(self.h, f, self.jets[-1]))
            self.t.append(self.t[-1]+self.h)
            
        self.t[-1]=self.t_F
        self.jets.append(method(self.t_F-self.t[-2], f, self.jets[-1]))
        
        return len(self.t)-2,self.t, np.array([d.getJet() for d in self.jets]).transpose()
        
    
    def Euler(self, f):
        return self.__iterate(f,self.__euler)
        
    def __euler(self,h,f,y):
        return y+f(y)*h
        