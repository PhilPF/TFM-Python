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
        if isinstance(self.t_0, int) or isinstance(self.t_0, float):
            self.t.append(self.t_0)
        if isinstance(self.jet_0, Jet):
            self.jets.append(self.jet_0)
                        
        self.t.append(self.t[-1]+self.h)
        while self.t[-1]<=self.t_F:
            self.jets.append(method(self.h, f, self.t[-self.r-1:-1], self.jets[-self.r:]))
            self.t.append(self.t[-1]+self.h)
         
        self.t.pop()
        #self.t[-1]=self.t_F
        #self.jets.append(method(self.t_F-self.t[-2], f, self.t[-1], self.jets[-1]))
                
        return len(self.t)-1,self.t, np.array([d.getJet() for d in self.jets]).transpose()
     
    a=[];b=[];c=[];s=0
    def __RK(self, h, f, t, y):
        if not (isinstance(t, int) or isinstance(t, float)): t=t[0]
        if not isinstance(y, Jet): y=y[0]
        k=[]
        k.append(f(t,y))
        for i in range(1,self.s):
            k.append(f(t+h*sum(self.c[1:i+1]),y+h*Jet.dot(self.a[i,0:i],k[0:i])))
            
        return y+h*Jet.dot(self.b,k)  
    
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
    
    u=[];v=[];r=1
    def __multistep(self, h, f, t, y):
        k=[]
        for i in range(self.r):
            k.append(f(t[i],y[i]))
        return -Jet.dot(self.u,y)+h*Jet.dot(self.v, k)
    
    def Euler_AB1(self, f):
        self.r=1
        self.u =[-1]
        self.v=[1]
        return self.__iterate(f, self.__multistep)
    
    def AB2(self, f):
        self.r=2
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0).Euler(f)
        self.u =[0,-1]
        self.v=[-1/2,3/2]
        return self.__iterate(f, self.__multistep)
    
    def AB3(self, f):
        self.r=3
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0).Euler(f)
        self.u =[0,0,-1]
        self.v=[5/12,-16/12,23/12]
        return self.__iterate(f, self.__multistep)
    
    def AB4(self, f):
        self.r=4
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0).Euler(f)
        self.u =[0,0,0,-1]
        self.v=[-9/24,37/24,-59/24,55/24]
        return self.__iterate(f, self.__multistep)
    
    def AB5(self, f):
        self.r=5
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0).Euler(f)
        self.u =[0,0,0,0,-1]
        self.v=[251/720,-1274/720,2616/720,-2774/720,1901/720]
        return self.__iterate(f, self.__multistep)
    
    def AB6(self, f):
        self.r=6
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0).Euler(f)
        self.u =[0,0,0,0,0,-1]
        self.v=[-475/1440,2877/1440,-7298/1440,9982/1440,-7923/1440, 4277/1440]
        return self.__iterate(f, self.__multistep)
    
class Implicit:
    
    t=[]; jets=[]
    
    def __init__(self, *, h, t_0, t_F, jet_0):
        self.h=h
        self.t_0=t_0
        self.t_F=t_F
        self.jet_0=jet_0
        
    def __iterate(self, f, method):
        if isinstance(self.t_0, int) or isinstance(self.t_0, float):
            self.t.append(self.t_0)
        if isinstance(self.jet_0, Jet):
            self.jets.append(self.jet_0)
                        
        self.t.append(self.t[-1]+self.h)
        while self.t[-1]<=self.t_F:
            self.jets.append(method(self.h, f, self.t[-self.r-1:-1], self.jets[-self.r:]))
            self.t.append(self.t[-1]+self.h)
         
        self.t.pop()
        #self.t[-1]=self.t_F
        #self.jets.append(method(self.t_F-self.t[-2], f, self.t[-1], self.jets[-1]))
                
        return len(self.t)-1,self.t, np.array([d.getJet() for d in self.jets]).transpose()
        
    def __Newton(self, *, F, f, t, y_0, tol=0.1):
        result=y_0
        N=10
        for i in range(N):
            result-=F(f,t,result,y_0)/Jet.naiveDerivative(result,f,y_0,F=F,t=t)
        return result
            
    def __euler(self,f,t,y,y_0):
        return y-y_0-self.h*f(t,y)
            
    def Euler(self, f):
        if isinstance(self.t_0, int) or isinstance(self.t_0, float):
            self.t.append(self.t_0)
        if isinstance(self.jet_0, Jet):
            self.jets.append(self.jet_0)
                        
        self.t.append(self.t[-1]+self.h)
        while self.t[-1]<=self.t_F:
            self.jets.append(self.__Newton(F=self.__euler, f=f, t=self.t[-1], y_0=self.jets[-1]))
            # self.jets.append(method(self.h, f, self.t[-self.r-1:-1], self.jets[-self.r:]))
            self.t.append(self.t[-1]+self.h)
         
        self.t.pop()
        #self.t[-1]=self.t_F
        #self.jets.append(method(self.t_F-self.t[-2], f, self.t[-1], self.jets[-1]))
                
        return len(self.t)-1,self.t, np.array([d.getJet() for d in self.jets]).transpose()
            
        
    