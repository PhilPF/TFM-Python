# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 19:23:50 2022

@author: ppf
"""

import numpy as np
from Jet import Jet

class Explicit:
    
    t=[]; jets=[]
    
    def __init__(self, *, h, t_0, t_F, jet_0, f):
        self.h=h
        self.t_0=t_0
        self.t_F=t_F
        self.jet_0=jet_0
        self.f=f
        
    def __iterate(self, method):
        if isinstance(self.t_0, (float, np.floating, int)):
            self.t.append(self.t_0)
        if isinstance(self.jet_0, Jet):
            self.jets.append(self.jet_0)
                        
        self.t.append(self.t[-1]+self.h)
        while self.t[-1]<=self.t_F:
            self.jets.append(method(self.h, self.t[-self.r-1:-1], self.jets[-self.r:]))
            self.t.append(self.t[-1]+self.h)
         
        self.t.pop()
        #self.t[-1]=self.t_F
        #self.jets.append(method(self.t_F-self.t[-2], f, self.t[-1], self.jets[-1]))
                
        return len(self.t)-1,self.t, np.array([d.getJet() for d in self.jets]).transpose()
     
    a=[];b=[];c=[];s=0
    def __RK(self, h, t, y):
        if not isinstance(t, (float, np.floating, int)): t=t[0]
        if not isinstance(y, Jet): y=y[0]
        k=[]
        k.append(self.f(t,y))
        for i in range(1,self.s):
            k.append(self.f(t+h*sum(self.c[1:i+1]),y+h*Jet.dot(self.a[i,0:i],k[0:i])))
            
        return y+h*Jet.dot(self.b,k)  
    
    def Euler(self):
        self.s=1
        self.a=np.zeros((self.s,self.s))
        # self.b=np.zeros(self.s)
        self.c=np.zeros(self.s)
        
        self.b=[1]
        
        return self.__iterate(self.__RK)
    
    def RK4(self):
        self.s=4
        self.a=np.zeros((self.s,self.s))
        # self.b=np.zeros(self.s)
        # self.c=np.zeros(self.s)
        
        self.a[1,0]=0.5
        self.a[2,1]=0.5
        self.a[3,2]=1
        self.b=[1/6,1/3,1/3,1/6]
        self.c=[0,1/2,1/2,1]
        
        return self.__iterate(self.__RK)
    
    u=[];v=[];r=1
    def __multistep(self, h, t, y):
        k=[]
        for i in range(self.r):
            k.append(self.f(t[i],y[i]))
        return -Jet.dot(self.u,y)+h*Jet.dot(self.v, k)
    
    def Euler_AB1(self):
        self.r=1
        self.u =[-1]
        self.v=[1]
        return self.__iterate(self.__multistep)
    
    def AB2(self):
        self.r=2
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0).Euler()
        self.u =[0,-1]
        self.v=[-1/2,3/2]
        return self.__iterate(self.__multistep)
    
    def AB3(self):
        self.r=3
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0).Euler()
        self.u =[0,0,-1]
        self.v=[5/12,-16/12,23/12]
        return self.__iterate(self.__multistep)
    
    def AB4(self):
        self.r=4
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0).Euler()
        self.u =[0,0,0,-1]
        self.v=[-9/24,37/24,-59/24,55/24]
        return self.__iterate(self.__multistep)
    
    def AB5(self):
        self.r=5
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0).Euler()
        self.u =[0,0,0,0,-1]
        self.v=[251/720,-1274/720,2616/720,-2774/720,1901/720]
        return self.__iterate(self.__multistep)
    
    def AB6(self):
        self.r=6
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0).Euler()
        self.u =[0,0,0,0,0,-1]
        self.v=[-475/1440,2877/1440,-7298/1440,9982/1440,-7923/1440, 4277/1440]
        return self.__iterate(self.__multistep)

class Implicit:
    
    t=[]; jets=[]
    
    def __init__(self, *, h, t_0, t_F, jet_0, f):
        self.h=h
        self.t_0=t_0
        self.t_F=t_F
        self.jet_0=jet_0
        self.f=f
        
    def __solve_implicit(self,method,t_0,y_0):
        
        print("Hey")        
        
        #Solve first the zero-order jet case
        y0_order0 = Jet([y_0.getJet()[0]])
        y=[y0_order0 for m in range(self.s)]
        
        print(y)
                
        Newt_it=1
        for i in range(Newt_it):
            method_y=method(k_0=y_0,k=y)
            print(method_y)
            multi_dmethod_y=[a-b for a,b in zip(method(k_0=y_0,k=[Jet([m.getJet()[0],1]) for m in y]),method_y)]  #Notice that F(y+s)=F(y)+F'(y_1)s
            dmethod_y=[Jet([m.getJet()[1]]) for m in multi_dmethod_y]
            y=[a-b for a,b in zip(y,[n/m for n,m in zip(method_y,dmethod_y)])]
            
        print(y)
        
        #Solve the general N-order jet case
        for g in range(1,y_0.getOrder()):
            
            #Notice F(a+y_0^[N]s^N, b)=F(a,b)+D_{y_0}F(a,b)y_0^[N]s^N
            #Notice F(a, b+s^N)=F(a,b)+D_{y}F(a,b)s^N
            #Notice y^[N]=-D_{y_0}F(a,b)y_0^[N]/D_{y}F(a,b)
                    
            y0_ncoeff=y_0.getJet()[g]
            multi_D_y0F_y0=[a-b for a,b in zip(method(k_0=Jet([*y0_order0.getJet(),y0_ncoeff]),k=y) , method(k_0=y_0,k=y))]
            D_y0F_y0=[m.getJet()[g] for m in multi_D_y0F_y0]
            multi_D_yF=[a-b for a,b in zip(method(k_0=y_0,k=[Jet([*m.getJet(),1]) for m in y]) , method(k_0=y_0,k=y)) ]
            D_yF=[m.getJet()[g] for m in multi_D_yF]
            
            y_n=[-n/m for n,m in zip(D_y0F_y0,D_yF)]
            y=[Jet([*n.getJet(),m]) for n,m in zip(y,y_n)]
                        
            y0_order0=Jet(y_0.getJet()[0:g]) 
                       
        return y
        
    def __iterate(self, method):
        if isinstance(self.t_0, (float, np.floating, int)):
            self.t.append(self.t_0)
        if isinstance(self.jet_0, Jet):
            self.jets.append(self.jet_0)
                        
        self.t.append(self.t[-1]+self.h)
        while self.t[-1]<=self.t_F:
            self.jets.append(method(self.t[-1], self.jets[-1]))
            self.t.append(self.t[-1]+self.h)
         
        self.t.pop()
        #self.t[-1]=self.t_F
        #self.jets.append(method(self.t_F-self.t[-2], f, self.t[-1], self.jets[-1]))
                
        return len(self.t)-1,self.t, np.array([d.getJet() for d in self.jets]).transpose()
     
    a=[];b=[];c=[];s=0
    
    def __implicit_RK(self, k_0, k):
            
        if self.s==1:
            return [self.f(self.t_0+self.h*self.c[0],k_0+self.h*self.a[0]*k[0])-k[0]]
        else:
            return [self.f(self.t_0+self.h*self.c[i],k_0+self.h*Jet.dot(self.a[i,0:self.s],k[0:self.s]))-k[i] for i in range(self.s)]
    
    def __RK(self, t, y):
        if not isinstance(t, (float, np.floating, int)): t=t[0]
        if not isinstance(y, Jet): y=y[0]
        
        k=self.__solve_implicit(self.__implicit_RK,t,y)
            
        return y+self.h*Jet.dot(self.b,k)  
    
    def Euler(self):
        self.s=1
        self.a=[1]
        self.b=[1]
        self.c=[1]
        
        return self.__iterate(self.__RK)
    
    def RK4(self):
        self.s=4
        self.a=np.zeros((self.s,self.s))
        # self.b=np.zeros(self.s)
        # self.c=np.zeros(self.s)
        
        self.a[1,0]=0.5
        self.a[2,1]=0.5
        self.a[3,2]=1
        self.b=[1/6,1/3,1/3,1/6]
        self.c=[0,1/2,1/2,1]
        
        return self.__iterate(self.__RK)
    
    u=[];v=[];r=1
    def __multistep(self, h, t, y):
        k=[]
        for i in range(self.r):
            k.append(self.f(t[i],y[i]))
        return -Jet.dot(self.u,y)+h*Jet.dot(self.v, k)
    
    def Euler_AB1(self):
        self.r=1
        self.u =[-1]
        self.v=[1]
        return self.__iterate(self.__multistep)
    
    def AB2(self):
        self.r=2
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0).Euler()
        self.u =[0,-1]
        self.v=[-1/2,3/2]
        return self.__iterate(self.__multistep)
    
    def AB3(self):
        self.r=3
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0).Euler()
        self.u =[0,0,-1]
        self.v=[5/12,-16/12,23/12]
        return self.__iterate(self.__multistep)
    
    def AB4(self):
        self.r=4
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0).Euler()
        self.u =[0,0,0,-1]
        self.v=[-9/24,37/24,-59/24,55/24]
        return self.__iterate(self.__multistep)
    
    def AB5(self):
        self.r=5
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0).Euler()
        self.u =[0,0,0,0,-1]
        self.v=[251/720,-1274/720,2616/720,-2774/720,1901/720]
        return self.__iterate(self.__multistep)
    
    def AB6(self):
        self.r=6
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0).Euler()
        self.u =[0,0,0,0,0,-1]
        self.v=[-475/1440,2877/1440,-7298/1440,9982/1440,-7923/1440, 4277/1440]
        return self.__iterate(self.__multistep)
 