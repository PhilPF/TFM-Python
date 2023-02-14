# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 19:23:50 2022

@author: ppf
"""

import numpy as np
from Jet import Jet

class Explicit:
        
    def __init__(self, *, h, t_0, t_F, jet_0, f):
        self.t=[]; self.jets=[]
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
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0,f=self.f).Euler()
        self.u =[0,-1]
        self.v=[-1/2,3/2]
        return self.__iterate(self.__multistep)
    
    def AB3(self):
        self.r=3
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0,f=self.f).Euler()
        self.u =[0,0,-1]
        self.v=[5/12,-16/12,23/12]
        return self.__iterate(self.__multistep)
    
    def AB4(self):
        self.r=4
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0,f=self.f).Euler()
        self.u =[0,0,0,-1]
        self.v=[-9/24,37/24,-59/24,55/24]
        return self.__iterate(self.__multistep)
    
    def AB5(self):
        self.r=5
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0,f=self.f).Euler()
        self.u =[0,0,0,0,-1]
        self.v=[251/720,-1274/720,2616/720,-2774/720,1901/720]
        return self.__iterate(self.__multistep)
    
    def AB6(self):
        self.r=6
        extra_step , self.t_0, self.jet_0 = Explicit(h=self.h,t_0=self.t_0,t_F=(self.r-1)*self.h+self.t_0,jet_0=self.jet_0,f=self.f).Euler()
        self.u =[0,0,0,0,0,-1]
        self.v=[-475/1440,2877/1440,-7298/1440,9982/1440,-7923/1440, 4277/1440]
        return self.__iterate(self.__multistep)

class Implicit:
        
    def __init__(self, exact_sol=None, *, h, t_0, t_F, jet_0, f, Newt_iter):
        self.t=[]; self.jets=[]
        self.h=h
        self.t_0=t_0
        self.t_F=t_F
        self.jet_0=jet_0
        self.f=f
        self.Newt_iter=Newt_iter
        
        self.exact_sol=exact_sol
        
    def __solve_implicit(self,method,t_0,y_0):
                
        #Solve first the zero-order jet case
        y0_order0 = Jet([y_0.getJet()[0]])
        
        y=[y0_order0 for m in range(self.s)]
        # y_approx_steps, y_approx_t, y_approx = Explicit(h=self.h,t_0=t_0,t_F=self.h+t_0,jet_0=y0_order0,f=self.f).Euler()
        # y=Jet([y_approx[0][1]])        
        
        if self.s==1:
            for i in range(self.Newt_iter):
                method_y=method(k_0=y0_order0,k=y)
                
                multi_dmethod_y=[a-b for a,b in zip(method(k_0=y0_order0,k=[Jet([m.getJet()[0],1]) for m in y]),method_y)]  #Notice that F(y+s)=F(y)+F'(y)s
                dmethod_y=[Jet([m.getJet()[1]]) for m in multi_dmethod_y]
                y=[a-b for a,b in zip(y,[n/m for n,m in zip(method_y,dmethod_y)])]
            
        else:
            for i in range(self.Newt_iter):
                method_y=method(k_0=y0_order0,k=y)
                
                Jac=[]
                for j in range(self.s):
                    k=y.copy()
                    k[j]=Jet([y[j].getJet()[0],1])
                    Jac.append([(n-m).getJet()[1] for n,m in zip(method(k_0=y0_order0, k=k),method_y)])
            
                #We solve a linear system to avoid matrix inverses UJac = Jac^-1*method_y
                Jac=np.array(Jac).transpose()
                UJac=np.linalg.solve(Jac,[m.getJet()[0] for m in method_y])

                y=[n-m for n,m in zip(y,UJac)]
        
        
        #Solve the general N-order jet case
        
        if self.s==1:
            for g in range(1,y_0.getOrder()):
                        
                y0_ncoeff=y_0.getJet()[g]
                D_y0F_y0=[(a-b).getJet()[g] for a,b in zip(method(k_0=Jet([*y0_order0.getJet(),y0_ncoeff]),k=y) , method(k_0=y0_order0,k=y))]
                    
                D_yF=[(a-b).getJet()[g] for a,b in zip(method(k_0=y0_order0,k=[Jet([*m.getJet(),1]) for m in y]) , method(k_0=y0_order0,k=y)) ]
                
                y_n=[-n/m for n,m in zip(D_y0F_y0,D_yF)]
                y=[Jet([*n.getJet(),m]) for n,m in zip(y,y_n)]
                            
                y0_order0=Jet(y_0.getJet()[0:g]) 
            
        else:
            for g in range(1,y_0.getOrder()):
                        
                y0_ncoeff=y_0.getJet()[g]
                D_y0F_y0=[(a-b).getJet()[g] for a,b in zip(method(k_0=Jet([*y0_order0.getJet(),y0_ncoeff]),k=y) , method(k_0=y0_order0,k=y))]
                
                Jac=[]
                for j in range(self.s):
                    k=y.copy()
                    k[j]=Jet([*y[j].getJet(),1])
                    Jac.append([(n-m).getJet()[g] for n,m in zip(method(k_0=y0_order0, k=k), method(k_0=y0_order0,k=y))])
     
                Jac = np.array(Jac).transpose()
                UJac=np.linalg.solve(Jac,D_y0F_y0)
                y_n=[-m for m in UJac]
                y=[Jet([*n.getJet(),m]) for n,m in zip(y,y_n)]
                            
                y0_order0=Jet(y_0.getJet()[0:g]) 
                       
        return y
        
    def __iterate(self, method):
        if isinstance(self.t_0, (float, np.floating, int)):
            self.t.append(self.t_0)
        if isinstance(self.jet_0, Jet):
            self.jets.append(self.jet_0)
            
        if self.exact_sol is not None:
            
            local_errors=[]
            global_errors=[]
            
            self.t.append(self.t[-1]+self.h)
            while self.t[-1]<=self.t_F:
                self.jets.append(method(self.t[-self.r-1:-1], self.jets[-self.r:]))
                
                exact = [sol(self.t[-1]) for sol in self.exact_sol]
                exact_last = [sol(self.t[-2]) for sol in self.exact_sol]
                exact_last_method=method(self.t[-2], Jet(exact_last))
                local_errors.append([exact[i]-exact_last_method.getJet()[i] for i in range(self.jets[-1].getOrder())])
                global_errors.append([self.jets[-1].getJet()[i]-exact[i] for i in range(self.jets[-1].getOrder())])
                
                self.t.append(self.t[-1]+self.h)

            self.t.pop()
            
            return len(self.t),self.t, np.array([d.getJet() for d in self.jets]).transpose(), np.array(local_errors).transpose(), np.array(global_errors).transpose()
            
        else:
            self.t.append(self.t[-1]+self.h)
            while self.t[-1]<=self.t_F:
                self.jets.append(method(self.t[-self.r-1:-1], self.jets[-self.r:]))
                self.t.append(self.t[-1]+self.h)
             
            self.t.pop()
            #self.t[-1]=self.t_F
            #self.jets.append(method(self.t_F-self.t[-2], f, self.t[-1], self.jets[-1]))
                
            return len(self.t),self.t, np.array([d.getJet() for d in self.jets]).transpose()
     
    a=[];b=[];c=[];s=0
    
    def __implicit_RK(self, k_0, k):
        if self.s==1:
            return [self.f(self.t[-1]+self.h*self.c[0],k_0+self.h*self.a[0]*k[0])-k[0]]
        else:
            return [self.f(self.t[-1]+self.h*self.c[i],k_0+self.h*Jet.dot(self.a[i,0:self.s],k[0:self.s]))-k[i] for i in range(self.s)]
    
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
    
    def Midpoint(self):
        self.s=1
        self.a=[1/2]
        self.b=[1]
        self.c=[1/2]
        
        return self.__iterate(self.__RK)
    
    def Crank_Nicolson(self):
        self.s=2
        self.a=np.zeros((self.s,self.s))
        self.a[1,0]=0.5
        self.a[1,1]=0.5
        self.b=[0.5,0.5]
        self.c=[0,1]
        
        return self.__iterate(self.__RK)
    
    def KS_DIRK(self):
        self.s=2
        self.a=np.array([[1/2,0],[-1/2,2]])
        self.b=[-1/2,3/2]
        self.c=[1/2,3/2]
        
        return self.__iterate(self.__RK)
    
    def QZ_DIRK(self):
        self.s=2
        self.a=np.array([[1/4,0],[1/2,1/4]])
        self.b=[1/2,1/2]
        self.c=[1/4,3/4]
        
        return self.__iterate(self.__RK)
    
    def Crouzeix_DIRK(self):
        self.s=2
        self.a=np.array([[1/2+np.sqrt(3)/6,0],[-np.sqrt(3)/3,1/2+np.sqrt(3)/6]])
        self.b=[1/2,1/2]
        self.c=[1/2+np.sqrt(3)/6,1/2-np.sqrt(3)/6]
        
        return self.__iterate(self.__RK)
    
    def RADAU_IA_3o(self):
        self.s=2
        self.a=np.array([[1/4,-1/4],[1/4,5/12]])
        self.b=[1/4,3/4]
        self.c=[0,2/3]
        
        return self.__iterate(self.__RK)
    
    def RADAU_IA_5o(self):
        sqrt6= np.sqrt(6)
        self.s=3
        self.a=np.array([[1/9,(-1-sqrt6)/18,(-1+sqrt6)/18],[1/9,11/45+7*sqrt6/360, 11/45-43*sqrt6/360],[1/9,11/45+43*sqrt6/360, 11/45-7*sqrt6/360]])
        self.b=[1/9,4/9+sqrt6/36,4/9-sqrt6/36]
        self.c=[0,3/5-sqrt6/10,3/5+sqrt6/10]
        
        return self.__iterate(self.__RK)
    
    def RADAU_IIA_3o(self):
        self.s=2
        self.a=np.array([[5/12,-1/12],[3/4,1/4]])
        self.b=[3/4,1/4]
        self.c=[1/3,1]
        
        return self.__iterate(self.__RK)
    


    u=[];v=[];r=1
    
    def __implicit_multistep(self, k_0, k):
        return [self.h*(self.v[-1]*self.f(self.t[-1],k)+Jet.dot(self.v[0:self.r-1],[self.f(self.t[i],k_0[i]) for i in range(self.r-1)]))-(self.u[-1]*k+Jet.dot(self.u[0:self.r-1],k_0))]
        
    def __multistep(self, t, y):
        return self.__solve_implicit(self.__implicit_multistep,t,y)[0]
    
    def Euler_AM1(self):
        self.r=1
        self.u =[1,1]
        self.v=[0,1]
        return self.__iterate(self.__multistep)