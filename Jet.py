# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 20:36:49 2022

@author: ppf
"""

import numpy as np

class Jet:
    def __init__(self,jet):
        self.jet=jet
        self.order=len(jet)
        
    def __repr__(self):
        s = ''
        for i in range(self.order):
            if i==0: s+=str(self.jet[i])
            elif i==1: s+='+'+str(self.jet[i])+'x' if self.jet[i]>=0 else str(self.jet[i])+'x'
            else: s+= '+'+str(self.jet[i])+'x^'+str(i) if self.jet[i]>=0 else str(self.jet[i])+'x^'+str(i)
        return str(s)
    
    def __add__(self, other):
        if isinstance(other, Jet) and self.order==other.getOrder():
            result = np.zeros(self.order)
            otherJet=other.getJet()
            for n in range(self.order):
                result[n]=self.jet[n]+otherJet[n]
            return Jet(result)
        if isinstance(other, float) or isinstance(other, int):
             return Jet([other+d for d in self.jet])
         
    def __radd__(self, val):
        if isinstance(val, float) or isinstance(val, int):
            return self+val

    def __sub__(self, other):
        if isinstance(other, Jet) and self.order==other.getOrder():
            result = np.zeros(self.order)
            otherJet=other.getJet()
            for n in range(self.order):
                result[n]=self.jet[n]-otherJet[n]
            return Jet(result)
        if isinstance(other, float) or isinstance(other, int):
             return Jet([other-d for d in self.jet])
         
    def __rsub__(self, val):
        if isinstance(val, float) or isinstance(val, int):
            return self-val
        
    def __mul__(self, other):
        if isinstance(other, Jet) and self.order==other.getOrder():
            result = np.zeros(self.order)
            otherJet=other.getJet()
            for n in range(self.order):
                for j in range(n+1):
                    result[n]+=self.jet[n-j]*otherJet[j]
            return Jet(result)
        if isinstance(other, float) or isinstance(other, int):
             return Jet([other*d for d in self.jet])
         
    def __rmul__(self, val):
        if isinstance(val, float) or isinstance(val, int):
            return self*val
        
    def __truediv__(self, other):
        if isinstance(other, Jet) and self.order==other.getOrder():
            result = np.zeros(self.order)
            otherJet=other.getJet()
            result[0]=self.jet[0]/other.jet[0]
            for n in range(1,self.order):
                result[n]+=self.jet[n]
                for j in range(1,n+1):
                    result[n]-=result[n-j]*otherJet[j]
                result[n]/=otherJet[0]
            return Jet(result)
        if isinstance(other, float) or isinstance(other, int):
             return Jet([d/other for d in self.jet])
         
    def __rtruediv__(self, val):
        if isinstance(val, float) or isinstance(val, int):
            return val*self.__pow__(-1)
   
    def __pow__(self, val):
        if isinstance(val, float) or isinstance(val, int):
            result = np.zeros(self.order)
            result[0]=self.jet[0]**val
            for n in range(1,self.order):
                for j in range(n):
                    result[n]+=(n*val-j*(val+1))*self.jet[n-j]*result[j]
                result[n]/=(n*self.jet[0])
            return Jet(result)
        
    def __neg__(self):
        return Jet([-d for d in self.jet])
    
    def getJet(self): 
        return self.jet
    
    def getOrder(self): 
        return self.order

    @classmethod
    def exp(cls, other):
        if isinstance(other, cls):
            otherJet=other.getJet()
            result = np.zeros(other.order)
            result[0]=np.exp(otherJet[0])
            for n in range(1,other.order):
                for j in range(n):
                    result[n]+=(n-j)*otherJet[n-j]*result[j]
                result[n]/=n
            return cls(result)

    @classmethod
    def ln(cls, other):
        if isinstance(other, cls):
            otherJet=other.getJet()
            result = np.zeros(other.order)
            result[0]=np.log(otherJet[0])
            for n in range(1,other.order):
                result[n]+=otherJet[n]
                for j in range(1,n):
                    result[n]-=(n-j)/n*otherJet[j]*result[n-j]
                result[n]/=otherJet[0]
            return cls(result)
    
    __sin=[]; __cos=[]
    def _computesincos(self):
        self.__sin = np.zeros(self.order)
        self.__cos = np.zeros(self.order)
        self.__sin[0]=np.sin(self.jet[0])
        self.__cos[0]=np.cos(self.jet[0])
        for n in range(1,self.order):
            for j in range(1,n+1):
                self.__sin[n]+=j*self.__cos[n-j]*self.jet[j]
                self.__cos[n]-=j*self.__sin[n-j]*self.jet[j]
            self.__sin[n]/=n
            self.__cos[n]/=n

    @classmethod
    def sin(cls, other):
        if len(other.__sin)==0: other._computesincos()
        return cls(other.__sin)
     
    @classmethod
    def cos(cls, other):
        if len(other.__cos)==0: other._computesincos()
        return cls(other.__cos)
    
    @classmethod
    def dot(cls, others, values):
        if len(others)==1 and len(values)==1:
            return others[0]*values[0]
        if len(values)==len(others):
            result=values[0]*others[0]
            for i in range(1,len(values)):
                result+=values[i]*others[i]
            return result
        
    @classmethod
    def naiveDerivative(cls, other, f, y_0, *, F, t, eps=0.001):
        if isinstance(other,cls):
            return (F(f,t,other+eps,y_0)-F(f,t,other,y_0))/eps
        