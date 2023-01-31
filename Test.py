# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 20:36:49 2022

@author: ppf
"""

import numpy as np
import matplotlib.pyplot as plt
from Jet import Jet
from Method import Explicit, Implicit

def f(t,y):
    return -10*y
    # return Jet.sin(y**(3/2))
    # return Jet.exp(-Jet.sin(y)**2)

def __Explicit_Euler(h,f,t_0,y_0):
    return y_0+h*f(t_0,y_0)
   
def __Implicit_Euler(h,f,t,y_0,y):
    return y_0+h*f(t+h,y)-y

def __Implicit_Trapezoid(h,f,t,y_0,y):
    return y_0+(h/2)*(f(t,y_0)+f(t+h,y))-y

def __Solve_Implicit(method,h,f,t,y_0):
            
    #Solve first the zero-order jet case
    y0_order0 = Jet([y_0.getJet()[0]])
    y=y0_order0
            
    Newt_it=5
    for i in range(Newt_it):
        method_y=method(h,f,t,y0_order0,y)
        dmethod_y=Jet([(method(h,f,t,y0_order0,Jet([y.getJet()[0],1]))-method_y).getJet()[1]])  #Notice that F(y+s)=F(y)+F'(y_1)s
        y=y-method_y/dmethod_y
                    
    #Solve the general N-order jet case
    for n in range(1,y_0.getOrder()):
        
        #Notice F(a+y_0^[N]s^N, b)=F(a,b)+D_{y_0}F(a,b)y_0^[N]s^N
        #Notice F(a, b+s^N)=F(a,b)+D_{y}F(a,b)s^N
        #Notice y^[N]=-D_{y_0}F(a,b)y_0^[N]/D_{y}F(a,b)
                
        y0_ncoeff=y_0.getJet()[n]
        D_y0F_y0=(method(h,f,t,Jet([*y0_order0.getJet(),y0_ncoeff]),y)-method(h,f,t,y0_order0,y)).getJet()[n]
        D_yF=(method(h,f,t,y0_order0,Jet([*y.getJet(),1]))-method(h,f,t,y0_order0,y)).getJet()[n]
        
        y_n=-D_y0F_y0/D_yF
        y=Jet([*y.getJet(),y_n])
        
        y0_order0=Jet(y_0.getJet()[0:n])
        
    return y
    

# steps, t, sols =  Explicit(h=0.1, t_0=0, t_F=2, jet_0=Jet([2.0,1.0]), f=f).Euler()
steps, t, sols =  Implicit(h=0.1, t_0=0, t_F=2, jet_0=Jet([2.0,1.0]), f=f).Euler()



for sol in sols:
      plt.plot(t,sol)
                 
print('steps: '+ str(steps))


# t=[]; jets=[]
# h=0.1
# t_0=0
# t_F=2
# # jet_0=2.0
# jet_0=Jet([2.0,1.0])
# # jet_0=Jet([2.0,1.0])

# t.append(t_0)
# jets.append(jet_0)

# t.append(t[-1]+h)
# while t[-1]<=t_F:
#     jets.append(__Solve_Implicit(__Implicit_Euler, h, f, t[-1], jets[-1]))
#     # jets.append(__Solve_Implicit(__Implicit_Trapezoid, h, f, t[-1], jets[-1]))
#     # jets.append(__Explicit_Euler(h,f,t[-1],jets[-1]))
#     t.append(t[-1]+h)
# t.pop()
        
# # steps, t, sols = len(t),t, [jets]
# steps, t, sols = len(t),t, np.array([d.getJet() for d in jets]).transpose()
# for sol in sols:
#       plt.plot(t,sol)
      
      
X=np.linspace(0,2,100)
plt.plot(X,[2*np.exp(-10*x) for x in X], 'm-')
plt.plot(X,[1*np.exp(-10*x) for x in X], 'r-')

plt.legend(["sol.","var. sol.", "exact", "var. exact"])
