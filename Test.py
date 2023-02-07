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
    
def exactsol(x): return 2*np.exp(-10*x)
def exactvarsol(x): return 1*np.exp(-10*x)

h=0.0125; t_0=0; t_F=2; jet_0=Jet([2.0,1.0]); Newt_iter=10

EXPLICIT = Explicit(h=h, t_0=t_0, t_F=t_F, jet_0=jet_0, f=f)
IMPLICIT = Implicit(h=h, t_0=t_0, t_F=t_F, jet_0=jet_0, f=f, Newt_iter=Newt_iter)



# steps, t, sols =  EXPLICIT.Euler()
steps, t, sols, local_errors, global_errors =  IMPLICIT.Euler(exact_sol=[exactsol, exactvarsol])
# steps, t, sols =  IMPLICIT.Midpoint()
# steps, t, sols =  IMPLICIT.Crank_Nicolson()
# steps, t, sols =  IMPLICIT.QZ_DIRK()
# steps, t, sols =  IMPLICIT.KS_DIRK()
# steps, t, sols =  IMPLICIT.RADAU_IA_3o()
# steps, t, sols =  IMPLICIT.RADAU_IIA_3o()

# for sol in sols:
#     plt.plot(t,sol)
    
for local_error in local_errors:
    plt.plot(t[1:],local_error)
for global_error in global_errors:
    plt.plot(t[1:],global_error)
    
plt.yscale('log')

plt.legend(["sol. local","diff. sol. local","sol. global","diff. sol. global"])

# plt.plot(t,sols[0])
# plt.plot(t,sols[1])
# print(sols[0][2]-exactsol(t[2]))
# print(sols[0][4]-exactsol(t[4]))

                 
print('steps:',steps)



X=np.linspace(t_0,t_F,100)
# plt.plot(X,exactsol(X), 'm-')
# plt.plot(X,exactvarsol(X), 'r-')

# plt.legend(["sol.","var. sol.","diff", "exact", "var. exact"])


print("sol. max error:",max([abs(sols[0][i]-exactsol(t[i])) for i in range(steps)]))
print("var. sol. max error:",max([abs(sols[1][i]-exactvarsol(t[i])) for i in range(steps)]))


















# def __Explicit_Euler(h,f,t_0,y_0):
#     return y_0+h*f(t_0,y_0)
   
# def __Implicit_Euler(h,f,t,y_0,y):
#     return y_0+h*f(t+h,y)-y

# def __Implicit_Trapezoid(h,f,t,y_0,y):
#     return y_0+(h/2)*(f(t,y_0)+f(t+h,y))-y

# def __Solve_Implicit(method,h,f,t,y_0):
            
#     #Solve first the zero-order jet case
#     y0_order0 = Jet([y_0.getJet()[0]])
#     y=y0_order0
            
#     Newt_it=5
#     for i in range(Newt_it):
#         method_y=method(h,f,t,y0_order0,y)
#         dmethod_y=Jet([(method(h,f,t,y0_order0,Jet([y.getJet()[0],1]))-method_y).getJet()[1]])  #Notice that F(y+s)=F(y)+F'(y_1)s
#         y=y-method_y/dmethod_y
                    
#     #Solve the general N-order jet case
#     for n in range(1,y_0.getOrder()):
        
#         #Notice F(a+y_0^[N]s^N, b)=F(a,b)+D_{y_0}F(a,b)y_0^[N]s^N
#         #Notice F(a, b+s^N)=F(a,b)+D_{y}F(a,b)s^N
#         #Notice y^[N]=-D_{y_0}F(a,b)y_0^[N]/D_{y}F(a,b)
                
#         y0_ncoeff=y_0.getJet()[n]
#         D_y0F_y0=(method(h,f,t,Jet([*y0_order0.getJet(),y0_ncoeff]),y)-method(h,f,t,y0_order0,y)).getJet()[n]
#         D_yF=(method(h,f,t,y0_order0,Jet([*y.getJet(),1]))-method(h,f,t,y0_order0,y)).getJet()[n]
        
#         y_n=-D_y0F_y0/D_yF
#         y=Jet([*y.getJet(),y_n])
        
#         y0_order0=Jet(y_0.getJet()[0:n])
        
#     return y

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
      
      

