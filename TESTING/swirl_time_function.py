#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 12:33:53 2018

@author: engels
"""


import numpy as np
import scipy
import matplotlib.pyplot as plt

def coefs( x1, x2, u1, u2):
    f1 = u1
    fp1 = 0.0
    f2 = u2
    fp2 = 0.0
    A = np.array( [ [x1**3,x1**2,x1,1.0], [3*x1**2,2*x1,1.0,0.0], [x2**3,x2**2,x2,1.0], [3*x2**2,2*x2,1.0,0.0]])
    coef = np.linalg.solve( A,  np.array([f1,fp1,f2,fp2]) )
    return coef

def velocity(u2, N):
    # end time
    t2 = 1.0
    # reversal time:
    t1 = 0.55
    # reversal duration:
    tau = 0.30

    t = np.linspace(0.0, t2, num=N, endpoint=True)

    u1 = 1.0

    u = np.zeros(t.shape) + 1000.0

    coef = coefs( 0.0, tau, 0.0, u1)
    coef2 = coefs(t1-tau, t1+tau, u1, u2)
    coef3 = coefs( t2-tau, t2, u2, 0.0)

    for it in range(len(t)):
        if ( t[it] <= t1-tau ):
            u[it]=u1

        elif t1-tau<=t[it] and t[it]<=t1+tau:
            u[it] = coef2[0]*t[it]**3 + coef2[1]*t[it]**2 + coef2[2]*t[it] + coef2[3]

        elif t[it]>t1+tau:
            u[it]=u2

        else:
            u[it] = 100.0

#     if ( t[it] <= tau ):
#            u[it] = coef[0]*t[it]**3 + coef[1]*t[it]**2 + coef[2]*t[it] + coef[3]
#        elif t[it]>tau and t[it]<t1-tau:
#            u[it]=u1
#        elif t1-tau<=t[it] and t[it]<=t1+tau:
#            u[it] = coef2[0]*t[it]**3 + coef2[1]*t[it]**2 + coef2[2]*t[it] + coef2[3]
#        elif t[it]>t1+tau and t[it]<t2-tau:
#            u[it]=u2
#        elif t[it]>t2-tau:
#            u[it]= coef3[0]*t[it]**3 + coef3[1]*t[it]**2 + coef3[2]*t[it] + coef3[3]

    return t, u

def cost(u2, N=10000):

    t, u = velocity(u2, N )

    return(np.mean( u))




u2_optimized  = scipy.optimize.bisect( cost, 0, -2.0, maxiter = 1000)

print( "%e" % (cost(u2_optimized, 10000)) )
t, u = velocity(u2_optimized, 2140 )
plt.plot(t,u, t, np.cos(2.0*np.pi*t))

uuu = scipy.integrate.cumtrapz(u, x=t, initial=0.0)
plt.figure()
plt.plot(t,uuu)
print(uuu[-1])

u2 = u2_optimized
N = 1000
# end time
t2 = 1.0
# reversal time:
t1 = 0.55
# reversal duration:
tau = 0.30

t = np.linspace(0.0, t2, num=N, endpoint=True)

u1 = 1.0

u = np.zeros(t.shape)

coef = coefs( 0.0, tau, 0.0, u1)
coef2 = coefs(t1-tau, t1+tau, u1, u2)
coef3 = coefs( t2-tau, t2, u2, 0.0)

for it in range(len(t)):
    if ( t[it] <= tau ):
        u[it] = coef[0]*t[it]**3 + coef[1]*t[it]**2 + coef[2]*t[it] + coef[3]
    elif t[it]>tau and t[it]<t1-tau:
        u[it]=u1
    elif t1-tau<=t[it] and t[it]<=t1+tau:
        u[it] = coef2[0]*t[it]**3 + coef2[1]*t[it]**2 + coef2[2]*t[it] + coef2[3]
    elif t[it]>t1+tau and t[it]<t2-tau:
        u[it]=u2
    elif t[it]>t2-tau:
        u[it]= coef3[0]*t[it]**3 + coef3[1]*t[it]**2 + coef3[2]*t[it] + coef3[3]