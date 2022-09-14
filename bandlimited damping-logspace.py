# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 14:53:41 2021

@author: 2182
"""


import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt


def f(x,x_i,x_j):
    g_i=2*(x/x_i)/(1+(x/x_i)**2)
    g_j=2*(x/x_j)/(1+(x/x_j)**2)
    return g_i*g_j

def g(x,x_i):
    g_i=2*(x/x_i)/(1+(x/x_i)**2)
    return g_i

def g_integrate(x,x_i):
    g_i=2*x*x_i/(x**2+x_i**2)
    return g_i
#print (v)



# fre_low=10
# fre_upper=100
fre_low=0.1
fre_upper=1

iterations=20
def uniform_damp(fre_data,fre_low,fre_upper,iterations=100):
    for i in range(iterations):
        print(i)
        N=12+i
        fre_cutoff=np.logspace(np.log10(fre_low),np.log10(fre_upper), num=N, base=10)  
        # print(fre_cutoff)
        X=np.zeros((N,N))
        Y=np.zeros((N))
        for i in range(N):
            v_y, err_y = integrate.quad(g, fre_low, fre_upper, args = (fre_cutoff[i]))
            # v_y=g_integrate(fre_upper,fre_cutoff[i])-g_integrate(fre_low,fre_cutoff[i])
            print(v_y)
            Y[i]=v_y
            for j in range(N):
                v_x, err_x = integrate.quad(f, fre_low, fre_upper, args = (fre_cutoff[i],fre_cutoff[j]))
                X[i,j]=v_x
        #print(X)        
        X_inv=np.linalg.inv(X)
        belta=np.dot(X_inv,X)
        alpha=np.linalg.solve(X,Y)

        damping_factor=0
        err=0
        for j in range(N):
            damping_factor=damping_factor+alpha[j]*g(fre_data,fre_cutoff[j])
            err=err+alpha[j]*g(fre_cutoff,fre_cutoff[j])
        if (abs(err-1.0)<=0.01).all():
            print(N)
            return damping_factor  ,Y ,alpha,fre_cutoff     
        


data_1=np.logspace(np.log10(fre_low/1e5),np.log10(fre_low),num=100, base=10)
data_2=np.logspace(np.log10(fre_low),np.log10(fre_upper),num=100, base=10)
data_3=np.logspace(np.log10(fre_upper),np.log10(fre_upper*1e5),num=100, base=10)
data_x=np.concatenate([data_1,data_2,data_3])

#data_x=np.linspace(0.01,200.0,101)
eta=0.05
data_y,y_temp,alpha,fre_cutoff=uniform_damp(data_x,fre_low,fre_upper,iterations)

plt.xscale('log')
plt.xlabel(u'frequency',fontproperties='SimHei')
plt.ylabel(u'damping ratio',fontproperties='SimHei')
plt.plot(data_x,data_y,label=str(fre_low)+'-'+str(fre_upper)) 
plt.legend(loc=0) 
#plt.xlim((1e-4, 100) )
plt.show()   




    


