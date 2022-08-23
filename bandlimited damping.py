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
#print (v)




fre_upper=1000
fre_low=10
N=200
fre_cutoff=np.linspace(fre_low,fre_upper,N+1)
#print(fre_cutoff)
X=np.zeros((N,N))
Y=np.zeros((N))
for i in range(N):
    v_y, err_y = integrate.quad(g, fre_low, fre_upper, args = (fre_cutoff[i]))
    Y[i]=v_y
    for j in range(N):
        v_x, err_x = integrate.quad(f, fre_low, fre_upper, args = (fre_cutoff[i],fre_cutoff[j]))
        X[i,j]=v_x
#print(X)        
X_inv=np.linalg.inv(X)
belta=np.dot(X_inv,X)


cond=np.linalg.cond(X)
print(np.finfo(X.dtype).eps)


alpha=np.linalg.solve(X,Y)


def o(x):
    sum=0
    for i in range(N):
        sum=sum+alpha[i]*g(x,fre_cutoff[i])
    return sum

data_1=np.linspace(0.00001,10,5000)
data_2=np.linspace(10,1000,5000)
data_3=np.linspace(1000,100000000.0,5000)
data_x=np.concatenate([data_1,data_2,data_3])

#data_x=np.linspace(0.01,200.0,101)
eta=0.05
data_y=o(data_x)*0.05

plt.xscale('log')
plt.xlabel(u'frequency',fontproperties='SimHei')
plt.ylabel(u'damping ratio',fontproperties='SimHei')
plt.plot(data_x,data_y)  
#plt.xlim((1e-4, 100) )
plt.show()   




    


