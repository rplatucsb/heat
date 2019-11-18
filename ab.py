# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 12:50:25 2019

@author: Adam
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
RPAname = "rpa300psithermo.txt"
f = open(RPAname,'r')
lines = f.readlines()
lines = lines[10:-2]
arr = []

for line in lines:
    line = line[:-1]
    vals = line.split('\t')
    row =[]
    for el in vals:
        try:
            row.append(float(el))
        except:
            continue
    arr.append(row)

arr = np.array(arr)

Hname = "X6002B.txt"
f = open(Hname,'r')
lines = f.readlines()
rarr = [[]]
for line in lines:
    line = line.replace("\n","")
    line = line.replace("x ","")
    line = line.replace("-","E-")
    vals = line.split(" ")
    rarr.append([float(i) for i in vals])
rarr.remove([])
matDat = np.array(rarr)
print(matDat)
qDat = matDat[:,2]
hDat = matDat[:,0]


 
#hDat = np.array([4300.,4415.,5700.,6670.,6990.,7650.,8270.])
#qDat = np.array([905.,865.,1170.,1550.,1240.,1326.,1130.])
qDat *= 11.356539E3
hDat *= 1055.066*2.2
tArr = np.zeros((np.size(qDat),2))
tArr[:,1] = hDat
tArr[:,0] = qDat
print(tArr)
tArr = tArr[np.argsort(tArr[:,0])]

print(tArr)
hOfQ = np.polyfit(tArr[:,0],tArr[:,1],1)
hOfQFunc = np.poly1d(hOfQ)

loc = arr[:,0]*1E-3
radius = arr[:,1]*1E-3
qflux = arr[:,3]*1000
h = arr[:,2] * 1000
qs = np.arange(min(qflux),max(qflux),100)

plt.figure(1)
plt.xlabel("q (w/m^2)")
plt.ylabel("H (J/kg)")
plt.plot(tArr[:,0],tArr[:,1],label = "Data")
plt.plot(qs,hOfQFunc(qs), label = "Extrapolation over engine qs")
plt.legend()

thicknesses = []
mTot = 0
QQ = 0
for lInd,l in enumerate(loc):
    h = loc[lInd]-loc[lInd-1]
    if h < 0: h = 0.013410000000000005
    r= radius[lInd]
    sA = 2*np.pi*r*h
    Q = qflux[lInd]*13*sA
    QQ += Q
    H = hOfQFunc(qflux[lInd])
    rho = 1778
    m = Q/H
    mTot += m
    V = m / rho
    t = np.sqrt(r**2+V/(np.pi*h))-r
    thicknesses.append(t)

print(mTot, "Estimated ablative mass (kg)")
thicknesses = np.array(thicknesses)

plt.figure(2)
plt.xlabel("Length along engine(cm)")
plt.ylabel("Ablative thickness(mm)")
plt.plot(loc*1E2,thicknesses*1E3)
plt.show()
f.close()
    