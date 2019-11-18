# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 13:14:57 2019

@author: Adam
"""

import numpy as np
import nolans_feeble as sim
import matplotlib.pyplot as plt
import scipy.optimize as opt

plt.close("all")
fig,ax = plt.subplots()
ax2 = ax.twinx()


RPAname =  "500psithermo.txt"
#RPAname = "rpa300psithermo.txt"
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

loc = arr[:,0]*1E-3
radius = arr[:,1]*1E-3
h = arr[:,2] * 1000
tAw = arr[:,6]

maxT = 313
def fMin(l,i):
    return abs(maxT - sim.runSim(l,h[i],tAw[i],radius[i])[0][-1])

#results = []
#for i in range(len(h)):
#    i=-1
#    results.append(opt.fminbound(fMin,0.002,0.02,args=(1,)))
#    print(results,h[i],tAw[i],radius[i])
#    
#results = np.array([0.00911312810007722, 0.009113663223527996, 0.009114142611088857, 0.009114194833689774, 0.009114194833689774, 0.009114194833689774, 0.009114194833689774, 0.009114194833689774, 0.009114194833689774, 0.009114194833689774, 0.009114194833689774, 0.009114194833689774, 0.009114194833689774, 0.009114194833689774, 0.009114194833689774, 0.009114194833689774, 0.009114194833689774, 0.009114131548339193, 0.009120386141606, 0.009138508391178794, 0.009169258801571852, 0.00921159404780107, 0.009250289011530324, 0.009268646543756231, 0.009254865791769214, 0.009241499549163945, 0.009109280485164278, 0.009008207845058317, 0.008904521240283656])

loc = loc*100/2.54
radius = radius*100/2.54
results = results*100/2.54

ln2, = ax2.plot(loc,radius,'b',label="Engine R")
ln, = ax.plot(loc,results,'r',label="Layer T")

m = 0
rho = 1000
for i in range(len(results)):
    he = loc[i]-loc[i-1]
    if he < 0: he = 0.013410000000000005
    m += np.pi * he * ((radius[i]+results[i])**2-(radius[i])**2) * rho
    print(m)
    
print(m, "kg ablative")
#ax.set_ylim([results,.01])
ax.set_ylabel("Ablative layer thickness (in)")
ax2.set_ylabel("Engine radius (in)")
ax.set_xlabel("Length along Engine (in)")
plt.legend((ln,ln2),("Layer T","Engine R"))

plt.show()


