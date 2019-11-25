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

mtoin = 39.3701
RPAname =  "350psi2.8MR11ContractThermo.txt"
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

chamberEndIndex = radius != max(radius)
loc,radius,h,tAw = loc[chamberEndIndex],radius[chamberEndIndex],h[chamberEndIndex],tAw[chamberEndIndex]
maxT = 313
def fMin(l,i):
    return abs(maxT - sim.runSim(l,h[i],tAw[i],radius[i])[0][-1])

results = []
for i in range(len(h)):
    results.append(opt.fminbound(fMin,0.002,0.03,args=(i,),disp=3))
    print(results,h[i],tAw[i],radius[i],i)
    
#results = np.array([0.008965688169238822, 0.008966559177213455, 0.00896751008122395, 0.008970438152918222, 0.008975959987590312, 0.008979287876320395, 0.008979853811019322, 0.008983469699604086, 0.008993852269194295, 0.008997867112812353, 0.009002243702006776, 0.009006706495569838, 0.009011230992000052, 0.009019178031146739, 0.009023770529651426, 0.00902838768976957, 0.009033774494744287, 0.009036989020721697, 0.009040462750671809, 0.009042375764882695, 0.00904723338847853, 0.009045088538341023, 0.009044255907142588, 0.009045034358887481, 0.009039452544342552, 0.009035762384702613, 0.009029217322208426, 0.009022826045050466, 0.009011522143028685, 0.008997706501592704, 0.008978712151546614, 0.008875388202501892, 0.008875388202501892, 0.008875388202501892, 0.008875388202501892, 0.008981134368388639, 0.009007518096835877])
#find only throat
#ind = np.where(radius==min(radius))
#ind = int(ind[0])
#throatThick = opt.fminbound(fMin,0.002,0.03,args=(ind,))
#results = np.array([throatThick])
#print(str(results*mtoin) + " thickness in inches at throat")

results = np.array(results)
loc = loc*100/2.54
radius = radius*100/2.54
results = results*100/2.54

ln2, = ax2.plot(loc,radius,'b',label="Engine R")
ln, = ax.plot(loc[:len(results)],results,'r',label="Layer T")

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


