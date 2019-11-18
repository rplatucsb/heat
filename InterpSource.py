# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 21:32:31 2019

@author: Adam
thermoclass p is in bar
"""

import numpy as np
import thermoClass as T
import matplotlib.pyplot as plt
methane = T.thermo()
#file = open("hfg.txt","w")
a = np.arange(110,150,.5)
b = []
k = []
c = []
#print(methane.gasTemp(30))
#file.write("%Temps\n") 
#for i in a:
#    bb=methane.eqState(1,i)
#    b.append(bb)
#    gg = methane.latentVap(bb,i)
#    k.append(gg)
#    c.append(methane.gasTemp(bb))
#    #file.write(str(i)+" ")
#file.write("\n%CPs\n")
#for i in k:
#    file.write(str(i) + " ")
#file.close
print(methane.latentVap())
#plt.subplot(2,1,1)
#z = np.array([8163,7911,7071,6619])*1000*62.332
#plt.plot([111.88,119.96,140.17,150],z)
#plt.subplot(2,1,2)
#plt.plot(a,k)

