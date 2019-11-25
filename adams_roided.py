# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 18:38:02 2019

@author: Adam
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from pylab import pcolor
# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

nr = 20
ntheta = 50
nr = ntheta
dt = .01
time = 10


###Propeties taken from https://www.researchgate.net/publication/260491141_Modeling_of_one-dimensional_thermal_response_of_silica-phenolic_composites_with_volume_ablation
sigma = 5.67E-8
R = 8.314
poly_frac = .51
fibre_frac = .38
char_frac = 1/5 #mass fraction of polymer that reacts and turns into char 
#phase 1 mat properties
k = .8 #w/mk
rho = 1100 #kg/m^3
cp = 1.09E3 #j/kgk
def cp(T):
    return (1.09 + 1.09E-3*T)*1E3
def k(T):
    return .803 + 2.76E-4 * T

#phase 2 mat properties
k2 = .9
rho2 = 2000
cp2 = .88E3
def cp2(T):
    return (.88 + 1.02E-3*T)*1E3
def k2(T):
    return .954+8.42E-4*T- 4.07E-6 * T**2 + 5.31E-9 * T**3
emmisivity = 0

rhof = 2200

def alpha(T):
    return 100 * k(T)/(cp(T)*rho)
def alpha2(T):
    return k2(T)/(cp2(T)*rho2)

#Chemical reaction properties
Ea = 383E3 / 4.5#j/mol
M = 0.018 #kg/mol
Ac = 3.63E6 # kg/m
dH = 418E3 * 2.5#J/kg

hg = 10.8377E3# 2.6E3  #w/m^2k
Taw = 3300 #k   





#def runSim(length,hg,Taw,ri,disp = True):

length,hg,Taw,ri,disp = .015,2593,2880,.0165,True
dr = length/nr
dtheta = 2*np.pi/ntheta
rarr = np.linspace(ri,ri+length,nr)
thetaarr = np.linspace(0,2*np.pi,ntheta)

T = np.zeros((2,ntheta,nr))



A1,A2,A3 = np.zeros((ntheta,nr)),np.zeros((ntheta,nr)),np.zeros((ntheta,nr))

#DEFINE BOUNDARY CONDITIONS
A1[0][0] = -1   
A1[0][1] = 1 
A1[-1][-1] = -1
A1[-1][-2] = 1

A2[0][0] = -1/2   
A2[0][1] = 1/2 
A2[-1][-1] = 1/2
A2[-1][-2] = -1/2
A2[0] = A2[0] * (1/rarr[0])
A2[-1] = A2[-1] * (1/rarr[-1])

A3[0][0] = -1   
A3[0][1] = 1 
A3[-1][-1] = -1
A3[-1][-2] = 1
A3[0] = A3[0] * (1/rarr[0])**2
A3[-1] = A3[-1] * (1/rarr[-1])**2

for i in range (1,ntheta-1):
    A1[i][i-1] = 1 
    A1[i][i+1] = 1
    A1[i][i] = -2
    
    A2[i][i-1] = -1/2 
    A2[i][i+1] = 1/2
    A2[i][i] = 0
    A2[i] = A2[i] * (1/rarr[i])
    
    A3[i][i-1] = 1 
    A3[i][i+1] = 1
    A3[i][i] = -2
    A3[i] = A3[i] * (1/rarr[i])**2


A1 = A1 / dr**2
A2 = A2 / (dr)
A3 = A3 / dtheta**2


T[0,:,:] = 285
T[0,0,:] = 280
dT = np.zeros((ntheta,nr))
for i in range(int(time/dt)):
    for count,Tr in enumerate(T[0]):
        dT[count,:] = alpha(Tr) * ((np.matmul(A1,Tr) + np.matmul(A2,Tr)))
    for count,Ttheta in enumerate(T[0].T):
        dT[:,count] += alpha(Ttheta) * (np.matmul(A3,Ttheta))
    T[1] = T[0] + dT*dt
    ttemp = T
    T[0] = ttemp[1]
    printProgressBar(i,int(time/dt))
    
if(disp):
    #plot the final value
    fig = plt.figure()
    ax = fig.add_subplot(111)
    R,Theta = np.meshgrid(rarr,thetaarr)
    X, Y = R*np.cos(Theta), R*np.sin(Theta)
    im = ax.pcolormesh(X, Y, T[0].T,cmap='plasma')
    fig.colorbar(im)
            
    
