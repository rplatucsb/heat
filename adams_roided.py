# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 18:38:02 2019

@author: Adam
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import sparse
import time
# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()
np.set_printoptions(precision=5)
#np.set_printoptions(suppress=True)

#sim mesh options
nr = 200
ntheta = 2
dt = .001
nettime = 13
height = 1

###Propeties taken from https://www.researchgate.net/publication/260491141_Modeling_of_one-dimensional_thermal_response_of_silica-phenolic_composites_with_volume_ablation
sigma = 5.67E-8
gascnst = 8.314
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
    return k(T)/(cp(T)*rho)
def alpha2(T):
    return k2(T)/(cp2(T)*rho2)

#Chemical reaction properties
Ea = 383E3 / 4.5#j/mol
M = 0.018 #kg/mol
Ac = 3.63E6 # kg/m
dH = 418E3 * 2.5#J/kg

hg = 10.8377E3# 2.6E3  #w/m^2k
Taw = 3300 #k   



def runSim(length,hg,Taw,ri,disp = True):
    timestart = time.time()
    dt = (length/nr)**2/(2*alpha(280)) * .05
    #length,hg,Taw,ri,disp = .015,2593,2880,.0165,True
    dr = length/nr
    dtheta = 2*np.pi/ntheta
    rarr = np.linspace(ri,ri+length,nr)
    thetaarr = np.linspace(0,2*np.pi,ntheta)
    vEl = np.array(dtheta * np.pi * height * (2 * rarr * dr + dr**2))
    
    #FORMAT is T(theta,radius)
    #IE iso theta is rows, iso R is columns
    T = np.zeros((2,ntheta,nr)) 
    V = np.ones((2,ntheta,nr))#volume fraction
    M = np.ones((2,ntheta,nr))#mass fraction
    mEl = np.zeros((2,ntheta,nr))#mass of individual elements
    
    #a1,a2 are second and first T derivative wrt R
    #a3 is second order T derivative with theta
    A1,A2,A3 = np.zeros((nr,nr)),np.zeros((nr,nr)),np.zeros((ntheta,ntheta))
    
    #DEFINE BOUNDARY CONDITIONS
    A1[0][0] = -2
    A1[0][1] = 1
    A1[-1][-1] = -1
    A1[-1][-2] = 1
    
    A2[0][1] = 1/2 * (1/rarr[0])
    A2[-1] = 0
    
    A3[0][0] = -2   
    A3[0][1] = 1 
    A3[0][-1] = 1
    A3[-1][-1] = -2
    A3[-1][-2] = 1
    A3[-1][0] = 1
    #central finite differences in r and theta
    for i in range (1,nr-1):
        A1[i][i-1] = 1 
        A1[i][i+1] = 1
        A1[i][i] = -2
        
        A2[i][i-1] = -1/2 
        A2[i][i+1] = 1/2
        A2[i][i] = 0
        A2[i] = A2[i] * (1/rarr[i])
        
    for i in range(1,ntheta-1):
        A3[i][i-1] = 1 
        A3[i][i+1] = 1
        A3[i][i] = -2
        A3[i] = A3[i]
    
    
    A1 = A1 / dr**2
    A2 = A2 / (dr)
    A3 = A3 / dtheta**2
    
    A12 = A1 + A2

    #A1 = sp.sparse.csr_matrix(A1)
    #A2 = sp.sparse.csr_matrix(A2)
    
    
    T[0,:,:] = 300
    mEl[0,0,:] = poly_frac*vEl*rho+fibre_frac*vEl*rhof
    dT = np.zeros((ntheta,nr))
    dV = np.zeros((ntheta,nr))

    for i in range(int(nettime/dt)):
        
        """ Original code, capable of heat flow in r and theta directions
        for count,Tr in enumerate(T[0]): #in a single iteration, we have a matrix of temps at a certain theta
            dT[count,:] = alpha(Tr) * ((np.matmul(A1,Tr) + np.matmul(A2,Tr)))
        for count,Ttheta in enumerate(T[0].T):#in a single iteration, we have a matrix of temps at a certain radius
            dT[:,count] += alpha(Ttheta) * ((np.matmul(A3,Ttheta))*(1/rarr[count]**2))
        """
        """ This code runs much faster, but can only do axisymmetric"""
        ''' Define mass averaged properties at nodes '''
        cpeff = (M[0,0,:]*cp(T[0,0,:])+(1-M[0,0,:])*cp2(T[0,0,:]))
        aeff = (M[0,0,:]*alpha(T[0,0,:])+(1-M[0,0,:])*alpha2(T[0,0,:]))
        k0 = (M[0,0,0]*k(T[0,0,0])+(1-M[0,0,0])*k(T[0,0,0]))
        '''define heat influx boundary condition '''
        imaginaryTemp = T[0,0,1] - (2*dr*hg/k0)*(T[0,0,0]-Taw)
        '''calculate conduction '''
        dTuniformconduction = dt*aeff * (np.matmul(A12,T[0,0,:])) 
        dTuniformconduction[0] = dTuniformconduction[0] + dt*aeff[0] * (imaginaryTemp/dr**2 - imaginaryTemp/(2*dr*rarr[0]))
        
        '''calculate ablation '''
        reactionRate = Ac * (poly_frac) * (V[0,0,:]) * np.exp(-Ea/(gascnst*T[0,0,:])) 
        dV = -reactionRate * dt / rho
        dV[V[0,0,:]-dV<0] = 0
        dE = dV * vEl*poly_frac * rho * dH 
        dTuniformablation = dE / (mEl[0,0,:] * cpeff)
        
        #sparse matrix implementation
    #    dtuniform = alpha(T[0,0,:]) * (A1 @ T[0,0,:] + A2 @ T[0,0,:])
    #    dT[:] = dtuniform
        dT[:] = dTuniformconduction + dTuniformablation
        dV[:] = dV

        T[1] = T[0] + dT
        V[1] = V[0] + dV
        M[:] = (V[1] * rho ) / (V[1]*rho+(1-V[1])*rho2)
        mEl[:] = poly_frac  * vEl * (V[0]*rho + char_frac*(1-V[0])*rho2) + rhof * fibre_frac * vEl
        
        Ttemp = T
        mElTemp = mEl
        Vtemp = V
        Mtemp = M
        T[0] = Ttemp[1]
        V[0] = Vtemp[1]
        M[0] = Mtemp[1]
        mEl[0] = mElTemp[1]
        
        printProgressBar(i,int(nettime/dt))
    if(disp):
        #plot the final T value in a donut
        fig = plt.figure()
        ax = fig.add_subplot(111)
        R,Theta = np.meshgrid(rarr,thetaarr)
        X, Y = R*np.cos(Theta), R*np.sin(Theta)
        im = ax.pcolormesh(X, Y, T[0],cmap='plasma')
        fig.colorbar(im)
        
        #plot the mass fraction and Temperature as a function of radius
        plt.figure(3)
        fig,ax = plt.subplots()
        ax2 = ax.twinx()
        ln, = ax.plot(rarr,T[0,0,:],'r')
        ln2, = ax2.plot(rarr,V[0,0,:])
        ax.set_ylabel("Temperature(K)")
        ax2.set_ylabel("Volume Fraction of original material")
        plt.xlabel("Radius (m)")
        print("Runtime: " + str(time.time() - timestart))
    return(T[0],V[0])
    
#tdat = []
#for i in np.linspace(600,1000,4):
#    nr = int(i)
#    T = runSim(.015,2593,2880,.0165,True)
#    tdat.append(np.max(T)-np.min(T))
