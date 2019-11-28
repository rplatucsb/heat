# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 12:27:03 2019

@author: Adam
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import pickle
import time

plt.close("all")
plt.rcParams["animation.ffmpeg_path"] = r"C:\Users\Adam\Desktop\RPL\Heet\ffmpeg"
fig,ax = plt.subplots()
ax2 = ax.twinx()
mspf = 20
animtime = 13

heatRate = .1

n = 500
dt = .01
length = .01 #depth of ab layer
dx = length/n
btime = 13 #* 10/(heatRate)
sA = np.array([.0001]) #m^2
ri = .03659 / 2
dy = .02
vEl = dx * sA

#def Q(t):
#    if t<60:
#        return 1E5
#    elif t<80:
#        return 3E5
#    else:
#        return 2E5



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

def a(T):
    return k(T)/(cp(T)*rho)
def a2(T):
    return k2(T)/(cp2(T)*rho2)

#Chemical reaction properties
Ea = 383E3 / 4.5#j/mol
M = 0.018 #kg/mol
Ac = 3.63E6 # kg/m
dH = 418E3 * 2.5#J/kg

hg = 10.8377E3# 2.6E3  #w/m^2k
Taw = 3300 #k   

def runSim(length,hg,Taw,ri,disp = False): 
    dx = length/n
    dt = dx**2/(2*a(285)) * .9
    rArr = np.linspace(ri,ri+length,n)
    vEl = np.pi*dy*((rArr+dx)**2-rArr**2)
    sA = 2 * np.pi * dy * (rArr)
#    
    t = np.zeros(n)
    T = np.zeros((2,n))
    V = np.zeros((2,n))
    M = np.zeros((2,n))
    mEl = np.zeros((2,n))
    qIn = np.zeros((2))
    
    T[0] = np.linspace(150,50,n)
    T[0].fill(300)
    V[0] = np.ones(n) #Resin is only 51% of the material initially, this fraction represents fraction of that 51%. 
    M[0] = np.ones(n) #resin mass fraction, again, see above. This will be used to evaluate virgin v char material properties
    mEl[0] = poly_frac*vEl*rho+fibre_frac*vEl*rhof
    
    dT = np.zeros(n)
    dV = np.zeros(n)
    
    A = np.zeros((n,n))


    #2n
    for i in range(1,n-1):
        A[i][i-1] = 1 #* (sA[i-1]/sA[i])
        A[i][i+1] = 1 #* (sA[i+1]/sA[i])
        A[i][i] = -2
        #old
        A[0][0] = -1   
        A[0][1] = 1 
        #new
#        A[0][0] = -11/6 
#        A[0][1] = 3
#        A[0][2] = -3/2
#        A[0][3] = 1/3
        
        A[-1][-1] = -1
        A[-1][-2] = 1
    
    for t in range(int(btime/dt)-1):
        #Heat transfer due to dT/dx
        dT = np.matmul(A,T[0])
        cpEff = (M[0]*cp(T[0])+(1-M[0])*cp2(T[0]))
        aEff = (M[0]*a(T[0])+(1-M[0])*a2(T[0]))
        dT *= (aEff*dt)/(dx*dx)
        #Heat used up in chemical reaction
        rR = Ac * (poly_frac) * (V[0]) * np.exp(-Ea/(R*T[0])) 
        dV = -rR * dt / rho
        dV[V[0]-dV<0] = 0
        dE = dV * vEl*poly_frac * rho * dH 
        dT2 = dE / (mEl[0] * cpEff)
        #boundary influx
        qIn[0] = hg * (Taw-T[0][0])
        #qIn[0] = 2E5
        dT[0] += (dt * qIn[0] * sA[0])/(mEl[0][0] * cpEff[0])
        #Radiation from outer surface
#        radPower = (emmisivity * sigma * sA[0]) * T[0][0] ** 4
#        dT[0] -= (radPower*dt)/(mEl[0][0] * cpEff[0])
        
        V[1] = V[0] + dV
        M[1] = (V[1] * rho ) / (V[1]*rho+(1-V[1])*rho2)
        mEl[1] = poly_frac  * vEl * (V[0]*rho + char_frac*(1-V[0])*rho2) + rhof * fibre_frac * vEl# 38% of the element is always fibre 
        T[1] = T[0] + dT + dT2
    #    #temperature averaging between virg and char
    #    cpEff2 = (M[1]*cp(T[t])+(1-M[1])*cp2(T[t]))
    #    dT3 = 0 #mEl[t]*cpEff*(dT+dT2)/(mEl[1]*cpEff2)-(dT-dT2)#-((dV*rho2 * poly_frac * vEl) * cp2(T[t]) * T[1])/(mEl[t]*cpEff+(dV*rho2 * poly_frac * vEl)*cp2(T[t]))
    #    T[1] = T[1] + dT3
        Ttemp = T
        mElTemp = mEl
        Vtemp = V
        Mtemp = M
        T[0] = Ttemp[1]
        V[0] = Vtemp[1]
        M[0] = Mtemp[1]
        mEl[0] = mElTemp[1]
    #print(T)
    if(disp):
        
        #ANIMATE
#        plt.figure(1)
#        x = np.arange(0,length,dx) * 100 #/ 2.54 #demon units
#        ln, = ax.plot(x,T[0],'r')
#        ln2, = ax2.plot(x,1-V[0])
#        
#        def anim(i):
#            ln.set_ydata(T[int(i*(btime/dt)*(mspf/1000)*(1/animtime))]) 
#            ln2.set_ydata(1-V[int(i*(btime/dt)*(mspf/1000)*(1/animtime))])
#            #plt.title("Throat region, t + "+str(i*(btime/animtime)*mspf/1000)+"s")
#            return ln,ln2,
#        
#        anim = animation.FuncAnimation(fig,anim,frames=int(animtime*1000/mspf-1),interval = mspf)
#        
#        ax2.set_ylim([-.01,1.01])
#        ax.set_ylim([np.min(T) - 50,np.max(T)+50])
#        ax.set_ylabel("Temperature(K)")
#        ax2.set_ylabel("Volume Fraction of Product")
#        ax.set_xlabel("Distance from inner surface (m)")
#        plt.show()
#        
        #PAPER COMPARISON
        #anim.save("Comparison.gif")
        #timearr = np.arange(0,btime,dt)
        #timearr *= heatRate
#        plt.figure(2)
#        plt.plot(timearr,V[:,0]*poly_frac,'r')
#        plt.plot(timearr,char_frac*(1-V[:,0])*poly_frac,'m')
#        plt.plot(timearr,(1-char_frac)*(1-V[:,0])*poly_frac,'g')
        #plt.plot(T[-1],cp(T[-1]),'r')
        
        #print(T)
        #print(V)
        #print(M)
        #print(mEl)
        print(A)
        print(dx**2/(2*a(min(T[0]))),">",dt,"?")
#        print("Max Temp dif", np.max(T[-1])-np.min(T[-1]))
#        print("Mass of polymer ablated is ", np.sum((1-V[-1])*( rho * vEl * poly_frac)), "kg")
#        print("Original polymer mass ", np.sum(mEl[0] - fibre_frac*vEl*rhof), "kg")
#        def dqEn(t):
#            return np.sum(mEl[t]*(T[t+1]-T[t])*(M[t]*cp(T[t]) + (1-M[t])*cp2(T[t])))
#        qEns = np.array([dqEn(a) for a in range(0,int(btime/dt)-2)])
#        #print(np.sum(mEl[-1]*(T[-1])*(M[-1]*cp(T[-1]) + (1-M[-1])*cp2(T[-1]))),"J inal thermal energy")
#        print("Energy into system is ", np.sum(qIn*sA[0]*dt), " J")
#        #print("Energy into system is ", np.sum(qIn*Area*dt), " J")
#        print("Thermal change of system is ", np.sum(qEns), " J")#np.trapz(qEns,range(0,int(btime/dt)-2)), " J")
#        #print(np.trapz(qEns,np.arange(0,btime,dt)))
#        print("Latent heat energy to vaporize material ", np.sum((1-V[-1])*vEl*poly_frac*rho*dH), " J")
    return T[-1],V[-1] 
    
#sweep different mesh sizes
#sweep = np.linspace(200,3000,30)
#xnet = []
#for i in sweep:
#    n = int(i)
#    X = runSim(.015,8593,2880,.0165,disp=False)
#    print(X)
#    print("\n Run " + str(i))
#    xnet.append(X)
#    


##Plot last value for this run

t = time.time()
T1,V1 = runSim(.015,8593,2880,.0165,disp=True)
print(time.time()-t)
X1 = np.linspace(0,.015,n)
T1 = T1 - 273.15

X1 = X1 * 10
plt.close('all')
fig,ax = plt.subplots()
ax2 = ax.twinx()
ln, = ax.plot(X1,T1,'r')
ln2, = ax2.plot(X1,1-V1)
ax2.set_ylim([-.01,1.01])
ax.set_ylim([np.min(T1) - 50,np.max(T1)+50])
ax.set_ylabel("Temperature(C)")
ax2.set_ylabel("Volume Fraction of Product")
ax.set_xlabel("Distance from inner surface (mm)")
plt.title("Engine Thermal Situation")

#x = np.linspace(0,.015,n)
#T2,F2 = runSim(.015,hg,Taw,0.0187,disp=True)
