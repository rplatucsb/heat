# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 15:29:48 2019

@author: Adam
"""
    
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

plt.rcParams["animation.ffmpeg_path"] = r"C:\Users\Adam\Desktop\RPL\Heet\ffmpeg"
k = 400 #W/m k
cp = 300 #J/kg k
rho = 9000 #kg/m^3
cpc = 1000
rhoc = 4000
kc = 10
H = 100000 #j/kg
elements = 200
dx = 1/elements #m
dt = .05#s
tNet = 2000#s
t = 0
area = 1
fps = 10 #really 1000/n

a = np.arange(50,150,100/elements)
#a = 273+100*np.sin(np.linspace(0,2*np.pi))
#a.fill(100)
#a = np.random.random_integers(100,high=200,size = (elements))
#a[25] = 200

elements = int(np.size(a))
mFracV = np.zeros((elements))
mFracV.fill(1)

rArr = np.zeros((int(tNet/dt)+1,elements,2))
rArr[0][:,0] = a
rArr[0][:,1] = mFracV
fig,ax = plt.subplots()

while (t < tNet):
    titCount = int(t/dt)
    #right dt array is shifting ts to the left
    darr = np.array(rArr[titCount][:,0])
    darr = np.delete(darr,0)
    darr = np.append(darr,rArr[titCount][:,0][-1])
    dTarr = darr - rArr[titCount][:,0]
    # left dt array is shifting ts to right
    darl = np.array(rArr[titCount][:,0])
    darl = np.delete(darl,-1)
    darl = np.insert(darl,0,rArr[titCount][:,0][0])
    dTarl = darl - rArr[titCount][:,0]
    
    Js = k * dt * (dTarl + dTarr) * area
    #Js[0] += 100000
    #test q that wants to flow to material, find m vap, stop at 0, change k based on m
#    tdTs = Js/(area * dx * rho * cp)
#    overT = 200-(rArr[titCount][:,0] + tdTs) 
#    overT[overT<0] = 0
#    jGone = overT*rho*dx*area*cp
#    dM = jGone/H
#    Js -= jGone
#    Js[rArr[titCount][:,1] - dM<0]=0
    
    dT = Js/(area * dx * rho * cp)
    
    dT = k/(rho*cp)*(darr+darl-2*rArr[titCount][:,0])/(dx*dx) * dt
    
    rArr[titCount+1][:,0] = rArr[titCount][:,0] + dT
    rArr[titCount+1][:,1] = rArr[titCount][:,1] #- dM
    rArr[rArr[:,:,0]<0]=0
#    rArr[:,0] = 100
#    rArr[:,-1] = 100
    t += dt

x = np.arange(0,elements*dx,dx)
ln, = ax.plot(x,rArr[0][:,0])

def anim(i):
    ln.set_ydata(rArr[int(i*(tNet/dt)*(fps/1000)*(1/5))][:,0]) #last denom is how long anim will last 
    return ln,

anim = animation.FuncAnimation(fig,anim,interval = fps)

plt.xlabel("X (m)")
plt.ylabel("Temp (K)")
plt.show()


#anim.save("snek.gif")


print(np.size(rArr))
print(rArr)
a= (k/(rho*cp))
print(dx**2/(2*a),">",dt,"?")