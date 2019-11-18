# -*- coding: utf-8 -*-

#T in Kelvin, density in g/ml ,, works only for liquid phase
import numpy as np
import sympy as sp
mmass = 16.042460
Tc = 190.55
Tt = 90.68
pc = 10
pt = 28.147
A = -0.178860165
C = -0.01848987
B = 0.04838475 
E = 0.36 
def density(T):
    w = (Tc-T)/(Tc-Tt)
    return mmass/1000 * (pc + (pt-pc)*(w**E*np.exp(A*(1-w**(2/3))+B*(1-w**(4/3))+C*(1-w**2))))



#returns cp in j/gk (only works for gas phase .-.) 
R = 0.831434*10**1
b5 = 4.7207907
b1 = 2.5998981 
b6 = 5.02288
b2 = 1.4449418
b3 = -1.8472716
b4 = 0.8211218 
def cp(T):
    q=T/400
    return R/mmass * (1+b1 + (4/3)*b2*q**(1/3) + (5/3)*b3*q**(2/3) + 2*b4*q + b5*np.exp(b6/q) * (b6/(q*(np.exp(b6/q)-1)))**2)


#returns cp in j/gk (only works for gas phase .-.) [innacurate where t is more than 220k]
m1 = -1.8044750507*10**6
m2 = 7.7426666393*10**4
m3 = -1.3241658754*10**3
m4 = 1.5438149595*10**1
m5 = -5.1479005257*10**-2
m6 = 1.0809172196*10**-4
m7 = -6.5501783437*10**-8
m8 = -6.7490056171
m9 = 3.0000000000*10**3
R = 0.831434*10**1
#T in kelvin and cp in k/mol k
def cp2(T):
    M = [m1,m2,m3,m4,m5,m6,m7]
    c = 0.0
    for i in range(7):
        c+= M[i]*T**(i-3)
    c += 7*(m8*(m9**2*np.exp(m9/T))/(T**2*((np.exp(m9/T))-1)**2))
    c*= R/mmass
    return c

n1 = 0.9898937956 *10**-5
n2 = 0.2199608275 *10**-1
n3 = -0.5322788000
n4 = 0.2021657962 *10**2
n5 = -0.2234398926 *10**4
n6 = 0.1067940280 *10**-4
n7 = 0.1457922469 *10**-3
n8 = -0.9265816666
n9 = 0.2915364732 *10**3
n10 = 0.2313546209 *10**-6
n11 = 0.1387214274 *10**-3
n12 = 0.4780467451 *10**-2
n13 = 0.1176103833 *10**-4
n14 = -0.1982096730 *10**-3
n15 = -0.2512887756 *10**-1
n16 = 0.9748899826 *10**-5
n17 = -0.1202192137 *10**-6
n18 = 0.4128353939 *10**-4
n19 = -0.7215842918 *10**-6
n20 = 0.5081738255 *10**3
n21 = -0.9198903192 *10**5
n22 = -0.2732264677 *10**1
n23 = 0.7499024351 *10**5
n24 = 0.1114060908 *10**-2
n25 = 0.1083955159 *10**1
n26 = -0.4490960312 *10**-4
n27 = -0.1380337847 *10**1
n28 = -0.2371902232 *10**-7
n29 = 0.3761652197 *10**-4
n30 = -0.2375166954 *10**-9
n31 = -0.1237640790 *10**-7
n32 = 0.6766926453 *10**-6
R = 0.831434*10**-2
rhoc = 10.139342719
y = 1/rhoc**2
#input in g/cm3 and K
def pressureP(T,p):
    pr = 0.0;
    p *= (1000/mmass)
    A = [R*T, n1*T+n2*T**.5+n3+n4/T+n5/T**2, n6*T+n7+n8/T+n9/T**2,
         n10*T+n11+n12/T,n13,n14/T+n15/T**2,n16/T,n17/T+n18/T**2,n19/T**2]
    B =  [n20/T**2+n21/T**3,n22/T**2+n23/T**4,n24/T**2+n25/T**3,n26/T**2+n27/T**4,
         n28/T**2+n29/T**3,n30/T**2+n31/T**3+n32/T**4]
    for i in range(9):
        pr+= A[i]*p**(i+1)
        #print(pr)
    #print("b")
    bt = 0.0
    for j in range(6):
        bt+= B[j]*(p**(2*(j+1)+1))
    pr += bt*(np.exp(-y*p**2))
    return pr
def denP(T,P):
    mn,mx = .3,1
    e = 1000
    g =  mn
    for i in range(mn,mx,.001):
        err = P
    return g