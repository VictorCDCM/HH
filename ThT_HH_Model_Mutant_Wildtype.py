#Emergence of ion channel mediated signaling in E. coli biofilm

# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 09:47:40 2022

@author: Victor , Emmanuel
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from numpy import savetxt
import csv

plt.rcParams['axes.labelsize'] = 25
plt.rcParams['xtick.labelsize'] = 25
plt.rcParams['ytick.labelsize'] = 25
#plt.rcParams['axes.facecolor'] = '#f2f2f2'
plt.rcParams['axes.facecolor'] = '#ffffff'


gk = 90                   #Potassium conductance 
gl = 0.2                  #Leak channel conductance 
gq = 7                    #New channel conductance 
Vk = -94                  #Nernst potential for potassium
Vq = -200                 #Nernst potential for the new channel
Vl = -156                 #Nernst potenatial for the leak channel 
Sth = 0.04                #Stress threshold for channel opening
Vth = -150                #voltage threshold for stress production
a0 = 2                    #Maximum rate of channel opening   
b0 = 1.3                  #channel rate decay constant
m = 1                     #cooperativity parameter 
sa = 0.2                  #stress threshold sharpness coefficient
ge = 10                   #extracellular potassium relaxation factor
dl = 8                    #leak slope coefficient
dq = 2			  #q channel slope coefficient
dk = 1                    #potassium slope coefficient
bs = 0.001                #Tht uptake rate coefficient                 
qs = 0.5                  #Tht uptake rate coefficient for the new channel 
gs = 0.1                  #stress decay rate
gt = 4                    #Tht decay constant                           
F = 5.6                   #membrane capacitance 
kshock = 200              #potassium shock




Sthm = Sth**m
Vss = Vk+ge*gl*(Vk-Vl)/(-ge*gl)

def f(vne,t):
    V, n, E, S, T, N, A = vne
    a = a0*S**m/(Sthm+S**m)
    b = b0
    Sin = (bs*(Vth-V)/(np.exp((Vth-V)/sa)-1)) 
    Qin = (qs*(Vq-V)/(np.exp((Vq-V)/sa)-1))
    Tin = bs*(Vss-V)
    Vls = Vl+dl*E
    Vqs = Vq+dq*E
    Vks = Vk+dk*E
    #dVdt = -gq*n**4*(Vq-Vks)-gk*n**4*(V-Vks) - gl*(V-Vls)
    dVdt = -gq*n**4*(V-Vqs)-gk*n**4*(V-Vks) - gl*(V-Vls)  
    dndt = a*(1-n)-b*n
    #dEdt = F*gq*gk*n**4*(V-Vks) - ge*E                      
    dEdt = F*gk*n**4*(V-Vks) - ge*E
    dSdt = Sin - gs*S
    dTdt = Tin - gt*T + Qin
    dNdt = bs*S - gs*N 
    dAdt = bs*E - gs*A
       
    return [dVdt, dndt, dEdt, dSdt, dTdt, dNdt, dAdt]


#wild-type Ecoli ion-channel mediated signaling 

Tmax = 350
dt = 0.01
vne0 = [1, 1, 1, 0, 0, 0, 0]
tvec0 = np.arange(0,Tmax,dt)
vne0_out = odeint(f, vne0, tvec0)

Tmax = 60
tvec = np.arange(0,Tmax,dt)
vne0 = [vne0_out[-1,0], vne0_out[-1,1], vne0_out[-1,2]+kshock, vne0_out[-1,3], vne0_out[-1,4], vne0_out[-1,5], vne0_out[-1,6]]
vne_out = odeint(f, vne0, tvec)

thtts = np.concatenate((vne0_out[:,4],vne_out[:,4]),axis=0)
normf = np.amax(abs(thtts[20000:-1]))
thtts = thtts*40/normf
tts = np.concatenate((tvec0[:],tvec[:]+tvec0[-1]),axis=0)-350 



# Kch_DH5alpha mutant  

a0 = 2         #maximum rate of channel opening quite high, faster than the kch wildtype 
b0 = 3       #reduced decay constant for the new channel 
gq = 0.01   #The channel conductance is quite low and typically habituates
gs = 0.1
qs = 0.0009
gk = 0 

Tmax = 350
dt = 0.01
vne0 = [1, 1, 1, 0, 0, 0, 0]
tvec0 = np.arange(0,Tmax,dt)
vne0_out = odeint(f, vne0, tvec0)

Tmax = 60
tvec = np.arange(0,Tmax,dt)
vne0 = [vne0_out[-1,0], vne0_out[-1,1], vne0_out[-1,2]+kshock, vne0_out[-1,3], vne0_out[-1,4], vne0_out[-1,5], vne0_out[-1,6]]
vne_out = odeint(f, vne0, tvec)

thtts2 = np.concatenate((vne0_out[:,4],vne_out[:,4]),axis=0)
normf = np.amax(abs(thtts2[20000:-1]))
thtts2 = thtts2*40/normf
tts2 = np.concatenate((tvec0[:],tvec[:]+tvec0[-1]),axis=0)   

ns =  np.concatenate((vne0_out[:,1],vne_out[:,1]),axis=0)
normn = np.amax(abs(ns[20000:-1]))
ns = ns*40/normn

fig = plt.figure(figsize=(9,6))
plt.plot(tts,thtts,'b',linewidth=8)
plt.plot(tts2,thtts2,'k',linewidth=8)
plt.xlabel('Time (mins)')
plt.ylabel('ThT (arb. units)')
plt.xlim(0,60)
plt.ylim(1,45)
#plt.grid()



#with open('##INSERT DIRECTORY PATH HERE##', 'w', encoding='UTF8', newline='') as f:
#    writer = csv.writer(f)

    # write the header
#    writer.writerow(thtts)

#   writer.writerow(tts)
