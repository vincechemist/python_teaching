#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Python code provided as is.
Made by Vincent Wieczny, from Chemistry Department, ENS de Lyon, France
This code is under licence CC-BY-NC-SA. It enables you to reuse the code by mentioning the orginal author and without making profit from it.
"""

#Librairies
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider

#Physical constants
F=96500.0 #Faraday number (C/mol)
R=8.314 #Gas constant (J/K/mol)
T=298.0 #Temperature (K)
nu=1.007e-6 #cinematic viscosity of water (m2/s)

#Electrochemical system
n=1.0 #number of exchanged electrons
Estd=0 #redox standard potential
Dox=1.0e-9 #oxidant diffusion coefficient (m2/s)
Dred=1.0e-9 #reductant diffusion coefficient (m2/s)

#Rotating disk electrode
A=0.00001 #electrode area (m2)
w_init=500 #initial rotational speed (rpm)
w_min=100 #minimal rotational speed (rpm)
w_max=1000 #maximal rotational speed (rpm)

#Bulk solution
Cox=1.0e-3 #initial oxidant bulk concentration (mol/L)
Cred=1.0e-3 #initial reductant bulk concentration (mol/L)
Cox_max=2.0e-3 #maximal oxidant bulk concentration (mol/L)
Cred_max=2.0e-3 #maximal reductant bulk concentration (mol/L)

#################
### Functions ###
#################

#Converting rpm to rad/s
def convert_w(w): 
  return w*2*np.pi/60

#Diffuse layer thickness calculation
def delta(w,D):
  return 1.61*w**(-1/2)*nu**(1/6)*D**(1/3)

#Anodic diffusive controlled current
def i_a(delta_red,Cred):
  return n*F*A*Dred*Cred/delta_red

#Cathodic diffusive controlled current
def i_c(delta_ox,Cox):
  return -n*F*A*Dox*Cox/delta_ox

#Half-wave potential
def Ehalfwave(delta_red,delta_ox,Dox,Dred):
  return Estd+(R*T)/(n*F)*np.log((Dred/delta_red)/(Dox/delta_ox))

#i-E curve
def i(w,Dred,Dox,Cred,Cox,E):
  w=convert_w(w)
  delta_red=delta(w,Dred)
  delta_ox=delta(w,Dox)
  ia=i_a(delta_red,Cred)
  ic=i_c(delta_ox,Cox)
  Ehv=Ehalfwave(delta_red,delta_ox,Dox,Dred)
  k=(E-Ehv)*(n*F)/(R*T)
  return (np.exp(k)*ia+ic)/(1+np.exp(k))

######################
### Graphical data ###
######################
  
#Graph definition
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,5))

#Potential domain
E=np.arange(Estd-0.25,Estd+0.26,0.01)

#Initial calculations

#Initial i-E curve
i_init=i(w_init,Dred,Dox,Cred,Cox,E)
#Nernst potential (i=0)
E_init=Estd+R*T/(n*F)*np.log(Cox/Cred) 
#Initial half-wave potential and current 
delta_red_init=delta(w_init,Dred)
delta_ox_init=delta(w_init,Dox)
E_hw_init=Ehalfwave(delta_red_init,delta_ox_init,Dox,Dred)
i_hw_init=(max(i_init)+min(i_init))/2
#Initial reductant concentration profile
k_red_init=Cred/max(i_init)
Cred_0_init=Cred-k_red_init*i(w_init,Dred,Dox,Cred,Cox,E_init)
#Initial oxidant concentration profile
k_ox_init=Cox/min(i_init)
Cox_0_init=Cox-k_ox_init*i(w_init,Dred,Dox,Cred,Cox,E_init)

#Graph initialisation

#Left graph

#x-axis definition
x_axis1,=ax1.plot([min(E),max(E)],[0,0],color='black')

#Diffusion-limited currents initialisation
imax_axis,=ax1.plot([min(E),max(E)],[max(i_init),max(i_init)],'--',color='black',label=r'$i_{lim}$')
imin_axis,=ax1.plot([min(E),max(E)],[min(i_init),min(i_init)],'--',color='black')

#i-E curve initialisation
iEcurve, = ax1.plot(E,i_init,lw=3,color='red',label=r'$i-E$ curve')

#Half-wage potential and current initialisation
xE_hw,=ax1.plot([E_hw_init,E_hw_init],[i_hw_init,0],':',color='black',label=r'$E_{1/2}$ and $i_{1/2}$')
yE_hw,=ax1.plot([min(E),E_hw_init],[i_hw_init,i_hw_init],':',color='black')

#Potential sweeping
ptE,=ax1.plot(E_init,i(w_init,Dred,Dox,Cred,Cox,E_init),'o',color='black')

#Axis label
ax1.set_xlabel(r'$E_{WE}$ $\mathrm{(V/SHE)}$')
ax1.set_ylabel(r'$i$ $\mathrm{(A)}$')

#Axis limits
ax1.set_xlim(min(E),max(E))
ax1.set_ylim(1.05*min(i(w_max,Dred,Dox,Cred_max,Cox_max,E)),1.05*max(i(w_max,Dred,Dox,Cred_max,Cox_max,E))),

ax1.legend()



#Right graph

#x-axis definition
x_axis2,=ax2.plot([0,5e-5],[0,0],color='black')

#Diffusive layer initialisation
xdelta_red,=ax2.plot([delta_red_init,delta_red_init],([0,Cred]),'--',color='blue',label=r'$Red$ $diffusive$ $layer$ $limit$ $\delta_{Red}$')
xdelta_ox,=ax2.plot([delta_ox_init,delta_ox_init],([0,Cox]),'--',color='green',label=r'$Ox$ $diffusive$ $layer$ $limit$ $\delta_{Ox}$')

#Convective concentration profile initialisation
profile_conv_red,=ax2.plot([delta_red_init,5e-5],[Cred,Cred],lw=3,color='blue',label=r'$c_{Red}(x)$')
profile_conv_ox,=ax2.plot([delta_ox_init,5e-5],[Cox,Cox],lw=3,color='green',label=r'$c_{Ox}(x)$')

#Diffusive concentration profile initialisation
profile_diff_red,=ax2.plot([0,delta_red_init],[Cred_0_init,Cred],color='blue',lw=3)
profile_diff_ox,=ax2.plot([0,delta_ox_init],[Cox_0_init,Cox],color='green',lw=3)

#Electrode (x<0)
ax2.add_artist(matplotlib.patches.Rectangle((-0.2*5e-5,-0.2*max(Cred_max,Cox_max)),0.2*5e-5,2.2*max(Cred_max,Cox_max), color = 'gray'))

#Axis label
ax2.set_xlabel(r'$Distance$ $x$ $from$ $the$ $electrode$ $\mathrm{(m)}$')
ax2.set_ylabel(r'$Concentration$ $\mathrm{(mol/L)}$')

#Axis limits
ax2.set_xlim(-0.2*max(delta(w_min,Dred),delta(w_min,Dox)),1.2*max(delta(w_min,Dred),delta(w_min,Dox)))
ax2.set_ylim(-0.2*max(Cred_max,Cox_max),2*max(Cred_max,Cox_max))

ax2.legend()

###############
### Cursors ###
###############

#Axis definitions

#Rotational speed axis
axw = plt.axes([0.20, 0, 0.2, 0.025])

#Potential sweeping axis
axESW=plt.axes([0.20,0.03, 0.2, 0.025])

#Bulk reductant concentration axis
axCred = plt.axes([0.65,0, 0.2, 0.025])

#Bulk oxidant concentration axis
axCox = plt.axes([0.65,0.03, 0.2, 0.025])

#Slider definitions

#Rotational speed slider
W=Slider(axw, r'$Rotational$ $speed$ $\mathrm{(rpm)}$', w_min, w_max, valinit=w_init,color='grey')

#Potential sweeping slider
ESW=Slider(axESW,r'$E_{WE}$ $\mathrm{(V/SHE)}$',min(E),max(E),valinit=E_init,color='grey')

#Bulk reductant concentration slider
CRED=Slider(axCred,r'$c_{Red}$ $\mathrm{(mol/L)}$',0,Cred_max,valinit=Cred,color='grey')

#Bulk oxidant concentration slider
COX=Slider(axCox,r'$c_{Ox}$ $\mathrm{(mol/L)}$',0,Cox_max,valinit=Cox,color='grey')

#Graph update
def update(val):
    #Updated values
    w=W.val
    Cred=CRED.val
    Cox=COX.val
    Esw=ESW.val
    
    #Updated calculations
    i_new=i(w,Dred,Dox,Cred,Cox,E)
    
    delta_red=delta(w,Dred)
    delta_ox=delta(w,Dox)
    
    E_hw=Ehalfwave(delta_red,delta_ox,Dox,Dred)
    i_hw=(max(i_new)+min(i_new))/2
    
    k_red=Cred/max(i_new)
    Cred_0=Cred-k_red*i(w,Dred,Dox,Cred,Cox,Esw)
    
    k_ox=Cox/min(i_new)
    Cox_0=Cox-k_ox*i(w,Dred,Dox,Cred,Cox,Esw)
    
    #Updated curves
    iEcurve.set_ydata(i_new)
    
    xE_hw.set_xdata([E_hw,E_hw])
    xE_hw.set_ydata([i_hw,0])
    yE_hw.set_xdata([min(E),E_hw])
    yE_hw.set_ydata([i_hw,i_hw])
    imax_axis.set_ydata([max(i_new),max(i_new)])
    imin_axis.set_ydata([min(i_new),min(i_new)])
    
    ptE.set_xdata(Esw)
    ptE.set_ydata(i(w,Dred,Dox,Cred,Cox,Esw))
    
    xdelta_red.set_xdata([delta_red,delta_red])
    xdelta_red.set_ydata([0,Cred])
    xdelta_ox.set_xdata([delta_ox,delta_ox])
    xdelta_ox.set_ydata([0,Cox])
    
    profile_conv_red.set_xdata([delta_red,5e-5])
    profile_conv_red.set_ydata([Cred,Cred])
    profile_conv_ox.set_xdata([delta_ox,5e-5])
    profile_conv_ox.set_ydata([Cox,Cox])
    
    profile_diff_red.set_xdata([0,delta_red])
    profile_diff_red.set_ydata([Cred_0,Cred])
    profile_diff_ox.set_xdata([0,delta_ox])
    profile_diff_ox.set_ydata([Cox_0,Cox])
    
    #Graph refresh
    fig.canvas.draw_idle()

#Call update function on slider value change
W.on_changed(update)
CRED.on_changed(update)
COX.on_changed(update)
ESW.on_changed(update)


plt.show()
pass
