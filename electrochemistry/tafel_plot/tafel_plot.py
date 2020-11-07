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


#Electrochemical system
n=1.0 #number of exchanged electrons
Estd=0 #redox standard potential
Dox=1.0e-9 #oxidant diffusion coefficient (m2/s)
Dred=1.0e-9 #reductant diffusion coefficient (m2/s)
i0=1e-11 #exchange current (A)
alpha=0.5 #charge transfer coefficient


#Electrode
A=0.00001 #electrode area (m2)
delta=1e-5 #diffusive layer (m)

#Bulk solution
Cox=1.0e-3 #initial oxidant bulk concentration (mol/L)
Cred=1.0e-3 #initial reductant bulk concentration (mol/L)


#################
### Functions ###
#################

#Nernst potential 
def E_Nersnt(Estd,Cred,Cox):
    return Estd+(R*T)/(n*F)*np.log(Cox/Cred)

#Anodic diffusive controlled current
def i_a(delta,Cred,Dred):
  return n*F*A*Dred*Cred/delta

#Cathodic diffusive controlled current
def i_c(delta,Cox,Dox):
  return -n*F*A*Dox*Cox/delta

#Diffusion-limited current
def i_diff(delta,Dred,Dox,Cred,Cox,eta):
    ia=i_a(delta,Cred,Dred)
    ic=i_c(delta,Cox,Dox)
    return (1-np.exp(-n*F*eta/(R*T)))/(1/ia-np.exp(-n*F*eta/(R*T))/ic)

#Charge-transfer limited current
def i_charge(i0,alpha,eta):
    return i0*(np.exp(alpha*n*F*eta/(R*T))-np.exp(-(1-alpha)*n*F*eta/(R*T)))

#Total current
def i_tot(idiff,icharge):
    return (idiff*icharge)/(idiff+icharge)

#Log(total current)
def log_i(i):
    return np.log(abs(i))



######################
### Graphical data ###
######################
  
#Graph definition
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,5))    

#Overpotential domain
ENersnt=E_Nersnt(Estd,Cred,Cox)
eta=np.arange(ENersnt-1,ENersnt+1.01,0.001)

#Initial calculations

#Initial i-E curves
idiff=i_diff(delta,Dred,Dox,Cred,Cox,eta)
icharge=i_charge(i0,alpha,eta)
itot=i_tot(idiff,icharge)

#Initial Tafel plots
logidiff=log_i(idiff)
logicharge=log_i(icharge)
logitot=log_i(itot)

#Initial overpotential spot
etavalue_init=0.001
idiff_eta=i_diff(delta,Dred,Dox,Cred,Cox,etavalue_init)
icharge_eta=i_charge(i0,alpha,etavalue_init)
itot_eta=i_tot(idiff_eta,icharge_eta)
logitot_eta=log_i(itot_eta)


#Graph initialisation

#Left graph

#Initial i-E curves
curve_idiff,=ax1.plot(eta,idiff,':',label=r'Diffusive current $i_{diff}$',color='grey',lw=2)
curve_icharge,=ax1.plot(eta,icharge,'--',label=r'Charge transfer limited current $i_{ct}$',color='grey',lw=2)
curve_itot,=ax1.plot(eta,itot,label=r'Total current $i_{tot}$',color='red',lw=3)

#Itinial overpotential spot
etapt,=ax1.plot(etavalue_init,itot_eta,'o',color='black',lw=3)

#Axis label
ax1.set_xlabel(r'$\eta$ $\mathrm{(V)}$')
ax1.set_ylabel(r'$i$ $\mathrm{(A)}$')

#Axis limits
ax1.set_xlim(min(eta),max(eta))
ax1.set_ylim(min(itot)*1.2,max(itot)*1.2)

ax1.legend()

#Right graph

#Initial Tafel plots
curve_logidiff,=ax2.plot(eta,logidiff,':',label=r'Diffusive current $i_{diff}$',color='grey',lw=2)
curve_logicharge,=ax2.plot(eta,logicharge,'--',label=r'Charge transfer limited current $i_{ct}$',color='grey',lw=2)
curve_logitot,=ax2.plot(eta,logitot,lw=3,color='red',label=r'Total current $i_{tot}$')

#Itinial overpotential spot
logetapt,=ax2.plot(etavalue_init,logitot_eta,'o',color='black')

#Axis label
ax2.set_xlabel(r'$\eta$ $(V)$')
ax2.set_ylabel(r'$\log{(i)}$ (with $i$ in $A$)')

#Axis limits
ax2.set_xlim(min(eta),max(eta))
ax2.set_ylim(0.6*min(logitot),max(logitot)*0.8)

ax2.legend()

###############
### Cursors ###
###############

#Axis definitions

#Overpotential axis
axETA = plt.axes([0.20, 0, 0.2, 0.025])

#Charge transfer coefficient cursor
axALPHA = plt.axes([0.65, 0, 0.2, 0.025])


#Slider definitions

#Overpotential slider
ETA=Slider(axETA, r'$\eta$ $(V)$', min(eta), max(eta), valinit=etavalue_init,color='grey')

#Charge transfer coefficient slider
ALPHA=Slider(axALPHA, r'$\alpha$', 0.1, 0.9, valinit=alpha,color='grey')


#Graph update
def update(val):
    #Updated values
    etavalue=ETA.val
    alphavalue=ALPHA.val
    
        #Updated calculations
    idiff_new=i_diff(delta,Dred,Dox,Cred,Cox,eta)
    icharge_new=i_charge(i0,alphavalue,eta)
    itot_new=i_tot(idiff_new,icharge_new)
    
    logidiff_new=log_i(idiff_new)
    logicharge_new=log_i(icharge_new)
    logitot_new=log_i(itot_new)    
    
    idiff_pt=i_diff(delta,Dred,Dox,Cred,Cox,etavalue)
    icharge_pt=i_charge(i0,alphavalue,etavalue)
    itot_pt=i_tot(idiff_pt,icharge_pt)
    logitot_pt=log_i(itot_pt)   
    
    #Updated curves
    curve_idiff.set_ydata(idiff_new)
    curve_icharge.set_ydata(icharge_new)
    curve_itot.set_ydata(itot_new)
    
    curve_logidiff.set_ydata(logidiff_new)
    curve_logicharge.set_ydata(logicharge_new)
    curve_logitot.set_ydata(logitot_new)
    
    etapt.set_xdata(etavalue)
    etapt.set_ydata(itot_pt)
        
    logetapt.set_xdata(etavalue)
    logetapt.set_ydata(logitot_pt)
    
    #Graph refresh
    fig.canvas.draw_idle()

#Call update function on slider value change
ETA.on_changed(update)
ALPHA.on_changed(update)

plt.show()
pass