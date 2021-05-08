#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Python code provided as is.
Made by Vincent Wieczny, from Chemistry Department, ENS de Lyon, France
This code is under licence CC-BY-NC-SA. It enables you to reuse the code by mentioning the orginal author and without making profit from it.
"""

#Librairies
import matplotlib.pyplot as plt
import numpy as np
import widgets
import scipy.constants as constants
from matplotlib import rc
import matplotlib.patches as patches
import matplotlib as mpl

################################
### Paramater initialization ###
################################

#Physical constants
h=6.62607004e-34 #Planck constant (m2.kg.s-1)
R=8.314 #Gas constant (J/K/mol)
T=298.0 #Temperature (K)
eV=1.60e-19 #1 eV in J

#EPR physical constants
g=2.0023 #Landé g-factor
muB=9.274009994e-24 #Bohr magneton (J.T-1)
nu=9388.2e6 #X-band frequency (Hz)

#EPR magnetic field
Bmin=0 #T
Bmax=0.5 #T
DeltaBmax=1e-1 #T

# Modulated parameters
parameters = {'DeltaB' : widgets.FloatSlider(value=0.05, description='$B_1$ $\mathrm{(T)}$', min=0.03, max=DeltaBmax),
              'B0' : widgets.FloatSlider(value=0.1, description='$B_0$ $\mathrm{(T)}$', min=Bmin, max=Bmax)}

#################
### Functions ###
#################

#Down-state energy 
def E_down(B0):
    return -0.5*g*muB*B0

#Up-state energy 
def E_up(B0):
    return 0.5*g*muB*B0

#Transition energy 
def E_trans():
    return h*nu

#Resonant magnetic field
def B_trans():
    return E_trans()/(g*muB)

#Sigma
def sigma(DeltaB):
    return DeltaB/6


def signal_abs(B0,DeltaB):
    return 1/(sigma(DeltaB)*np.sqrt(2*3.1416))*np.exp(-(B0-B_trans())**2/(2*sigma(DeltaB)**2))

#Derivative signal
def signal_der(B0,DeltaB):
    return 1/(sigma(DeltaB)*np.sqrt(2*3.1416))*-(B0-B_trans())/sigma(DeltaB)**2*signal_abs(B0,DeltaB)

#===========================================================
# --- Plot of the updated curves ---------------------------
#===========================================================


## This function is called when the sliders are changed 
def plot_data(B0,DeltaB):
    
    lines['Absorption spot'].set_data(B0,signal_abs(B0,DeltaB))
    lines['First derivative spot'].set_data(B0,signal_der(B0,DeltaB))
    truc['$Abs_courbe$'].set_data(B,signal_abs(B,DeltaB))
    truc['$Der_courbe$'].set_data(B,signal_der(B,DeltaB))
    truc['$E_\mathrm{trans}$'].set_data([B0,B0],[-1,1])
    r1.set_transform(mpl.transforms.Affine2D().translate(B0-DeltaB/2,-E_trans()/2)+ax1.transData)
    r1.set_width(DeltaB)
       
    fig.canvas.draw_idle()


##===========================================================
## --- Initialization of the plot ---------------------------
##===========================================================

#Plot definition
fig=plt.figure(figsize=(18,8))

ax1 = fig.add_axes([0.2, 0.2, 0.35, 0.7])
ax2 = fig.add_axes([0.60, 0.6, 0.35, 0.3])
ax3 = fig.add_axes([0.60, 0.2, 0.35, 0.3])


#Plot comments
fig.suptitle(r'Simulation of an EPR spectrum at X waveband for a free electron',weight='bold')

fig.text(0.01,0.9,r'EPR magnetic field', multialignment='left', verticalalignment='top',weight='bold')
fig.text(0.01,0.85,r'$B=B_0+B_1 \ \cos{(2 \, \pi \, \nu \, t)}$', multialignment='left', verticalalignment='top')
fig.text(0.01,0.82,r'with $\nu=100 \ \mathrm{kHz}$', multialignment='left', verticalalignment='top')
fig.text(0.01,0.77,r'EPR X band frequency', multialignment='left', verticalalignment='top',weight='bold')
fig.text(0.01,0.72,r'$\nu_\mathrm{X}=9388.2 \ \mathrm{MHz}$', multialignment='left', verticalalignment='top')
fig.text(0.01,0.67,r'EPR spin level energies', multialignment='left', verticalalignment='top',weight='bold')
fig.text(0.01,0.62,r'Up-state', multialignment='left', verticalalignment='top')
fig.text(0.01,0.59,r'$E_\mathrm{up}=\frac{1}{2} \, g \, \mu_\mathrm{B} \, B$', multialignment='left', verticalalignment='top')
fig.text(0.01,0.54,r'Down-state', multialignment='left', verticalalignment='top')
fig.text(0.01,0.51,r'$E_\mathrm{down}=-\frac{1}{2} \, g \, \mu_\mathrm{B} \, B$', multialignment='left', verticalalignment='top')


B=np.arange(Bmin,Bmax,0.0005)




if __name__=='__main__':
    
    ax1.plot(B,E_up(B),lw=2,color='red',label='Up-state energy')
    ax1.plot(B,E_down(B),lw=2,color='blue',label='Down-state energy')
    ax1.plot([Bmin,Bmax],[E_trans()/2,E_trans()/2],':',lw=2,color='grey',label='X-band energy')    
    ax1.plot([Bmin,Bmax],[-E_trans()/2,-E_trans()/2],':',lw=2,color='grey')    
    
    
    ax1.set_xlim(Bmin,Bmax)
    ax1.set_xlabel('$B_0$ $\mathrm{(T)}$')
    ax1.set_ylabel('$E$ $\mathrm{(J)}$')
    
    
    ax2.set_xlim(Bmin,Bmax)
    ax2.set_ylim(-10,150)
    ax2.set_yticklabels([])
    ax2.set_xlabel('$B_0$ $\mathrm{(T)}$')
    ax2.set_ylabel('$Absorption \ intensity$')
    
    
    ax3.set_xlim(Bmin,Bmax)
    ax3.set_ylim(-1000000,1000000)
    ax3.set_yticklabels([])
    ax3.set_xlabel('$B_0$ $\mathrm{(T)}$')
    ax3.set_ylabel('$First \ derivative \ intensity$')
    
    truc={}
    
    truc['$E_\mathrm{trans}$'], = ax1.plot([], [],'--',lw=2, color='gray',label='$B_0$')
    truc['$Abs_courbe$'], = ax2.plot([],[],lw=2,color='red',label='Absorption signal')
    truc['$Der_courbe$'], = ax3.plot([],[],lw=2,color='red',label='First derivative signal')
    r1 = ax1.add_patch(patches.Rectangle((0, 0),DeltaBmax,E_trans(), edgecolor = '#000000', facecolor = '#dddddd', fill=True,label='Excitation band'))

    lines = {}

    lines['Absorption spot'], = ax2.plot([],[],'o',color='black',lw=2)
    lines['First derivative spot'], = ax3.plot([],[],'o',color='black',lw=2)
    
    ax1.legend()
    ax2.legend()
    ax3.legend()
 
    param_widgets = widgets.make_param_widgets(parameters, plot_data, slider_box=[0.20, 0.05, 0.35, 0.05])
    choose_widget = widgets.make_choose_plot(lines,box=[0.01,0.2,0.12, 0.1])
    reset_button = widgets.make_reset_button(param_widgets,box=[0.85, 0.05, 0.10, 0.05])
    
    plt.show()