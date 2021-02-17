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

################################
### Paramater initialization ###
################################

#Physical constants
F=96500.0 #Faraday number (C/mol)
R=8.314 #Gas constant (J/K/mol)
T=298.0 #Temperature (K)

# Modulated parameters
parameters = {'C' : widgets.FloatSlider(value=10, description='$C$ $\mathrm{(\mu F)}$', min=1, max=100),
              'Re' : widgets.FloatSlider(value=100, description='$R_\mathrm{e}$ $\mathrm{(\Omega)}$', min=0, max=1000),
              'Rtc' : widgets.FloatSlider(value=1000, description='$R_\mathrm{tc}$ $\mathrm{(\Omega)}$', min=100, max=10000)}

#################
### Functions ###
#################

def ReZ(logf,C,Re,Rtc):
    C2=C*1e-6
    f=10**logf
    omega=2*np.pi*f
    return Re+Rtc/(1+(omega**2)*(Rtc**2)*(C2**2))

def ImZ(logf,C,Re,Rtc):
    C2=C*1e-6
    f=10**logf
    omega=2*np.pi*f
    return omega*(Rtc**2)*C2/(1+(omega**2)*(Rtc**2)*(C2**2))

#===========================================================
# --- Plot of the updated curves ---------------------------
#===========================================================

# This function is called when the sliders are changed 
def plot_data(C,Re,Rtc):
    lines['$EIS$'].set_data(ReZ(logf,C,Re,Rtc),ImZ(logf,C,Re,Rtc))
    
    Zmax=max(ReZ(logf,C,Re,Rtc).max(),ImZ(logf,C,Re,Rtc).max())
    ax1.set_xlim(0,Zmax)
    ax1.set_ylim(0,Zmax)
    fig.canvas.draw_idle()


#===========================================================
# --- Initialization of the plot ---------------------------
#===========================================================


fig=plt.figure(figsize=(10,10))



logf=np.arange(-3,8,0.1)




ax1=fig.add_axes([0.2, 0.2, 0.7, 0.7])




ax1.set_xlim(0,1500)
ax1.set_ylim(0,1500)

ax1.set_xlabel('$\mathrm{Re(Z)}$ $(\Omega)$')
ax1.set_ylabel('$-\mathrm{Im(Z)}$ $(\Omega)$')




lines = {}
lines['$EIS$'], = ax1.plot([], [],'o',color='red',lw=3)

ax1.legend()




param_widgets = widgets.make_param_widgets(parameters, plot_data, slider_box=[0.20, 0.05, 0.35, 0.05])
#choose_widget = widgets.make_choose_plot(lines, box=[0.01,0.2,0.12, 0.2])
reset_button = widgets.make_reset_button(param_widgets,box=[0.85, 0.05, 0.10, 0.05])

if __name__=='__main__':
    plt.show()
