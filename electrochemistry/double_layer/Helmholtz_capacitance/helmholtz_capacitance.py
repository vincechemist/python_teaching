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
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons

################################
### Paramater initialization ###
################################

#Physical constants
F=96500.0 #Faraday number (C/mol)
R=8.314 #Gas constant (J/K/mol)
T=298.0 #Temperature (K)

#Physical parameters
epsilon_0=8.85418782e-12 #vacuum permittivity (F/m)
epsilon_r=78.5

#Electrolyte
z=1




## Modulated parameters

parameters = {'deltaphi' : widgets.FloatSlider(value=0, description='$\phi_\mathrm{m}-\phi_\mathrm{sol}$ $\mathrm{(V)}$', min=-0.5, max=0.5),
              'xH' : widgets.FloatSlider(value=2, description='$x_\mathrm{H}$ $\mathrm{(nm)}$', min=1.5, max=5),
         
              }

#################
### Functions ###
#################

#Helmholtz capacitor (in µF/cm²)
def C_H(E,xH):
    return epsilon_0*epsilon_r/(xH*1e-9)*1e2


#Potential at a distance x from the electrode
def phi(x,xH,deltaphi):
    return deltaphi*(1-x/xH)

#===========================================================
# --- Plot of the updated curves ---------------------------
#===========================================================

# This function is called when the sliders are changed 
def plot_data(xH,deltaphi):
    lines['$C_\mathrm{H}$'].set_data(E_domain,C_H(E_domain,xH))
    lines['$phi_m$'].set_data([-2e8,0],[deltaphi,deltaphi])
    lines['$phi1$'].set_data([0,xH],[deltaphi,0])
    lines['$phi2$'].set_data([xH,15],[0,0])
    lines['$xH$'].set_data([xH,xH],[-1,1])
    fig.canvas.draw_idle()


#===========================================================
# --- Initialization of the plot ---------------------------
#===========================================================


fig=plt.figure(figsize=(18,6))

fig.suptitle(r'Helmholtz electrochemical double-layer',weight='bold')


E_domain=np.arange(-0.5000001,0.5000001,0.0011)
x_domain=np.arange(0,15,0.1)



fig.text(0.01,0.9,r'Helmholtz model', multialignment='left', verticalalignment='top',weight='bold')
fig.text(0.01,0.85,r'Helmholtz capacitance', multialignment='left', verticalalignment='top')
fig.text(0.01,0.78,r'$C_\mathrm{H}=\dfrac{\epsilon}{x_\mathrm{H}}$')
fig.text(0.01,0.73,r'Compact double-layer characterised by', multialignment='left', verticalalignment='top')
fig.text(0.01,0.70,r'the Helmholtz length $x_\mathrm{H}$', multialignment='left', verticalalignment='top')
fig.text(0.01,0.67,r'where the potential of the solution $\phi(x)$ is', multialignment='left', verticalalignment='top')
fig.text(0.01,0.60,r'$\phi(x)=\left( \, \phi_\mathrm{m}-\phi_\mathrm{sol} \, \right) \, \left( \, 1 - \dfrac{x}{x_\mathrm{H}} \, \right) + \phi_\mathrm{sol}$ for $x<x_\mathrm{H}$', multialignment='left', verticalalignment='top')
fig.text(0.01,0.55,r'$\phi(x)=\phi_\mathrm{sol}$ for $x>x_\mathrm{H}$', multialignment='left', verticalalignment='top')
fig.text(0.01,0.50,r'with', multialignment='left', verticalalignment='top')
fig.text(0.03,0.47,r'- $\phi_\mathrm{m}$, the metal phase potential', multialignment='left', verticalalignment='top')
fig.text(0.03,0.43,r'- $\phi_\mathrm{sol}$, the bulk solution potential', multialignment='left', verticalalignment='top')




ax1 = fig.add_axes([0.23, 0.2, 0.35, 0.7])
ax2 = fig.add_axes([0.63, 0.2, 0.35, 0.7])







    



ax1.set_xlim(E_domain.min(), E_domain.max())
ax1.set_ylim(0,50)

ax1.set_xlabel('$\phi_\mathrm{m}-\phi_\mathrm{sol}$ $\mathrm{(V)}$')
ax1.set_ylabel('$C_\mathrm{H}$ $\mathrm{(\mu F/cm^2)}$')


ax2.add_artist(matplotlib.patches.Rectangle((-0.2*100,-1),0.2*100,2, color = 'grey'))


ax2.set_xlim(-0.2*x_domain.max(),x_domain.max())
ax2.set_ylim(-0.55,0.55)

ax2.plot([0,100],[0,0],':',color='gray')

ax2.set_xlabel('$x$ $\mathrm{(nm)}$')
ax2.set_ylabel('$\phi(x)-\phi_\mathrm{sol}$ $\mathrm{(V)}$')


lines = {}

lines['$C_\mathrm{H}$'], = ax1.plot([], [],color='red',lw=3)
lines['$phi1$'], = ax2.plot([], [],color='red',lw=3)
lines['$phi_m$'], = ax2.plot([], [],color='red',lw=3)
lines['$phi2$'], = ax2.plot([], [],color='red',lw=3)
lines['$xH$'], = ax2.plot([], [],'--',color='gray',lw=2)
ax1.legend()





param_widgets = widgets.make_param_widgets(parameters, plot_data, slider_box=[0.23, 0.05, 0.35, 0.05])
#choose_widget = widgets.make_choose_plot(lines, box=[0.01,0.2,0.12, 0.2])
reset_button = widgets.make_reset_button(param_widgets,box=[0.88, 0.05, 0.10, 0.05])

if __name__=='__main__':
    plt.show()
