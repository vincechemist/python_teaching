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
z=3


## Modulated parameters

parameters = {'deltaphi' : widgets.FloatSlider(value=0, description='$\phi_\mathrm{m}-\phi_\mathrm{sol}$ $\mathrm{(V)}$', min=-0.5, max=0.5),
              'logc' : widgets.FloatSlider(value=-1, description='$\log{(c)}$ $\mathrm{(c \ in \ mol/L)}$', min=-6, max=0),
              'z' : widgets.IntSlider(value=1, description='$z:z$', min=1, max=3)             
              }

#################
### Functions ###
#################

#Gouy-Chapman capacitor (in µF/cm²)
def C_GC(E,logc,z):
    c=10**logc
    return z*F*np.sqrt(2*epsilon_0*epsilon_r*c/R/T)*np.cosh(z*F/(2*R*T)*E)*1e2


#Debye length (in m)
def lD(c,z):
    return 1/(z*F*np.sqrt(2*c/(epsilon_0*epsilon_r*R*T)))


#Potential at a distance x from the electrode
def phi(x,logc,deltaphi,z):
    c=10**logc
    X=x*1e-9
    lDebye=lD(c,z)
    return deltaphi*exp(-X/lDebye)

#===========================================================
# --- Plot of the updated curves ---------------------------
#===========================================================

# This function is called when the sliders are changed 
def plot_data(logc,deltaphi,z):
    lines['$C_\mathrm{GC}$'].set_data(E_domain,C_GC(E_domain,logc,z))
    lines['$phi$'].set_data(x_domain,phi(x_domain,logc,deltaphi,z))
    
    lines['$phi_m$'].set_data([-2e8,0],[deltaphi,deltaphi])
    lines['spot'].set_data(deltaphi,C_GC(deltaphi,logc,z))
    lines['Debye'].set_data([0,lD(10**logc,z)*1e9],[deltaphi,0])
    fig.canvas.draw_idle()


#===========================================================
# --- Initialization of the plot ---------------------------
#===========================================================


fig=plt.figure(figsize=(18,6))

fig.suptitle(r'Gouy-Chapman electrochemical double-layer',weight='bold')


E_domain=np.arange(-0.5000001,0.5000001,0.0011)
x_domain=np.arange(0,100,0.1)



fig.text(0.01,0.9,r'Gouy-Chapman model', multialignment='left', verticalalignment='top',weight='bold')
fig.text(0.01,0.85,r'Gouy-Chapman capacitance', multialignment='left', verticalalignment='top')
fig.text(0.01,0.75,r'$C_\mathrm{GC}=z \, F \, \sqrt{\dfrac{2 \, \epsilon \, c}{R \, T}} \ \cosh{ \, \left[ \, \dfrac{z \, F}{R \, T} \, \left( \phi_\mathrm{m} - \phi_\mathrm{sol} \right)\, \right]}$')
fig.text(0.01,0.70,r'Diffusive double-layer characterised by', multialignment='left', verticalalignment='top')
fig.text(0.01,0.67,r'the Debye length $l_\mathrm{D}$', multialignment='left', verticalalignment='top')
fig.text(0.01,0.57,r'$l_\mathrm{D}=\dfrac{1}{z \, F} \, \sqrt{ \dfrac{\epsilon \, R \, T}{2 \, c}}$')
fig.text(0.01,0.52,r'where the potential of the solution $\phi(x)$ is', multialignment='left', verticalalignment='top')
fig.text(0.01,0.47,r'$\phi(x)=\left( \, \phi_\mathrm{m}-\phi_\mathrm{sol} \, \right) \, \exp{\left( \, - \, \dfrac{x}{l_\mathrm{D}} \, \right)}$', multialignment='left', verticalalignment='top')
fig.text(0.01,0.40,r'with', multialignment='left', verticalalignment='top')
fig.text(0.03,0.40,r'- z, the electrolyte charge', multialignment='left', verticalalignment='top')
fig.text(0.03,0.37,r'- c, the electrolyte concentration', multialignment='left', verticalalignment='top')
fig.text(0.03,0.34,r'- $\phi_\mathrm{m}$, the metal phase potential', multialignment='left', verticalalignment='top')
fig.text(0.03,0.31,r'- $\phi_\mathrm{sol}$, the bulk solution potential', multialignment='left', verticalalignment='top')

ax1 = fig.add_axes([0.23, 0.2, 0.35, 0.7])
ax2 = fig.add_axes([0.63, 0.2, 0.35, 0.7])


ax1.set_xlim(E_domain.min(), E_domain.max())
ax1.set_ylim(0,50)

ax1.set_xlabel('$\phi_\mathrm{m}-\phi_\mathrm{sol}$ $\mathrm{(V)}$')
ax1.set_ylabel('$C_\mathrm{GC}$ $\mathrm{(\mu F/cm^2)}$')


ax2.add_artist(matplotlib.patches.Rectangle((-0.2*100,-1),0.2*100,2, color = 'grey'))


ax2.set_xlim(-0.2*x_domain.max(),x_domain.max())
ax2.set_ylim(-0.55,0.55)

ax2.plot([0,100],[0,0],':',color='gray')

ax2.set_xlabel('$x$ $\mathrm{(nm)}$')
ax2.set_ylabel('$\phi(x)-\phi_\mathrm{sol}$ $\mathrm{(V)}$')


lines = {}
lines['Debye'], = ax2.plot([], [],'--',color='gray',lw=2)
lines['$C_\mathrm{GC}$'], = ax1.plot([], [],color='red',lw=3)
lines['$phi$'], = ax2.plot([], [],color='red',lw=3)
lines['$phi_m$'], = ax2.plot([], [],color='red',lw=3)
lines['spot'], = ax1.plot([], [],'o',color='black',lw=2)

ax1.legend()





param_widgets = widgets.make_param_widgets(parameters, plot_data, slider_box=[0.23, 0.02, 0.35, 0.075])
#choose_widget = widgets.make_choose_plot(lines, box=[0.01,0.2,0.12, 0.2])
reset_button = widgets.make_reset_button(param_widgets,box=[0.88, 0.05, 0.10, 0.05])

if __name__=='__main__':
    plt.show()
