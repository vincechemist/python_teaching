#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Python code provided as is.
Made by Vincent Wieczny, from Chemistry Department, ENS de Lyon, France
This code is under licence CC-BY-NC-SA. It enables you to reuse the code by mentioning the orginal author and without making profit from it.
This code is using widgets.py that you need to download in the same diretory as the main Python file.
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

#Redox system data

#Ce4+/Ce3+
E_std_Ce=1.44 #standard potential (V/SHE)
D_Ce3=1e-9
D_Ce4=1e-9
c0_Ce4=1e-3 #titrant concentration (mol/L)

#Fe3+/Fe2+
E_std_Fe=0.77 #standard potential (V/SHE)
D_Fe3=1e-9
D_Fe2=1e-9
c0_Fe2=1e-3 #titrated concentration (mol/L)

#Electrochemical setup and system
delta=1e-5 #diffuse layer thickness (m)
n=1 #number of exchanged electrons
A=1e-5 #electrode area (m²) 

#Titration parameters
V_Fe2=10 #titrated volume (mL)
V0=10 #total volume (mL)

## Modulated parameters

parameters = {'V' : widgets.FloatSlider(value=0.000001, description='$V$ $\mathrm{(mL)}$', min=0.00000001, max=20)}

#################
### Functions ###
#################

# Titration volume determination
def Veq(c0_Fe2,c0_Ce4,V_Fe2):
    return c0_Fe2*V_Fe2/c0_Ce4

# Fe2+ concentration as a function of the added volume V
def c_Fe2(V,Veq): 
    if V<=Veq:
        return (c0_Fe2*V_Fe2-c0_Ce4*V)/(V0)
    else:
        return 0.0000000001 # non-zero value due to a divergent behaviour

# Fe3+ concentration as a function of the added volume V
def c_Fe3(V,Veq):
    if V<Veq:
        return c0_Ce4*V/(V0)
    else:
        return c0_Ce4*Veq/(V0)
    
# Ce4+ concentration as a function of the added volume V
def c_Ce4(V,Veq):
    if V<Veq:
        return 0.0000000001
    else:
        return c0_Ce4*(V-Veq)/(V0)

# Ce3+ concentration as a function of the added volume V    
def c_Ce3(V,Veq):
    if V<Veq:
        return c0_Ce4*V/(V0)
    else:
        return c0_Ce4*Veq/(V0)

#Anodic diffusive controlled current
def i_a(delta,Cred,Dred):
  return n*F*A*Dred*Cred/delta

#Cathodic diffusive controlled current
def i_c(delta,Cox,Dox):
  return -n*F*A*Dox*Cox/delta

#Diffusion-limited current
def i_diff(delta,Dred,Dox,Cred,Cox,E,E_std):
    ia=i_a(delta,Cred,Dred)
    ic=i_c(delta,Cox,Dox)
    k=(E-E_std)*(n*F)/(R*T)
    return (np.exp(k)*ia+ic)/(1+np.exp(k))

#i-E Fe data at the added volume V:
def iE_data_Fe(V,Veq,E):
    #Concentration calculation
    Fe2=c_Fe2(V,Veq)
    Fe3=c_Fe3(V,Veq)
    
    #Curve calculation
    i_diff_Fe=i_diff(delta,D_Fe2,D_Fe3,Fe2,Fe3,E,E_std_Fe)

    return i_diff_Fe

#i-E data at the added volume V:
def iE_data_Ce(V,Veq,E):
    #Concentration calculation
    Ce4=c_Ce4(V,Veq)
    Ce3=c_Ce3(V,Veq)
    
    #Curve calculation
    i_diff_Ce=i_diff(delta,D_Ce3,D_Ce4,Ce3,Ce4,E,E_std_Ce)

    
    return i_diff_Ce

#i-E data at the added volume V:
def iE_data_tot(V,Veq,E):
    #Concentration calculation
    Fe2=c_Fe2(V,Veq)
    Fe3=c_Fe3(V,Veq)
    Ce4=c_Ce4(V,Veq)
    Ce3=c_Ce3(V,Veq)
    
    #Curve calculation
    i_diff_Fe=i_diff(delta,D_Fe2,D_Fe3,Fe2,Fe3,E,E_std_Fe)
    i_diff_Ce=i_diff(delta,D_Ce3,D_Ce4,Ce3,Ce4,E,E_std_Ce)
    i_diff_tot=i_diff_Fe+i_diff_Ce
    
    return i_diff_tot

#Titration spot
def titration_spot(V):
    i=iE_data_tot(V,Veq,E)
    counter=0
    for k in range(0,len(i)):
        if i[k]<0:
            counter=counter+1
    return E[counter]

#Titration curve
def titration_curve(V_domain):
    E_titr=[]
    for v in V_domain:
        E_titr.append(titration_spot(v))
    return E_titr






#===========================================================
# --- Initialization of the plot ---------------------------
#===========================================================

#fig,(ax1,ax2)=plt.subplots(1,2,figsize=(16,6))
fig=plt.figure(figsize=(18,6))

fig.suptitle(r'Potentiometric titration $(i=0)$  of a $\mathbf{Fe^{2+}}$ solution by a $\mathbf{Ce^{4+}}$ solution',weight='bold')


Veq=Veq(c0_Fe2,c0_Ce4,V_Fe2)

E=np.arange(0.0001,2.001,0.0011)
V_domain=np.arange(0.0000001,2*Veq+0.00001,0.01)



fig.text(0.01,0.9,r'Titration conditions', multialignment='left', verticalalignment='top',weight='bold')
fig.text(0.01,0.85,r'Titrated solution', multialignment='left', verticalalignment='top')
fig.text(0.01,0.80,r'$c_0=${:.3f} mol/L'.format(c0_Fe2), multialignment='left', verticalalignment='top')
fig.text(0.01,0.77,r'$V_0=${:.2f} mL'.format(V_Fe2), multialignment='left', verticalalignment='top')
fig.text(0.01,0.72,r'Titrant solution', multialignment='left', verticalalignment='top')
fig.text(0.01,0.67,r'$c_1=${:.3f} mol/L'.format(c0_Ce4), multialignment='left', verticalalignment='top')
fig.text(0.01,0.62,r'End point', multialignment='left', verticalalignment='top')
fig.text(0.01,0.57,r'$V_e=${:.2f} mL'.format(Veq), multialignment='left', verticalalignment='top')
fig.text(0.01,0.47,r'Dilution is not taken into account.', multialignment='left', verticalalignment='top')


ax1 = fig.add_axes([0.2, 0.2, 0.35, 0.7])
ax2 = fig.add_axes([0.60, 0.2, 0.35, 0.7])
#ax.axhline(0, color='k')



ax1.text(0.77,1.25e-7,'$\mathrm{Fe^{3+}_{(aq)} + e^- \leftrightarrows \, Fe^{2+}_{(aq)})}}$',horizontalalignment='center',
     verticalalignment='center')

ax1.text(1.44,1.25e-7,'$\mathrm{Ce^{4+}_{(aq)} + e^- \leftrightarrows \, Ce^{3+}_{(aq)})}}$',horizontalalignment='center',
     verticalalignment='center')    
 
ax1.plot([E.min(), E.max()],[0,0],':',lw=1,color='grey')    

ax1.set_xlim(E.min(), E.max())
ax1.set_ylim(-2e-7,2e-7)

ax1.set_xlabel('$E$ $\mathrm{(V/ESH)}$')
ax1.set_ylabel('Current $i$ $\mathrm{(A)}$')

ax2.plot([Veq,Veq],[0.25,1.75,],':',lw=1,color='grey')





ax2.set_xlim(V_domain.min(), V_domain.max())
ax2.set_ylim(0.25,1.75)

ax2.set_xlabel('$V$ $\mathrm{(mL)}$')
ax2.set_ylabel('$E$ $\mathrm{(V/ESH)}$')

#===========================================================
# --- Plot of the updated curves ---------------------------
#===========================================================

spot_list_V=[]
spot_list_E=[]
titration_curve_0=titration_curve(V_domain)



# This function is called when the sliders are changed 
def plot_data(V):
    spot_list_V.append(V)
    spot_list_E.append(titration_spot(V))
    lines['$i_\mathrm{Fe}$'].set_data(E,iE_data_Fe(V,Veq,E))
    lines['$i_\mathrm{Ce}$'].set_data(E,iE_data_Ce(V,Veq,E))
    lines['$i_\mathrm{tot}$'].set_data(E,iE_data_tot(V,Veq,E))
    lines['$Titration \ step \ by \ step$'].set_data(spot_list_V,spot_list_E)
    lines['$Titration \ curve$'].set_data(V_domain,titration_curve_0)
    truc['$Titration \ spot \ (left)$'].set_data(titration_spot(V),0)
    truc['$Titration \ spot \ (right)$'].set_data(V,titration_spot(V))
    fig.canvas.draw_idle()


#########################

lines = {}
lines['$i_\mathrm{Fe}$'], = ax1.plot([], [],color='green',lw=2,label='$i_\mathrm{Fe}$')
lines['$i_\mathrm{Ce}$'], = ax1.plot([], [],color='blue',lw=2,label='$i_\mathrm{Ce}$')
lines['$i_\mathrm{tot}$'], = ax1.plot([], [], lw=3, color='red',label='$i_\mathrm{tot}$')
lines['$Titration \ step \ by \ step$'],=ax2.plot([],[],'s',color='black',lw=1)
lines['$Titration \ curve$'],=ax2.plot(V_domain,titration_curve_0,color='red',lw=3,label='$Titration \ curve$')
truc={}
truc['$Titration \ spot \ (left)$'], = ax1.plot([], [],'o',color='black',lw=2,label='$Titration \ spot$')
truc['$Titration \ spot \ (right)$'], = ax2.plot([], [],'o',color='black',lw=2,label='$Titration \ spot$')


ax1.legend()

ax2.legend()

param_widgets = widgets.make_param_widgets(parameters, plot_data, slider_box=[0.20, 0.05, 0.35, 0.05])
choose_widget = widgets.make_choose_plot(lines, box=[0.01,0.2,0.12, 0.2])
reset_button = widgets.make_reset_button(param_widgets,box=[0.85, 0.05, 0.10, 0.05])

if __name__=='__main__':
    plt.show()