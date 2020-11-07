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

#Electrochemical system
n=1.0 #number of exchanged electrons

####################
###Initialisation###
####################

E_Oxe_init=0 #Initial Gibbs free energy of oxidant + n electrons
E_Red_init=0 #Initial Gibbs free energy of reductant
E_tr_init=50 #Initial activation Gibbs free energy (symmetric case)
alpha_init=0.5 #Initial charge transfer coeffient
DeltaE_init=-0.25 #Initial electrode potential difference (relative to standard potential)
pente_init=E_tr_init/(15*2*alpha_init) 


DeltaGm_init=-n*F*DeltaE_init/1000

x_tr_init=(E_tr_init-(1-alpha_init)*DeltaGm_init)/pente_init-15
y_tr_init=E_Oxe_init+E_tr_init+alpha_init*DeltaGm_init


fig,ax=plt.subplots()
plt.xlim(-25,25)

niv_E_Oxe_init,=plt.plot([-25,-15],[E_Oxe_init,E_Oxe_init],lw=3,color='black')

plt.text(-20,-20,r'$Ox+ne^-$',horizontalalignment='center',verticalalignment='center',color='black',fontsize=12)
plt.text(20,-20,r'$Red$',horizontalalignment='center',verticalalignment='center',color='orange',fontsize=12)

niv_E_Red_init,=plt.plot([25,15],[E_Red_init,E_Red_init],lw=3,color='orange')

niv_E_Oxe,=plt.plot([-25,-15],[E_Oxe_init+DeltaGm_init,E_Oxe_init+DeltaGm_init],'--',lw=3,color='black')

ch_Oxe_tr_std,=plt.plot([-15,15*(2*alpha_init-1)],[E_Oxe_init,E_Oxe_init+E_tr_init],lw=2,color='black')
ch_Red_tr_std,=plt.plot([15,15*(2*alpha_init-1)],[E_Red_init,E_Red_init+E_tr_init],lw=2,color='orange')

ch_Oxe_tr,=plt.plot([-15,x_tr_init],[E_Oxe_init+DeltaGm_init,y_tr_init],'--',lw=2,color='black')

ch_Red_tr,=plt.plot([15,x_tr_init],[E_Red_init,y_tr_init],'--',lw=2,color='orange')


ecart_1,=plt.plot([-18.33,-18.33],[E_Oxe_init,E_Oxe_init+DeltaGm_init],color='black',lw=2,label=r'$\Delta G_m =-nF(E-E^\circ)$')

ecart_2,=plt.plot([-21.66,-21.66],[E_Oxe_init+DeltaGm_init,E_Oxe_init+DeltaGm_init+E_tr_init-(1-alpha_init)*DeltaGm_init],'-.',color='black',lw=2,label=r'$\Delta^\ddag G_c$')

ecart_3,=plt.plot([20,20],[E_Red_init,E_tr_init+alpha_init*DeltaGm_init],'-.',color='orange',lw=2,label=r'$\Delta^\ddag G_a$')

ecart_4,=plt.plot([-15,-15],[E_tr_init,E_tr_init-(1-alpha_init)*DeltaGm_init],':',color='black',lw=2,label=r'$(1-\alpha) \Delta G_m$')

ecart_5,=plt.plot([15,15],[E_tr_init,E_tr_init+alpha_init*DeltaGm_init],':',color='orange',lw=2,label=r'$\alpha \Delta G_m$')



bar_act_c,=plt.plot([-25,25],[E_tr_init-(1-alpha_init)*DeltaGm_init,E_tr_init-(1-alpha_init)*DeltaGm_init],':',color='gray')
bar_act_std,=plt.plot([-25,25],[E_tr_init,E_tr_init],':',color='black')
bar_act_a,=plt.plot([-25,25],[E_tr_init+alpha_init*DeltaGm_init,E_tr_init+alpha_init*DeltaGm_init],':',color='gray')

plt.ylim(-100,150)
ax.set_xticklabels([])


plt.ylabel(r'$G_\mathrm{m}$ $\mathrm{(kJ \cdot mol^{-1})}$')
plt.legend()

axE = plt.axes([0.25, 0.03, 0.2, 0.025])
axalpha = plt.axes([0.60, 0.03, 0.2, 0.025])
E=Slider(axE, r'$E-E^\circ$ $\mathrm{(V)}$',-1,1, valinit=DeltaE_init,color='red')
Alpha=Slider(axalpha, r'$\alpha$',0,1, valinit=alpha_init,color='red')

def update(val):
    # amp is the current value of the slider
    pot=E.val
    alpha=Alpha.val
    DeltaGm=-n*F*pot/1000
    pente=E_tr_init/(15*2*alpha)
    x_tr=(E_tr_init-(1-alpha)*DeltaGm)/pente-15
    y_tr=E_Oxe_init+E_tr_init+alpha*DeltaGm
    
    
    niv_E_Oxe.set_ydata([E_Oxe_init+DeltaGm,E_Oxe_init+DeltaGm])
    
    
    bar_act_c.set_ydata([E_tr_init-(1-alpha)*DeltaGm,E_tr_init-(1-alpha)*DeltaGm])
    bar_act_a.set_ydata([E_tr_init+alpha*DeltaGm,E_tr_init+alpha*DeltaGm])
  
    ch_Oxe_tr_std.set_xdata([-15,15*(2*alpha-1)])
    ch_Red_tr_std.set_xdata([15,15*(2*alpha-1)])
 
    ch_Oxe_tr.set_xdata([-15,x_tr])
    ch_Oxe_tr.set_ydata([E_Oxe_init+DeltaGm,y_tr])
    
    ch_Red_tr.set_xdata([15,x_tr])
    ch_Red_tr.set_ydata([E_Red_init,y_tr])
    
 
    ecart_1.set_ydata([E_Oxe_init,E_Oxe_init+DeltaGm])
    
    ecart_2.set_ydata([E_Oxe_init+DeltaGm,E_Oxe_init+DeltaGm+E_tr_init-(1-alpha)*DeltaGm])

    ecart_3.set_ydata([E_Red_init,E_tr_init+alpha*DeltaGm])
    
    ecart_4.set_ydata([E_tr_init,E_tr_init-(1-alpha)*DeltaGm])
    
    ecart_5.set_ydata([E_tr_init,E_tr_init+alpha*DeltaGm])
    # redraw canvas while idle
    fig.canvas.draw_idle()

# call update function on slider value change
E.on_changed(update)
Alpha.on_changed(update)


plt.show()

