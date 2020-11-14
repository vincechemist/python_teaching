#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Python code provided as is.
Made by Vincent Wieczny, from Chemistry Department, ENS de Lyon, France.
This code is under licence CC-BY-NC-SA. It enables you to reuse the code by mentioning the orginal author and without making profit from it.

Objective. The aim of this code is to illustrate what a hard pulse means. Even if the pulse is non-resonant due to a non-zero offset, the chemical shift modulation results in a slighty modification of the resultant magnetization.

How to.
Modulate theta angle to change the angle of rotation about Beff field
Modulate delta to change Beff from a resonant to a non-resonant pulse.
"""

#Librairies
import qutip as qt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider

####################
### Initial data ###
####################

v0=500e6 #spectrometer 1H frequency (Hz)
delta_init=5 #resonant pulse chemical shift (ppm)
t_pulse=12e-6 #Pi/2 resonant pulse duration (s)
theta_init=90 #initial rotation angle (°)

#################
### Functions ###
#################

#Offset pulsation
def offset(v0,delta):
    return 2*np.pi*v0*(delta-delta_init)*1e-6

#B1 field pulsation
def B1(theta,t_pulse):
    return (theta*np.pi/180)/t_pulse

#Beff angle (°) with z-axis with Beff in (Oyz)
def Beff_angle(B1,offset):
    if offset==0:
        return 90
    else:
        if B1/offset>0:
            return np.arctan(B1/offset)*180/np.pi
        if B1/offset<0:
            return np.arctan(B1/offset)*180/np.pi+180

#B1 axis definition
def B1_axis(angle):
    return np.array([0,np.sin(angle*np.pi/180),0])

#Delta Omega axis definition
def Delta_Omega_axis(angle):
    return np.array([0,0,np.cos(angle*np.pi/180)])

#Beff axis definition in (Oyz)
def Beff_axis(angle):
    return np.array([0,np.sin(angle*np.pi/180),np.cos(angle*np.pi/180)])

#Magnetization rotation
def rotation(rot_axis,vec_init,theta):
    c=np.cos(theta*np.pi/180)
    s=np.sin(theta*np.pi/180)
    x=rot_axis[0]
    y=rot_axis[1]
    z=rot_axis[2]
    rot_mat=np.array([[x**2*(1-c)+c,x*y*(1-c)-z*s,x*z*(1-c)+y*s],
	              [x*y*(1-c)+z*s,y**2*(1-c)+c,y*z*(1-c)-x*s],
	              [x*z*(1-c)-y*s,y*z*(1-c)+x*s,z**2*(1-c)+c]])
    return np.dot(rot_mat,vec_init)

################
### Graphics ###
################

#Graphic initialization
b=qt.Bloch()
b.xlabel=[r'$x^\prime$', '']
b.ylabel=[r'$y^\prime$', '']
b.zlabel=[r'$z$', '']
b.xlpos = b.ylpos = b.zlpos = [1.3, 1.3]
b.view=[-50,10]
b.figsize=[10,10] 
b.sphere_alpha=0.05 #sphere transparency

#Vector initialization
M0=np.array([0,0,1]) #Initial magnetization

#Beff vector initial calculation
offset_init=offset(v0,delta_init)
B1_init=B1(theta_init,t_pulse)
angle_init=Beff_angle(B1_init,offset_init)

Beff_init=Beff_axis(angle_init) #Initial Beff 


M_init=rotation(Beff_init,M0,theta_init) #Initial resulting magnetization

#Add initial vectors
b.add_vectors(M0)
b.add_annotation(M0*1.1,r'$\overrightarrow{M}_0$',color='black')

b.add_vectors(Beff_init)
b.add_annotation(Beff_init*1.1,r'$\overrightarrow{B}_\mathrm{eff}$',color='orange')

b.add_vectors(M_init)
b.add_annotation(M_init*1.1,r'$\overrightarrow{M}$',color='black')

b.vector_color=['black','orange','black']
b.vector_width=5


#Graph display
b.show()

#################
### Animation ###
#################

#Cursor definition

#Impulsion angle axis
axtheta= plt.axes([0.1, 0.05, 0.3, 0.025])

#Chemical shift delta
axdelta=plt.axes([0.6, 0.05, 0.3, 0.025])

#Implusion angle slider
THETA=Slider(axtheta, r'$\theta$ $(^\circ )$', 0, 180, valinit=theta_init,color='grey')

DELTA=Slider(axdelta, r'$\delta$ $(ppm)$', 0, 10, valinit=delta_init,color='grey')


def update(val):
    #Updated values
    theta=THETA.val
    delta=DELTA.val

    #Graph reinitialization
    b.clear()
    
    #Rotation axis new coordinates
    offset_update=offset(v0,delta)
    B1_update=B1(theta,t_pulse)
    angle_update=Beff_angle(B1_update,offset_update)
    
    B1_vec_update=B1_axis(angle_update)
    Delta_Omega_vec_update=Delta_Omega_axis(angle_update)
    Beff_update=Beff_axis(angle_update)
        
    #Magnetization new coordinates
    M=rotation(Beff_update,M0,theta)

    #Updated vectors
    b.add_vectors(M0)
    b.add_annotation(M0*1.1,r'$\overrightarrow{M}_0$',color='black')
    
    b.add_vectors(Beff_update)
    b.add_annotation(Beff_update*1.1,r'$\overrightarrow{B}_\mathrm{eff}$',color='orange')
    
    b.add_vectors(B1_vec_update)
    if np.abs(delta-delta_init)>2.5:
        b.add_annotation(B1_vec_update*1.1,r'$\overrightarrow{B}_1$',color='orange')
    else:
        pass
    
    b.add_vectors(Delta_Omega_vec_update)
    if np.abs(delta-delta_init)>2.5:
        b.add_annotation(Delta_Omega_vec_update*1.5,r'$\overrightarrow{\Delta B}$',color='orange')
    else:
        pass
    
    b.add_vectors(M)
    b.add_annotation(M*1.1,r'$\overrightarrow{M}$',color='black')
    
    b.vector_color=['black','orange','orange','orange','black']
    b.vector_width=5
    
    #Graph refresh
    b.show()

#Call update function on slider value change
THETA.on_changed(update)
DELTA.on_changed(update)



