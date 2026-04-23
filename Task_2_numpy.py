# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 21:19:08 2026

@author: jeanb
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

### Constant

Lx=60
Ly=60
#Nx=200
#Ny=200
dx=dy=0.3
Nx = round(Lx/dx)
Ny = round(Ly/dy)

dt=0.1
tfinal=10000
Nt=int(tfinal//dt)

T0=10
Q=0.2
bc=[0,0,0,0] # border (right, up, left, down)

borehole=[30,30]
x0, y0 = borehole

ix = int(x0 / dx)
iy = int(y0 / dy)

bc=[0,0,0,0] # border (right, up, left, down)

def vx(x):
    return 0.02*(x +20)+0.01 

def vy(y):
    return 0.02*(y-20) 

x_coords = np.arange(Nx) * dx
y_coords = np.arange(Ny) * dy

vx_field = vx(x_coords)[:, None]
vy_field = vy(y_coords)[None, :]

### Initialisation


Told=np.zeros((Nx,Ny))+T0
T=np.zeros((Nx,Ny))

target_times = [2000, 5000, 10000]
results = [None, None, None]
 
#Calculating the temperature profile timestep by timestep
for k in range(Nt):

    
    # upwind x
    dTdx = np.where(
        vx_field > 0,
        (Told - np.roll(Told, 1, axis=0)) / dx,
        (np.roll(Told, -1, axis=0) - Told) / dx
    )
    
    # upwind y
    dTdy = np.where(
        vy_field > 0,
        (Told - np.roll(Told, 1, axis=1)) / dy,
        (np.roll(Told, -1, axis=1) - Told) / dy
    )
    
    T = Told - dt * (vx_field * dTdx + vy_field * dTdy)
            
    T[ix, iy] += Q * dt
    
    T[0,  :]  = T[1,  :]    # bord gauche
    T[-1, :]  = T[-2, :]    # bord droit
    T[:,  0]  = T[:,  1]    # bord bas
    T[:, -1]  = T[:, -2]    # bord haut

    Told=T.copy()

    if k%10000==0:
            print("t =", k*dt, "s")
        #print(T)
    
    t_current = int(round(k * dt))

    if t_current in target_times:
        idx = target_times.index(t_current)
        results[idx] = T.copy()


# plots
x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
Y, X = np.meshgrid(y, x)

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

levels = None

for i in range(3):
    T_plot = results[i]

    if levels is None:
        vmin = min(r.min() for r in results)
        vmax = max(r.max() for r in results)
        levels = np.linspace(vmin, vmax, 21)

    ax = axes[i]
    surf = ax.contourf(X, Y, T_plot, levels=levels, cmap=cm.coolwarm)

    ax.set_title(f"t = {target_times[i]} s")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

fig.colorbar(surf, ax=axes.ravel().tolist(),
             orientation='horizontal', pad=0.15, shrink=0.9)

plt.show()
