# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 16:18:11 2026

@author: jeanb
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
 
### Constant

Lx=60
Ly=60
#Nx=200
#Ny=200
dx=0.3
dy=0.3
Nx = round(Lx/dx)
Ny = round(Ly/dy)

dt=0.15
tfinal=10000
Nt=int(tfinal//dt)

T0=10
alpha=5e-3

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

K = alpha*dt/(dx*dx)
vx2 = vx_field*dt/dx
vy2 = vy_field*dt/dy

vx_pos = np.maximum(vx2, 0) 
vx_neg = np.minimum(vx2, 0)
vy_pos = np.maximum(vy2, 0)
vy_neg = np.minimum(vy2, 0)

Q= 0.3*np.random.randn(Nt)

### Initialisation
target_times = [2000, 5000, 10000]
results = []

Told=np.zeros((Nx,Ny))+T0
T=np.zeros((Nx,Ny))
 
#Calculating the temperature profile timestep by timestep
for k in range(Nt):

    
    T = Told*(1-4*K) + K*(np.roll(Told, 1, axis=0) + np.roll(Told, -1, axis=0) 
                          + np.roll(Told, 1, axis=1) + np.roll(Told, -1, axis=1))
    
    
    T = T - vx_pos * (Told - np.roll(Told, 1, axis=0))  # Dérivée arrière si vx > 0
    T = T - vx_neg * (np.roll(Told, -1, axis=0) - Told)  # Dérivée avant si vx < 0
    

    T = T - vy_pos * (Told - np.roll(Told, 1, axis=1))
    T = T - vy_neg * (np.roll(Told, -1, axis=1) - Told)
    
  
    T[ix, iy] += Q[k] * dt
    
    T[0,  :]  = T[1,  :]    # bord gauche
    T[-1, :]  = T[-2, :]    # bord droit
    T[:,  0]  = T[:,  1]    # bord bas
    T[:, -1]  = T[:, -2]    # bord haut

    Told=T.copy()

    if int(k*dt) in target_times:
        results.append(T.copy())    

    if k%100==0:
        print(k*dt)
        #print(T)

        
###Plotting the result
"""
x=np.linspace(0,Lx,Nx)
y=np.linspace(0,Ly,Ny)
Y, X = np.meshgrid(y, x)

T_levels=np.linspace(9.9,np.amax(T),21)
fig, ax=plt.subplots() 
surf=ax.contourf(X,Y,T,levels=T_levels,cmap=cm.coolwarm)
bar=fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.show()
"""
###     Echelle de couleur logarithmique
import numpy as np

x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
Y, X = np.meshgrid(y, x)

titles = ["t = 2 000 s", "t = 5 000 s", "t = 10 000 s"]

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for i, t_val in enumerate(target_times):
    dT = results[i] - T0
    surf = axes[i].contourf(X, Y, dT,
                            levels=15,
                            cmap='coolwarm')
    axes[i].set_title(f't = {t_val} s')
    axes[i].set_xlabel('x [m]')
    axes[i].set_ylabel('y [m]')
cbar = fig.colorbar(surf, ax=axes.ravel().tolist(),
                    orientation='horizontal',
                    pad=0.21,
                    shrink=0.9)

cbar.set_label('T - T₀ [°C]')
cbar.ax.xaxis.set_label_position('top')
cbar.ax.xaxis.set_ticks_position('top')
plt.tight_layout()
plt.show()

