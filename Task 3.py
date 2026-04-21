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

dt=0.1
tfinal=5000
Nt=int(tfinal//dt)

T0=10
Q=5
alpha=5e-3

bc=[0,0,0,0] # border (right, up, left, down)

borehole=[30,30]
x0, y0 = borehole

ix = int(x0 / dx)
iy = int(y0 / dy)

def vx(x):
    return 0.02*(x +20)+0.01 

def vy(y):
    return 0.02*(y-20) 

x_coords = np.arange(Nx) * dx
y_coords = np.arange(Ny) * dy

vx_field = vx(x_coords)[:, None]
vy_field = vy(y_coords)[None, :]

K = alpha*dt/(dx*dy)
vx2 = vx_field*dt/dx
vy2 = vy_field*dt/dy

vx_pos = np.maximum(vx2, 0) 
vx_neg = np.minimum(vx2, 0)
vy_pos = np.maximum(vy2, 0)
vy_neg = np.minimum(vy2, 0)


### Initialisation

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
            
    T[ix, iy] += Q * dt
    
    for x in range(1,Nx-1):
        T[x,0] = (Told[x-1,0]+Told[x+1,0]+ 2*Told[x,1] - 2*dx*bc[2])/4 
        T[x,-1] = (Told[x-1,-1]+Told[x+1,-1]+ 2*Told[x,-2] - 2*dx*bc[0])/4
        
    for y in range(1,Ny-1):
        T[0,y] = (Told[0,y+1]+Told[0,y-1]+ 2*Told[1,y] - 2*dx*bc[3])/4 
        T[-1,y] = (Told[-1,y+1]+Told[-1,y-1]+ 2*Told[-2,y] - 2*dx*bc[1])/4 
        
    T[0,0] = (Told[1,0]+Told[0,1])/2
    T[-1,0] = (Told[-1,1]+Told[-2,0])/2
    T[0,-1] = (Told[0,-2]+Told[1,-1])/2
    T[-1,-1] = (Told[-1,-2]+Told[-2,-1])/2

    Told=T.copy()

    if k%100==0:
        print(k*dt)
        #print(T)

        
###Plotting the result

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

###     Echelle de couleur logarithmique
"""
#Creating vexctors with the x and y values of the grid and plotting the result
x=np.linspace(0,Lx,Nx)
y=np.linspace(0,Ly,Ny)
Y, X = np.meshgrid(y, x)

fig, ax=plt.subplots() 
surf = ax.contourf(X, Y, T,
                   levels=100,
                   norm=LogNorm(vmin=T.min()+1e-6, vmax=T.max()),
                   cmap=cm.coolwarm)
bar=fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.show()
"""
