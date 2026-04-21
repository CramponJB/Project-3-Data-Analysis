import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

### Constant

Lx=60
Ly=60
#Nx=200
#Ny=200
dx=0.5
dy=0.5
Nx = round(Lx/dx)
Ny = round(Ly/dy)

dt=0.1
tfinal=1500
Nt=int(tfinal//dt)

T0=10
Q=1
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

### Initialisation


Told=np.zeros((Nx,Ny))+T0
T=np.zeros((Nx,Ny))
result=np.zeros((Nx,Ny,3))
 
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
    if abs(k*dt - 200) < dt:
        result[:,:,0]=T
        #print("T200=",T)
    if abs(k*dt - 500) < dt:
        result[:,:,1]=T
        #print("T500=",T)
    if abs(k*dt - 1000) < dt:
        result[:,:,2]=T
        #print("T1000=",T)

        
###Plotting the result
##

fig=plt.figure()

    
#Creating vexctors with the x and y values of the grid and plotting the result
x=np.linspace(0,Lx,Nx)
y=np.linspace(0,Ly,Ny)
Y, X = np.meshgrid(y, x)

T_levels=np.linspace(0,np.amax(T),21)
fig, ax=plt.subplots() 
surf=ax.contourf(X,Y,T,levels=T_levels,cmap=cm.coolwarm)
bar=fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.show()

fig=plt.figure()

    
#Creating vexctors with the x and y values of the grid and plotting the result
x=np.linspace(0,Lx,Nx)
y=np.linspace(0,Ly,Ny)
Y, X = np.meshgrid(y, x)


vmin = np.min(result[:, :, :3])
vmax = np.max(result[:, :, :3])
levels = np.linspace(vmin, vmax, 21)

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

t=[200,500,1000]

for i in range(3):
    ax = axes[i]
    surf = ax.contourf(X, Y, result[:, :, i], levels=levels, cmap=cm.coolwarm)
    

    
    ax.set_title("t = " + str(t[i]) + " s")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

cbar = fig.colorbar(surf, ax=axes.ravel().tolist(),
                    orientation='horizontal', pad=0.15, shrink=0.9)


plt.tight_layout()
plt.show()
