import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

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
tfinal=3000
Nt=int(tfinal//dt)

alpha=5e-3
T0=10
Q=0.2
bc=[0,0,0,0] # border (right, up, left, down)

borehole=[30,30]

def vx(x):
    return 0.02*(x +20)+0.01 

def vy(y):
    return 0.02*(y-20) 

### Initialisation


Told=np.zeros((Nx,Ny))+T0
T=np.zeros((Nx,Ny))

 
#Calculating the temperature profile timestep by timestep
for k in range(1,Nt+1):
    
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

    for x in range(1,Nx-1):
        for y in range(1,Ny-1):
            
            if vx(x) > 0:
                dTdx = (Told[x,y] - Told[x-1,y]) / dx
            else:
                dTdx = (Told[x+1,y] - Told[x,y]) / dx
            

            if vy(y) > 0:
                dTdy = (Told[x,y] - Told[x,y-1]) / dy
            else:
                dTdy = (Told[x,y+1] - Told[x,y]) / dy
                
            T[x,y] = Told[x,y] - dt*( vx(x)*dTdx + vy(y)*dTdy )
            
           
            if abs(x-borehole[0])<dx and abs(y-borehole[1])<dy:
                T[x,y] += Q*dt

    Told=T.copy()

    if k%100==0:
        print(k)
        #print(T)

        
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