import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from math import sqrt

### Constant

Lx=60 #m
Ly=60 #m
#Nx=200
#Ny=200
dx=1
dy=1
Nx = round(Lx/dx)#nb d'itérations selon x
Ny = round(Ly/dy)# // selon y

dt=1
tfinal=10000
Nt=int(tfinal//dt)#// dans le temps

alpha=5e-3
T0=10#degres C
Q=0.2
bc=[0,0,0,0] # border (right, up, left, down)

borehole=[30,30]

### Initialisation


Told=np.zeros((Nx,Ny))+T0
T=np.zeros((Nx,Ny))

result=np.zeros((Nx,Ny,3))

 
#Calculating the temperature profile timestep by timestep
for k in range(1,Nt+1): #pour chaque itérations de temps (100 000)
    
    for x in range(1,Nx-1):
        T[x,0] = (Told[x-1,0]+Told[x+1,0]+ 2*Told[x,1] - 2*dx*bc[2])/4 #ok
        T[x,-1] = (Told[x-1,-1]+Told[x+1,-1]+ 2*Told[x,-2] - 2*dx*bc[0])/4 #ok
        
    for y in range(1,Ny-1):
        T[0,y] = (Told[0,y+1]+Told[0,y-1]+ 2*Told[1,y] - 2*dx*bc[3])/4 #ok
        T[-1,y] = (Told[-1,y+1]+Told[-1,y-1]+ 2*Told[-2,y] - 2*dx*bc[1])/4 #ok
        
    T[0,0] = (Told[1,0]+Told[0,1])/2 #ok
    T[-1,0] = (Told[-1,1]+Told[-2,0])/2 #ok
    T[0,-1] = (Told[0,-2]+Told[1,-1])/2 #ok
    T[-1,-1] = (Told[-1,-2]+Told[-2,-1])/2 #ok

    for x in range(1,Nx-1):
        for y in range(1,Ny-1):
            
            K =T[x+1,y]+T[x-1,y]+T[x,y+1]+T[x,y-1] - 4*T[x,y] #ok
            T[x,y]= Told[x,y] + alpha*dt*K/(dx*dy) 
            
            if abs(x-borehole[0])<dx and abs(y-borehole[1])<dy:
                T[x,y] += Q*dt

    Told=T.copy()
    if round(k*dt)==2000:
        result[:,:,0]=T
    if round(k*dt)==5000:
        result[:,:,1]=T
    if round(k*dt)== 10000:
        result[:,:,2]=T
    
    if k%100==0:
        print(k)
        
###Plotting the result
##

fig=plt.figure()

    
#Creating vexctors with the x and y values of the grid and plotting the result
x=np.linspace(0,Lx,Nx)
y=np.linspace(0,Ly,Ny)
Y, X = np.meshgrid(y, x)


vmin = np.min(result[:, :, :3])
vmax = np.max(result[:, :, :3])
levels = np.linspace(vmin, vmax, 21)

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

t=[2000,5000,10000]

for i in range(3):
    ax = axes[i]
    surf = ax.contourf(X, Y, result[:, :, i], levels=levels, cmap=cm.coolwarm)
    
    # cercle
    theta = np.linspace(0, 2*np.pi, 400)
    L= sqrt(4*alpha*t[i])
    xc = borehole[0] + L * np.cos(theta)
    yc = borehole[1] + L * np.sin(theta)
    ax.plot(xc, yc, 'k-', linewidth=2)
    
    ax.set_title("t = " + str(t[i]) + " s")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

cbar = fig.colorbar(surf, ax=axes.ravel().tolist(),
                    orientation='horizontal', pad=0.15, shrink=0.9)


plt.tight_layout()
plt.show()
