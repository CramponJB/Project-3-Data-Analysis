import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from math import sqrt

### Constant

Lx=60 # m
Ly=60 # m
dx=dy=0.3 # m
Nx = round(Lx/dx) # 200m (nb d'itérations selon x)
Ny = round(Ly/dy) # 200m (// selon y)
dt=0.1
tfinal=10000
Nt=int(tfinal//dt) # 100 000m (// dans le temps)

alpha=5e-3 # m²/s
T0=10 # °C
Q=0.2 # °C/s

borehole=[30,30] # m
ix = int(borehole[0] / dx) # 100m
iy = int(borehole[1] / dy) # 100m

bc=[0,0,0,0] # border (right, up, left, down)

### Verification

r = alpha * dt / dx**2
print(f"Nombre de diffusion r = {r:.4f}  (doit être ≤ 0.25 en 2D)")
assert r <= 0.25, "Schéma instable !"

### Initialisation

Told=np.zeros((Nx,Ny))+T0
T=np.zeros((Nx,Ny))

result=np.zeros((Nx,Ny,3))

#Calculating the temperature profile timestep by timestep
for k in range(1,Nt+1): #pour chaque itérations de temps (100 000)
    
    K = (
        Told[2:,1:-1] +
        Told[:-2,1:-1] +
        Told[1:-1,2:] +
        Told[1:-1,:-2] -
        4*Told[1:-1,1:-1]
    )
    
    T[1:-1,1:-1] = Told[1:-1,1:-1] + alpha*dt*K/(dx**2)
    
    T[ix, iy] += Q * dt  
    
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
    
    Told=T.copy()
        
    if round(k*dt)==2000:
        result[:,:,0]=T
    if round(k*dt)==5000:
        result[:,:,1]=T
    if round(k*dt)== 10000:
        result[:,:,2]=T
    if k%10000==0:
        print("t =", k*dt, "s")
    
        
###Plume radius
eps = 0.05  # seuil de température

X, Y = np.meshgrid(np.linspace(0, Lx, Nx),
                   np.linspace(0, Ly, Ny))

dist = np.sqrt((X - borehole[0])**2 + (Y - borehole[1])**2) # distance au borehole
mask = (T - T0) > eps # zones "chaudes"

if np.any(mask):
    plume_radius = np.max(dist[mask]) # plume radius = distance max dans la zone chaude
else:
    plume_radius = 0

print("Plume radius:", plume_radius)
    

###Plotting the result

fig=plt.figure()

#Creating vexctors with the x and y values of the grid and plotting the result
x=np.linspace(0,Lx,Nx)
y=np.linspace(0,Ly,Ny)
Y, X = np.meshgrid(y, x)

vmin = np.min(result[:, :, :3])
vmax = np.max(result[:, :, :3])
levels = np.linspace(vmin, vmax, 21)

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle("Task 1 – Diffusion", fontsize=13)

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

plt.show()
