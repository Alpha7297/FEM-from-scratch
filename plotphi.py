import numpy as np
import matplotlib.pyplot as plt
title=input("请输入图表名")
points=[]
with open("points.txt","r") as f:
    for line in f:
        x,y=map(float,line.split())
        points.append([x,y])
points=np.array(points)
triangles=[]
with open("triangles.txt","r") as f:
    for line in f:
        a,b,c=map(int,line.split())
        triangles.append([a,b,c])
triangles=np.array(triangles)
phi=[]
with open("phi.txt","r") as f:
    for line in f:
        phi.append(float(line))
phi=np.array(phi)

plt.figure(figsize=(8,8))
plt.tripcolor(points[:,0],points[:,1],triangles,phi,shading='gouraud',cmap='viridis')
plt.colorbar(label='phi')
plt.triplot(points[:,0],points[:,1],triangles,color='k',linewidth=0.5,alpha=0.5)
plt.xlabel('x')
plt.ylabel('y')
plt.title(title)
plt.axis('equal')
plt.grid(False)
plt.show()