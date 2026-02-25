import matplotlib.pyplot as plt
import numpy as np
points=np.loadtxt("points.txt")
px=np.array(points[:,0])
py=np.array(points[:,1])
np.append(px,px[0])
np.append(py,py[0])
plt.plot(px,py,'o-',color='blue')
plt.axis('equal')
plt.show()