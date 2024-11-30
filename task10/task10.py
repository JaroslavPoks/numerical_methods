import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


values1 = np.loadtxt("1.txt")
N_1 = len(values1)
M_1 = len(values1[0])

values2 = np.loadtxt("2.txt")
N_2 = len(values2)
M_2 = len(values2[0])

x_1 = np.linspace(0, 1, M_1)
x_2 = np.linspace(0, 1, M_2)

y_1 = np.linspace(0, 1, N_1)
y_2 = np.linspace(0, 1, N_2)

xgrid_1, ygrid_1 = np.meshgrid(x_1, y_1)
xgrid_2, ygrid_2 = np.meshgrid(x_2, y_2)

fig = plt.figure(figsize=(14, 7))
ax_3d = fig.add_subplot(1, 2, 1, projection='3d')
ax_3d.set(xlabel="$x$", ylabel="$y$", zlabel="$Explicit$")
ax_3d.plot_surface(xgrid_1, ygrid_1, values1, rstride=1, cstride=5, cmap='plasma')
ax_3d = fig.add_subplot(1, 2, 2, projection='3d')
ax_3d.plot_surface(xgrid_2, ygrid_2, values2, rstride=1, cstride=5, cmap='plasma')
ax_3d.set(xlabel="$x$", ylabel="$y$", zlabel="$Implicit$")
plt.savefig('1.png')
plt.show()

