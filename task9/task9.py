import numpy as np
import matplotlib.pyplot as plt


values1 = np.loadtxt("1.txt")
x1 = values1[:, 0]
y1 = values1[:, 1]

values2 = np.loadtxt("2.txt")
x2 = values2[:, 0]
y2 = values2[:, 1]

values3 = np.loadtxt("3.txt")
x3 = values3[:, 0]
y3 = values3[:, 1]

values4 = np.loadtxt("4.txt")
x4 = values4[:, 0]
y4 = values4[:, 1]

fig, ax = plt.subplots(ncols=2, figsize=(12, 5))

ax[0].grid()
ax[0].plot(x1, y1, label='y - метод пристрелки')
ax[0].plot(x2, y2, label='y - метод конечных разностей')
ax[0].legend()

ax[1].grid()
ax[1].plot(x3, y3, label='y - метод Ньютона')
ax[1].plot(x4, y4, label='y - метод секущих')
ax[1].legend()

plt.savefig('1.png')
plt.show()

