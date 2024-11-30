import numpy as np
import matplotlib.pyplot as plt


values1 = np.loadtxt("1.txt")
x1 = values1[:, 0]
y1 = values1[:, 1]
z1 = values1[:, 2]

values2 = np.loadtxt("2.txt")
x2 = values2[:, 0]
y2 = values2[:, 1]
z2 = values2[:, 2]

values3 = np.loadtxt("3.txt")
x3 = values3[:, 0]
y3 = values3[:, 1]
z3 = values3[:, 2]

values4 = np.loadtxt("4.txt")
x4 = values4[:, 0]
y4 = values4[:, 1]
z4 = values4[:, 2]

values5 = np.loadtxt("5.txt")
x5 = values5[:, 0]
y5 = values5[:, 1]
z5 = values5[:, 2]

fig, ax = plt.subplots(ncols=2, figsize=(12, 5))

ax[0].grid()
ax[0].plot(x1, y1, label='y - метод Эйлера')
ax[0].plot(x2, y2, label='y - метод Эйлера с пересчетом')
ax[0].plot(x3, y3, label='y - метод Рунге-Кутта 2 порядка')
ax[0].plot(x4, y4, label='y - метод Рунге-Кутта 4 порядка')
ax[0].plot(x5, y5, label='y - метод Адамса')
ax[0].legend()

ax[1].grid()
ax[1].plot(x1, z1, label="y' - метод Эйлера")
ax[1].plot(x2, z2, label="y' - метод Эйлера с пересчетом")
ax[1].plot(x3, z3, label="y' - метод Рунге-Кутта 2 порядка")
ax[1].plot(x4, z4, label="y' - метод Рунге-Кутта 4 порядка")
ax[1].plot(x5, z5, label="y' - метод Адамса")
ax[1].legend()

plt.savefig('1.png')
plt.show()

