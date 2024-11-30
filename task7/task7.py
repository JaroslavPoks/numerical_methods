import matplotlib.pyplot as plt
import numpy as np


def f1(x, y):
    return x**2 - y - 1


def f2(x, y):
    return y - np.arctan(x)


x = np.linspace(-5, 5, 200)
y = np.linspace(-5, 5, 200)
X, Y = np.meshgrid(x, y)
Z = f1(X, Y)**2 + f2(X, Y)**2
C = np.arange(0.5, 5, 0.5)

fig, ax = plt.subplots(ncols=2, figsize=(12, 5))

ax[0].contour(X, Y, Z, C)
ax[0].grid()
ax[1].plot(x, x**2 - 1)
ax[1].plot(x, np.arctan(x))
ax[1].grid()
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.savefig('1.png')
plt.show()
