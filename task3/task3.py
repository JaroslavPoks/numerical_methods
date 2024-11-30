import re
import matplotlib.pyplot as plt
import numpy as np
import math


def func(x):
    #return 0.15 * x * x * x * x + 2 * x * x + 1
    #return np.sin(x) * x - x**2 * np.cos(x)
    #return x**3
    return np.sin(x) * x


f = open("1.txt", 'r')

line1 = list(f.readline().split(' '))
a = float(line1[0])
b = float(line1[1])
Z = int(line1[2])
KL = int(line1[3])
tmp = list(re.split("\n", f.read().rstrip()))
tmp = list(map(float, tmp))

nods = tmp[:KL]
values = tmp[KL:]

x = np.linspace(a, b, Z)
y1 = list(map(func, x))
values_nods = list(map(func, nods))

fig, ax = plt.subplots(figsize=(10, 7))
ax.grid()
ax.plot(x, y1, label='Функция')
ax.plot(x, values, label='Интерполянт')
ax.scatter(nods, values_nods, label='Узлы', color='r')
ax.legend()
plt.xlim([x[0], round(x[-1] + 0.5)])
#plt.ylim([round(min(y) - 0.51), round(max(y) + 0.51)])
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.savefig('1.png')
plt.show()

f.close()
