import re
import matplotlib.pyplot as plt
import numpy as np
import math


def func(x):
    #return 0.15 *  x * x * x * x + 2 * x * x + 1
    return np.sin(x)


f = open("1.txt", 'r')

line1 = list(f.readline().split(' '))
a = float(line1[0])
b = float(line1[1])
N = int(line1[2])
M = int(line1[3])
values_nods = list(re.split(" ", f.readline().rstrip()))
values_nods = list(map(float, values_nods))
nods = list(re.split(" ", f.readline().rstrip()))
nods = list(map(float, nods))
values = list(re.split(" ", f.readline().rstrip()))
values = list(map(float, values))
# print(values_nods, "\n", nods)

x = np.linspace(a, b, M)
y = list(map(func, x))

fig, ax = plt.subplots(figsize=(10, 7))
ax.grid()
ax.plot(x, y, label='Функция')
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
