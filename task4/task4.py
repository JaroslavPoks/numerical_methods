import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math


def func(x):
    #return 0.15 * x * x * x * x + 2 * x * x + 1
    #return np.sin(x) * x - x**2 * np.cos(x)
    #return 1 / (1 + x * x)
    #return x * x * x
    return np.sin(x)


def d_func(x):
    return np.cos(x)
    #return 3 * x * x


f1 = open("1.txt", 'r')
f2 = open("2.csv", 'r')

line1 = np.array(f1.readline().split(' '))
a, b, Z, M = float(line1[0]), float(line1[1]), int(line1[2]), int(line1[3])
values = np.array(re.split("\n", f1.read().rstrip()))
values = np.array(list(map(float, values)))

x = np.linspace(a, b, Z)
nods2 = np.linspace(a, b, M)
nods1 = np.linspace(a, b, 2 * M - 1)
y1 = np.array(list(map(d_func, x)))

runge = values[2 * M::2]
values_nods1 = values[:2 * M - 1]
values_nods2 = values[2 * M - 1::2]
values_nods3 = np.add(values_nods2, runge)

fig, ax = plt.subplots(figsize=(10, 7))
ax.grid()
ax.plot(x, y1, label='Производная')
#ax.plot(x, values, label='Интерполянт')
ax.scatter(nods1, values_nods1, label='Узлы h/2', color='r')
ax.scatter(nods2, values_nods2, label='Узлы h', color='b')
ax.scatter(nods2, values_nods3, label='Рунге', color='g')
ax.legend()
plt.xlim([x[0], round(x[-1] + 0.5)])
#plt.ylim([round(min(y) - 0.51), round(max(y) + 0.51)])
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.savefig('1.png')
#plt.show()


L_inf_r = max(abs(runge))
L_1_r = sum(abs(runge))
L_2_r = math.sqrt(sum(runge**2))
pd.set_option('display.float_format', '{:.6E}'.format)
df = pd.read_csv("2.csv", delimiter=',',
                 names=['h', 'h/2', 'Runge'])
row = {0: 'L_1 relative', 1: 'L_2 relative', 2: 'L_inf relative',
       3: 'L_1 absolute', 4: 'L_2 absolute', 5: 'L_inf absolute'}
df = df.rename(index=row)
df_r = pd.DataFrame([[L_1_r, L_2_r, L_inf_r],
                    [df.at['L_1 absolute', 'h'], df.at['L_2 absolute', 'h'], df.at['L_inf absolute', 'h']]],
                    columns=['L_1', 'L_2', 'L_inf'],
                    index=['Нормы главного члена погрешности по Рунге', 'Нормы на сетке с шагом h'])
print(df, '\n')
print(df_r)

f1.close()
f2.close()
