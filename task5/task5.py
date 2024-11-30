import pandas as pd
import  numpy as np

f = open('1.csv', 'r')

K = int(f.readline())
row = {1: 'Rectangle', 2: 'Trapezoid', 3: 'Simpson', 4: 'Newton-Cotes', 5: 'Gauss'}

df = pd.read_csv("1.csv", delimiter=',',
                 names=["K = {}".format(K), "K = {}".format(2 * K)]).drop([0], axis=0)
df = df.rename(index=row)
print(df)
f.close()
