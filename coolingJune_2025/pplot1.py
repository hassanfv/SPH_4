
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv("out.csv")

x = df["x"]
y = df["y"]

#plt.figure(figsize = (16, 12))

plt.plot(x, y, c = 'k')
plt.scatter(x, y, c = 'r', s = 30)
plt.show()

