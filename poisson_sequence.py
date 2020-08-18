# -*- coding: utf-8 -*-
"""
Spyder Editor

Dies ist eine temporäre Skriptdatei.
"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set(context = "poster", style = "white")
sns.despine()
a = np.random.exponential(size = 10)
b = np.cumsum(a)

c = np.zeros_like(b)

fig, ax = plt.subplots()
ax.scatter(b, c, c = "k")
ax.axhline(c = "k")

fig.savefig("poisson_sequence.svg")