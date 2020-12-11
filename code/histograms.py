import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(context = "poster", style = "ticks")
n = 100
y1 = np.random.uniform(size=n)
print(sum(y1))

y2 = np.copy(y1)*np.random.normal(loc = 1, scale = 0.1, size =n)
y2[y2>1] = 1

low = 30
hi  = 60
y2[low:hi] = y2[low:hi]+(1-y2[low:hi])/2

colors = ["tab:grey"]*n
for i in range(low,hi):
    colors[i] = "tab:blue"


x = np.arange(1,n+1)
fig, (ax1, ax2) = plt.subplots(1,2, figsize = (10,2.0))
ax1.bar(x, width = 0.2, height = y1, color = colors, edgecolor = colors)
ax2.bar(x, width = 0.2, height = y2, color = colors, edgecolor = colors)

for ax in (ax1, ax2):
    ax.set_ylim(0,1.2)
    ax.set_xlabel("genes")

    ax.set_xticks([])
ax1.set_ylabel("")
ax2.set_ylabel("")
plt.show()

fig.savefig("histograms.svg")