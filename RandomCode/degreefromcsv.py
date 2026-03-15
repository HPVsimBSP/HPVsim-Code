import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

all_deg = pd.read_csv("degree_distribution.csv")
unique_deg, counts = np.unique(all_deg, return_counts=True)
np.append(0, unique_deg)
props = counts / sum(counts)
cum_prop = []
dummy = 0
for prop in props:
    dummy += prop
    cum_prop.append(1 - dummy)
fig, ax = plt.subplots()
ax.plot(np.arange(0, len(unique_deg)), cum_prop)
ax.set_xscale("log")
ax.set_yscale("log")
plt.ylim(top = 1)
plt.title('log-log Degree Distribution- Random Network')
plt.xlabel('Number of Partners (log scale)')
plt.ylabel('proportion')
plt.xlim(right = len(unique_deg))
#ax.fill_between(np.arange(1, 58), 0, cum_prop, color = 'tab:blue', alpha =0.31)
plt.show()