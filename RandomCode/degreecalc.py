import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from pathlib import Path

all_deg = []
edges = pd.read_csv(r"NetworkCSVs\combined_edges.csv", usecols=["source", "target", "seed"])
nodes = pd.read_csv(r"NetworkCSVs\combined_nodes.csv", usecols=["node", "bipartite", "seed"])
seeds = [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 24, 25, 26, 27, 28, 30, 31, 32, 33, 34, 36, 37, 38, 39, 40, 42, 43, 44, 45, 46, 48, 49, 50, 51, 52, 54, 55, 56, 57, 58] #The original seeds misses every 6th seed due to multsim
#builds a graph for each seed
def build_graph(seed: int) -> nx.Graph:
    nodes_s = nodes[nodes["seed"] == seed][["node", "bipartite"]]
    edges_s = edges[edges["seed"] == seed][["source", "target"]]

    G = nx.Graph()
    G.add_nodes_from((r.node, {"bipartite": int(r.bipartite)})
                     for r in nodes_s.itertuples(index=False))
    G.add_edges_from(edges_s.itertuples(index=False, name=None))
    return G

for seed_i in seeds:
    G = build_graph(seed_i)
    print(f'Seed {seed_i} done')
    # Degree added
    degrees = np.array([d for _, d in G.degree()], dtype=int)
    all_deg.append(degrees)

all_deg = np.concat(all_deg)
np.savetxt("degree_distribution.csv", all_deg, delimiter=",")
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