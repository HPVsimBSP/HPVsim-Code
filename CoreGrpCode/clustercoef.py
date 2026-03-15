import pandas as pd
import numpy as np
import networkx as nx

edges = pd.read_csv("NetworkCSVs\combined_edges.csv", usecols=["source", "target", "seed"])
nodes = pd.read_csv("NetworkCSVs\combined_nodes.csv", usecols=["node", "bipartite", "seed"])

#builds a graph for each seed like in degree calc
def build_graph(seed: int) -> nx.Graph:
    nodes_s = nodes[nodes["seed"] == seed][["node", "bipartite"]]
    edges_s = edges[edges["seed"] == seed][["source", "target"]]

    G = nx.Graph()
    G.add_nodes_from((r.node, {"bipartite": int(r.bipartite)})
                     for r in nodes_s.itertuples(index=False))
    G.add_edges_from(edges_s.itertuples(index=False, name=None))
    return G

all_rows = []

for seed in sorted(edges["seed"].unique()):
    G = build_graph(int(seed))

    # 4-cycle clustering coefficient C4
    c4 = nx.square_clustering(G)  # dict: node -> C4

    bip = nx.get_node_attributes(G, "bipartite")
    df = pd.DataFrame({
        "node": list(c4.keys()),
        "seed": int(seed),
        "bipartite": [bip.get(n, np.nan) for n in c4.keys()],
        "C4": list(c4.values()),
    })
    all_rows.append(df)
    print(f"seed {seed} done")

out = pd.concat(all_rows, ignore_index=True)
out.to_csv("c4_clustering_by_node.csv", index=False)

# Example: overall mean C4
print(out["C4"].mean())