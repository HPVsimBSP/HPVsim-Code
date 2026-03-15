import pandas as pd
import numpy as np

nodes = pd.read_csv("c4_clustering_by_node.csv", usecols=["node", "C4"])

non_zero = nodes[nodes["C4"] != 0]["C4"] #store all non zero values of C4
print(f"Average Clustering Coeff (non-zero)= {non_zero.mean()}")
allmean = nodes["C4"].mean()
print(f"Average C4 across all= {allmean}")
prop = len(non_zero)/len(nodes)
print(f"prop of non-zero C4 nodes is: {prop}")