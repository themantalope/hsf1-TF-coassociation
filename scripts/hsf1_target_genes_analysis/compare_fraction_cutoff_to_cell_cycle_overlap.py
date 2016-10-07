import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import json

data_dir = os.path.join("..", "..", "data_files")
mapping_file = os.path.join(data_dir,"hsf1_mendillo_genes_map.json")


mendillo_genes_binding_file = os.path.join(data_dir, "mendillo_hsf1_binding_genes.xlsx")
oscillations_file = os.path.join(data_dir, "dominguez_wang_periodic_expression_modified.xlsx")

oscillations = pd.read_excel(open(oscillations_file))
mendillo_genes = pd.read_excel(open(mendillo_genes_binding_file))
mapping = json.load(open(mapping_file))



cutoffs = np.linspace(0.05, 0.95,num=19)
proportions = np.zeros(cutoffs.shape)

for i, c in enumerate(cutoffs):
    cur_group = mendillo_genes[mendillo_genes["frac"] >= c]
    genes = [mapping[gid] for gid in cur_group.index if gid in mapping]
    cur_osc = oscillations[oscillations["ENSG"].isin(genes)]
    proportions[i] = float(len(cur_osc)) / float(len(cur_group))

fig = plt.figure()
plt.plot(cutoffs, proportions)
plt.title("Proportion of genes reported with cell-cycle oscillation vs stringency of cutoff")
plt.xlabel("stringency of cutoff")
plt.ylabel("proportion of genes with cell cycle oscillatory patterns")
fig.savefig("proportion_vs_cutoff_all.pdf")
plt.show()


oscillations = oscillations[~(oscillations["Phase "] == "None")]


cutoffs = np.linspace(0.05, 0.95,num=19)
proportions = np.zeros(cutoffs.shape)

for i, c in enumerate(cutoffs):
    cur_group = mendillo_genes[mendillo_genes["frac"] >= c]
    genes = [mapping[gid] for gid in cur_group.index if gid in mapping]
    cur_osc = oscillations[oscillations["ENSG"].isin(genes)]
    proportions[i] = float(len(cur_osc)) / float(len(cur_group))

fig = plt.figure()
plt.plot(cutoffs, proportions)
plt.title("Proportion of genes reported with cell-cycle oscillation vs stringency of cutoff")
plt.xlabel("stringency of cutoff")
plt.ylabel("proportion of genes with cell cycle oscillatory patterns")
fig.savefig("proportion_vs_cutoff_1182_reported_genes.pdf")
plt.show()