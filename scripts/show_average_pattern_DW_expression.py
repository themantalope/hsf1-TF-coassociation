import pandas as pd
import os
import matplotlib.pyplot as plt
import maakclusterutils
import pydendroheatmap as pdh
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler, RobustScaler


data_dir = os.path.join("..", "data_files")
dw_expression_file = os.path.join(data_dir, "hsf1_DW_reg_oscillatory_genes.xlsx")

dw_expression = pd.read_excel(open(dw_expression_file))

dw_expression = dw_expression[["ENSG", "Symbol"] + [c for c in dw_expression.columns if "FPKM" in c and "Normalized" not in c]]


# print dw_expression

# print [c for c in dw_expression.columns if "FPKM" in c]


expr_data = dw_expression[[c for c in dw_expression.columns if "FPKM" in c]]
# expr_data["mean"] = expr_data.apply(lambda x: x.mean(), axis=1)
expr_data.sort_values("Avg FPKM", axis=0, ascending=False, inplace=True)
print expr_data

# expr_data = expr_data[expr_data["mean"] < 300]
expr_data = expr_data[expr_data["Avg FPKM"] > 10]
expr_data = expr_data[[c for c in expr_data.columns if "FPKM" in c and "Avg" not in c]]

mean_expr = expr_data.mean()
std_expr = expr_data.std()
print mean_expr.values.tolist()
print len(mean_expr)


mmsc = MinMaxScaler()
rscl = RobustScaler()

# expr_norm = rscl.fit_transform(expr_data.values.T)
# expr_norm = expr_norm.T

expr_norm = maakclusterutils.max_scale_matrix_along_axis(expr_data.values)

print expr_norm[0,:]

#          1      2     3       4   5       6       7     8   9   10     11   12    13    14      15   16
xticks = ["G1/S","S","S/G2", "G2", "G2/M", "G1", "G1", "G1/S", "S","S","S/G2","G2","G2/M","M", "G1"]

plt.plot(expr_norm.mean(axis=0))
# plt.plot((mean_expr+std_expr).values.tolist(), linestyle="-")
# plt.plot((mean_expr-std_expr).values.tolist(), linestyle="-")
plt.xticks(range(len(mean_expr)),xticks)
plt.xlabel("Cell Cycle Phase")
plt.ylabel("Mean FPKM")
plt.show()




print expr_norm

# corrmat, corrmatpv = maakclusterutils.calculatePearsonCorrelationMatrix(expr_norm, getpvalmat=True)
# corrmat, pvmat = maakclusterutils.calculateSpearmanRankCorrelationMatrix(expr_data.values, getpvalmat=True)
# print corrmat

link, origidx, clustered = maakclusterutils.UPGMACluster(maakclusterutils.calculateDistanceMatrix(expr_norm), get_original_indexing=True,get_reorderd_matrix=True)
# link, origidx, clustered = maakclusterutils.cluster_correlation_matrix(corrmat, get_reordered_matrix=True, get_original_indexing=True)

# print origidx

plt.matshow(clustered, interpolation="nearest")
plt.show()



plt.matshow(expr_norm, interpolation="nearest")
plt.show()

inertia = []
for k in range(2,10):
    km = KMeans(n_clusters=k)
    km.fit(expr_norm)
    inertia.append(km.inertia_)

plt.plot(range(2,10),inertia)
plt.show()


km = KMeans(n_clusters=2)
km.fit(expr_data)
print km.labels_

g1_idx = []
g2_idx = []

for i, l in enumerate(km.labels_):
    if l == 0:
        g1_idx.append(i)
    else:
        g2_idx.append(i)

g1_expr = expr_norm[g1_idx, :]
g2_expr = expr_norm[g2_idx, :]

g1_mean = g1_expr.mean(axis=0)
g2_mean = g2_expr.mean(axis=0)

plt.plot(g1_mean, color="r")
plt.plot(g2_mean, color="b")
plt.xticks(range(len(mean_expr)),xticks)
plt.show()

# dhm = pydendroheatmap.DendroHeatMap(expr_data.values)
# dhm.show()

dw_expr_cut = dw_expression[dw_expression.index.isin(expr_data.index)]

dw_g1 = dw_expr_cut.iloc[g1_idx, :]
dw_g2 = dw_expr_cut.iloc[g2_idx, :]

# print dw_g1
# print dw_g2

g1w = pd.ExcelWriter(os.path.join(data_dir, "dw_group_1.xlsx"))
dw_g1.to_excel(g1w)
g1w.save()

g2w = pd.ExcelWriter(os.path.join(data_dir, "dw_group_2.xlsx"))
dw_g2.to_excel(g2w)
g2w.save()

cctf_pos_file = os.path.join(data_dir, "hsf1_cctf_position_analysis.xlsx")
cctf_pos = pd.read_excel(open(cctf_pos_file), sheetname=None)
e2f1_pos = cctf_pos["E2F1"]
e2f1_syms = e2f1_pos["Symbol"].values.tolist()

g1_dw_syms = dw_g1["Symbol"].values.tolist()
g2_dw_syms = dw_g2["Symbol"].values.tolist()

g1overlap = set(e2f1_syms).intersection(g1_dw_syms)
print "overlap :", len(g1overlap)
print "len g1: ", len(dw_g1)

g2overlap = set(e2f1_syms).intersection(g2_dw_syms)
print "g2 overlap: ",len(g2overlap)
print "len g2: ", len(dw_g2)

print expr_norm.shape
print len(xticks)

rowlink, rowreindex = maakclusterutils.UPGMACluster(maakclusterutils.calculateDistanceMatrix(expr_norm),get_original_indexing=True)

expr_norm = expr_norm[rowreindex, :]

hm = pdh.DendroHeatMap(expr_norm,
                       left_dendrogram=rowlink,
                       col_labels=[str(x+1) for x in range(expr_norm.shape[0])] )
hm.show()

g1_norm_expr = maakclusterutils.normalize_matrix_along_axis(dw_g1[[c for c in dw_g1.columns if "FPKM" in c and "Normalized" not in c and "Avg" not in c]].values)

g1rowlink, g1rowreindex = maakclusterutils.UPGMACluster(maakclusterutils.calculateDistanceMatrix(g1_norm_expr), get_original_indexing=True)
g1_norm_expr =g1_norm_expr[g1rowreindex, :]

hm.heat_map_data = g1_norm_expr
hm.left_dendrogram = g1rowlink
hm.show()

g2_norm_expr = maakclusterutils.normalize_matrix_along_axis(dw_g2[[c for c in dw_g1.columns if "FPKM" in c and "Normalized" not in c and "Avg" not in c]].values)

g2rowlink, g2rowreindex = maakclusterutils.UPGMACluster(maakclusterutils.calculateDistanceMatrix(g2_norm_expr), get_original_indexing=True)

g2_norm_expr = g2_norm_expr[g2rowreindex, :]


hm.heat_map_data = g2_norm_expr
hm.left_dendrogram = g2rowlink
hm.show()
