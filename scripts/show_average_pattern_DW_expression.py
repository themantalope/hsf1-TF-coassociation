import pandas as pd
import os
import matplotlib.pyplot as plt
import maakclusterutils
import pydendroheatmap as pdh
import re
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler, RobustScaler

data_dir = os.path.join("..", "data_files")
dw_expression_file = os.path.join(data_dir, "hsf1_DW_reg_oscillatory_genes.xlsx")

dw_expression = pd.read_excel(open(dw_expression_file))

dw_expression = dw_expression[["ENSG", "Symbol"] + [c for c in dw_expression.columns if "FPKM" in c and "Normalized" not in c]]

expr_data = dw_expression[[c for c in dw_expression.columns if "FPKM" in c]]
# expr_data["mean"] = expr_data.apply(lambda x: x.mean(), axis=1)
expr_data.sort_values("Avg FPKM", axis=0, ascending=False, inplace=True)
print expr_data

# expr_data = expr_data[expr_data["mean"] < 300]
expr_data = expr_data[expr_data["Avg FPKM"] > 25]
expr_data = expr_data[[c for c in expr_data.columns if "FPKM" in c and "Avg" not in c]]

rolling_df = pd.DataFrame(columns=expr_data.iloc[0,:].index.values)

for i in range(0, len(expr_data)):
    rolling_df.loc[expr_data.index[i],expr_data.iloc[i,:].index.values] = expr_data.iloc[i,:].rolling(window=2).mean().values


print rolling_df
rolling_df.dropna(axis=1,how="all", inplace=True)
print rolling_df
mean_expr = expr_data.mean()
std_expr = expr_data.std()
print mean_expr.values.tolist()
print len(mean_expr)


mmsc = MinMaxScaler()
rscl = RobustScaler()

# expr_norm = rscl.fit_transform(expr_data.values.T)
# expr_norm = expr_norm.T
expr_data = rolling_df
# print type(expr_data)
# print type(expr_data.values)
# raw_input()
# expr_data.values.astype(float).std(axis=1)
expr_norm = maakclusterutils.normalize_matrix_along_axis(expr_data.values.astype(float))

# expr_norm = maakclusterutils.max_scale_matrix_along_axis(expr_data.astype(float).values)

# expr_norm = mmsc.fit_transform(expr_data.astype(float))

print expr_norm[0,:]

#          1      2     3       4   5       6       7     8     9   10     11   12    13    14      15   16
xticks = ["G1/S","S","S/G2", "G2", "G2/M", "G1", "G1", "G1/S", "S","S/G2","G2","G2/M","M","M/G1", "G1"]

plt.plot(expr_norm.mean(axis=0))
plt.xticks(range(len(mean_expr)),xticks)
plt.xlabel("Cell Cycle Phase")
plt.ylabel("Mean FPKM")
# plt.show()



link, origidx, clustered = maakclusterutils.UPGMACluster(maakclusterutils.calculateDistanceMatrix(expr_norm), get_original_indexing=True,get_reorderd_matrix=True)
# link, origidx, clustered = maakclusterutils.cluster_correlation_matrix(corrmat, get_reordered_matrix=True, get_original_indexing=True)

# print origidx

plt.matshow(clustered, interpolation="nearest")
# plt.show()



plt.matshow(expr_norm, interpolation="nearest")
# plt.show()

inertia = []
for k in range(2,10):
    km = KMeans(n_clusters=k)
    km.fit(expr_norm)
    inertia.append(km.inertia_)

plt.plot(range(2,10),inertia)
# plt.show()


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
# plt.show()

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
roll = maakclusterutils.normalize_matrix_along_axis(rolling_df.astype(float).values[rowreindex, :])
print type(roll)

hm = pdh.DendroHeatMap(roll,
                       left_dendrogram=rowlink,
                       col_labels=[str(x+1) for x in range(expr_norm.shape[0])] )
# hm.show()

g1_norm_expr = maakclusterutils.normalize_matrix_along_axis(dw_g1[[c for c in dw_g1.columns if "FPKM" in c and "Normalized" not in c and "Avg" not in c]].values)

g1rowlink, g1rowreindex = maakclusterutils.UPGMACluster(maakclusterutils.calculateDistanceMatrix(g1_norm_expr), get_original_indexing=True)
g1_norm_expr =g1_norm_expr[g1rowreindex, :]

hm.heat_map_data = g1_norm_expr
hm.left_dendrogram = g1rowlink
# hm.show()

g2_norm_expr = maakclusterutils.normalize_matrix_along_axis(dw_g2[[c for c in dw_g1.columns if "FPKM" in c and "Normalized" not in c and "Avg" not in c]].values)

g2rowlink, g2rowreindex = maakclusterutils.UPGMACluster(maakclusterutils.calculateDistanceMatrix(g2_norm_expr), get_original_indexing=True)

g2_norm_expr = g2_norm_expr[g2rowreindex, :]


hm.heat_map_data = g2_norm_expr
hm.left_dendrogram = g2rowlink
# hm.show()

e2f1_expr = dw_expression[dw_expression["Symbol"].isin(e2f1_syms)]




bmyb_syms = cctf_pos["BMYB"]["Symbol"].values
bmyb_expr = dw_expression[dw_expression["Symbol"].isin(bmyb_syms)]

overlap = list(set(e2f1_syms).intersection(set(bmyb_syms)))
print "e2f1 bmyb overlap: ", len(overlap)

e2f1_expr = e2f1_expr[~e2f1_expr["Symbol"].isin(overlap)]
bmyb_expr = bmyb_expr[~bmyb_expr["Symbol"].isin(overlap)]

print "e2f1 unique: ", len(e2f1_expr)
print "bmyb unique: ", len(bmyb_expr)

for i in range(0, len(e2f1_expr)):
    row = e2f1_expr.loc[e2f1_expr.index[i], [c for c in e2f1_expr.columns if "FPKM" in c and "Normalized" not in c and "Avg" not in c]]
    row = row.rolling(window=2).mean().values
    e2f1_expr.loc[e2f1_expr.index[i], [c for c in e2f1_expr.columns if "FPKM" in c and "Normalized" not in c and "Avg" not in c]] = row

e2f1_expr.dropna(inplace=True,axis=1)

e2f1_norm = maakclusterutils.normalize_matrix_along_axis(e2f1_expr[[c for c in e2f1_expr.columns if "FPKM" in c and "Normalized" not in c and "Avg" not in c]].values)
e2f1_link , e2f1_reindex = maakclusterutils.UPGMACluster(maakclusterutils.calculateDistanceMatrix(e2f1_norm), get_original_indexing=True)


e2f1_argmax = e2f1_expr[[c for c in e2f1_expr.columns if "FPKM" in c and "Normalized" not in c and "Avg" not in c]].idxmax(axis=1).values
print e2f1_argmax

e2f1_cols = [c for c in e2f1_expr.columns if "FPKM" in c and "Normalized" not in c and "Avg" not in c]

e2f1_phase_max = [xticks[e2f1_cols.index(c)+1] for c in e2f1_argmax]
g1s_patt = re.compile("(?P<g1s>G1|S)")
g2m_patt = re.compile("(?P<g2m>G2|M)")

e2f1_n_g1s = sum([1 for x in e2f1_phase_max if bool(re.match(g1s_patt, x))])
e2f1_n_g2m = sum([1 for x in e2f1_phase_max if bool(re.match(g2m_patt, x))])

print e2f1_n_g1s, "e2f1 n g1s"
print e2f1_n_g2m, "e2f1 n g2m"

bmyb_argmax = bmyb_expr[e2f1_cols].idxmax(axis=1).values
bmyb_phase_max = [xticks[e2f1_cols.index(c)+1] for c in bmyb_argmax]

bmyb_n_g1s = sum([1 for x in bmyb_phase_max if bool(re.match(g1s_patt, x))])
bmyb_n_g2m = sum([1 for x in bmyb_phase_max if bool(re.match(g2m_patt, x))])

print bmyb_n_g1s, "bmyb n g1s"
print bmyb_n_g2m, "bmyb n g2m"

hm.heat_map_data = e2f1_norm[e2f1_reindex, :]
hm.row_labels = e2f1_expr["Symbol"].values[e2f1_reindex].tolist()
hm.title = "HSF1 E2F1 Expr"
hm.left_dendrogram = e2f1_link
hm.color_legend_title = "Normalized FPKM"
hm.export(filename=os.path.join(data_dir, "cctf_hsf1_expression_figures", "hsf1_e2f1_rnaseq_hm.pdf"))


for i in range(0, len(bmyb_expr)):
    row = bmyb_expr.loc[bmyb_expr.index[i], [c for c in dw_g1.columns if "FPKM" in c and "Normalized" not in c and "Avg" not in c]]
    row = row.rolling(window=2).mean().values
    bmyb_expr.loc[bmyb_expr.index[i], [c for c in dw_g1.columns if "FPKM" in c and "Normalized" not in c and "Avg" not in c]] = row


bmyb_expr.dropna(axis=1, inplace=True)

bmyb_norm = maakclusterutils.normalize_matrix_along_axis(bmyb_expr[[c for c in bmyb_expr.columns if "FPKM" in c and "Normalized" not in c and "Avg" not in c]].values)
bmyb_link, bmyb_reindex = maakclusterutils.UPGMACluster(maakclusterutils.calculateDistanceMatrix(bmyb_norm), get_original_indexing=True)

# plt.plot(bmyb_norm.mean(axis=0)); plt.title("BMYB Mean");plt.show()



hm.heat_map_data = bmyb_norm[bmyb_reindex, :]
hm.title = 'HSF1 BMYB'
hm.row_labels = bmyb_expr["Symbol"].values[bmyb_reindex].tolist()
hm.left_dendrogram = bmyb_link
hm.color_legend_title = "Normalized FPKM"
hm.export(filename=os.path.join(data_dir, "cctf_hsf1_expression_figures", "hsf1_bmyb_rnaseq_hm.pdf"))
