import os
import matplotlib.pyplot as plt
import shutil
import pandas as pd
import maakclusterutils as cu
import pydendroheatmap as pdh

data_dir = os.path.join("..", "data_files")
gw_expression_file = os.path.join(data_dir, "grant_whitfield_cell_cycle_expression.xlsx")
# gw_annotations_file = os.path.join(data_dir, "grant_whitfield_probe_annotations.xlsx")
cctf_dist_file = os.path.join(data_dir, "hsf1_cctf_position_analysis.xlsx")

print "Loading"
cctf_dist = pd.read_excel(open(cctf_dist_file), sheetname=None)
gw_expression = pd.read_excel(open(gw_expression_file), sheetname=None)
print "Data loaded"
good_cols = ['0.1', 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38]
time_labels = ["S"] * 3 + ["S/G2"] + ["G2"] * 4 + ["M"] * 1 + ["G1"] * 3 + ["S"] * 3 + ["G2"] * 3 + ["M"] * 2
orig = gw_expression["original"]

# for k in cctf_dist.keys():
#     e2f1 = cctf_dist[k]
#     syms = [s for s in e2f1["Symbol"].values.tolist() if isinstance(s, basestring)]
#
#     idxs = []
#     # print syms
#
#     for s in syms:
#         found = orig[~pd.isnull(orig["NAME"])]
#         found = found[found["NAME"].str.contains(s)]
#         idxs += found.index.values.tolist()
#
#     if len(idxs) == 0 : continue
#     gw_data = orig.loc[idxs]
#     # print gw_data.columns.tolist()
#
#     means= gw_data[good_cols]
#     print len(means)
#     if len(means) < 5: continue
#     try:
#         means.interpolate(method="cubic", axis=0, inplace=True)
#     except:
#         continue
#     means.dropna(axis=0, inplace=True)
#     norms = cu.normalize_matrix_along_axis(means.astype(float).values)
#     plt.plot(means.mean().tolist()); plt.show()
#     link, rowreindex = cu.UPGMACluster(cu.calculateDistanceMatrix(norms, axis=0), get_original_indexing=True)
#
#     norms = norms[rowreindex, :]
#
#     hm = pdh.DendroHeatMap(heat_map_data=norms)
#     hm.left_dendrogram = link
#     hm.title = k
#     hm.show()

# print norms








e2f1_s = cctf_dist["E2F1"]["Symbol"].values.tolist()
bmyb_s = cctf_dist["BMYB"]["Symbol"].values.tolist()
foxm1_s = cctf_dist["FOXM1"]["Symbol"].values.tolist()

print "overlap: ", len(set(e2f1_s).intersection(set(bmyb_s)))
print "e2f1 num, bmyb num: ", len(e2f1_s), len(bmyb_s)

fandb = list(set(foxm1_s + bmyb_s))

e2f1_vs_fb_overlap = len(set(e2f1_s).intersection(set(fandb)))
e2f1_only = set(e2f1_s).difference(set(fandb))
fb_only = set(fandb).difference(set(e2f1_s))

print len(fb_only)
print len(e2f1_only)

e2f1_exclusive_idx = []
fb_exclusive_idx = []

for s in e2f1_only:
    found = orig[~pd.isnull(orig["NAME"])]
    found = found[found["NAME"].str.contains(s)]
    e2f1_exclusive_idx += found.index.values.tolist()

for s in fb_only:
    found = orig[~pd.isnull(orig["NAME"])]
    found = found[found["NAME"].str.contains(s)]
    fb_exclusive_idx += found.index.values.tolist()

print "e2f1 exclusive: ", len(e2f1_exclusive_idx)
print "fb_exclusive: ", len(fb_exclusive_idx)

e2f1_data = orig.loc[e2f1_exclusive_idx]

fb_data = orig.loc[fb_exclusive_idx]

e2f1_syms = e2f1_data["NAME"].values
fb_syms = fb_data["NAME"].values[[0,1,4,5,6,7]]


e2f1_data = e2f1_data[good_cols]
e2f1_data.interpolate(axis=0, method="cubic", inplace=True)
e2f1_data.dropna(inplace=True, axis=0)
fb_data = fb_data[good_cols]
fb_data = fb_data.iloc[[0,1,4,5,6,7], :]
fb_data.interpolate(axis=0, method="cubic", inplace=True)
fb_data.dropna(inplace=True, axis=0)


e2f1_data_norm = cu.normalize_matrix_along_axis(e2f1_data.astype(float).values)
fb_data_norm = cu.normalize_matrix_along_axis(fb_data.astype(float).values)

e2f1link, e2f1reindex = cu.UPGMACluster(cu.calculateDistanceMatrix(e2f1_data_norm, axis=0), get_original_indexing=True)
fblink, fbreindex = cu.UPGMACluster(cu.calculateDistanceMatrix(fb_data_norm, axis=0), get_original_indexing=True)

if not os.path.isdir(os.path.join(data_dir, "cctf_hsf1_expression_figures")):
    os.mkdir(os.path.join(data_dir, "cctf_hsf1_expression_figures"))
else:
    shutil.rmtree(os.path.join(data_dir, "cctf_hsf1_expression_figures"))
    os.mkdir(os.path.join(data_dir, "cctf_hsf1_expression_figures"))


hm = pdh.DendroHeatMap(heat_map_data=e2f1_data_norm[e2f1reindex, :])
hm.left_dendrogram = e2f1link
hm.title = "E2F1 + HSF1"
hm.color_legend_title = "Microarray Expression Normalized Over Time"
hm.row_labels = [str(r.split("^")[1]) for r in e2f1_syms[e2f1reindex].tolist()]
hm.col_labels = time_labels
hm.export(filename=os.path.join(data_dir, "cctf_hsf1_expression_figures", "hsf1_e2f1_microarray.pdf"))
hm.show()

plt.plot(e2f1_data.mean())
plt.title("E2F1 HSF1 Mean Expression")
plt.xlabel("Time (h)")
plt.ylabel("Relative Expression")
plt.xticks([i*2 for i in range(len(e2f1_data.mean()))], time_labels)
plt.show()
plt.savefig(os.path.join(data_dir, "cctf_hsf1_expression_figures", "hsf1_e2f1_microarray_avg.pdf"))


hm.heat_map_data = fb_data_norm[fbreindex, :]
hm.left_dendrogram = fblink
hm.title = "FOXM1/BMYB + HSF1"
hm.row_labels = [str(r.split("^")[1]) for r in fb_syms[fbreindex].tolist()]
hm.col_labels = [str(c) for c in good_cols]
hm.color_legend_title = "Microarray Expression Normalized Over Time"
hm.export(filename=os.path.join(data_dir, "cctf_hsf1_expression_figures", "hsf1_foxm1bmyb_microarray.pdf"))
hm.show()

plt.plot(fb_data.mean())
plt.title("FOXM1/BMYB HSF1 Mean Expression")
plt.xlabel("Time (h)")
plt.ylabel("Relative Expression")
plt.show()
plt.savefig(os.path.join(data_dir, "cctf_hsf1_expression_figures", "hsf1_foxm1bmyb_microarray_avg.pdf"))


e2f1_data_argmax = e2f1_data.idxmax(axis=1)
print e2f1_data_argmax

e2f1_phases = [time_labels[x+1] for x in e2f1_data_argmax.values]
print e2f1_phases