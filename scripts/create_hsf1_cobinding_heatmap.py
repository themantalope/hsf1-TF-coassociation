import pandas as pd
import os
import pydendroheatmap as pdh
import maakclusterutils as cu
import numpy as np

data_dir = os.path.join("..", "data_files")
data_file = os.path.join(data_dir, "hsf1_HS_nonHS_feature_comparisons_combined.v2.xlsx")

data = pd.read_excel(open(data_file))

data.set_index("TF", inplace=True)

row_labels = data.index.values.tolist()
col_labels = [c for c in data.columns if "OR" in c]

# print row_labels
# print col_labels

pv_data = data[[c for c in data.columns if "OR" in c]]

log_pv_data = pv_data.values

# print log_pv_data

rowlink, row_reindex = cu.UPGMACluster(cu.calculateDistanceMatrix(log_pv_data,axis=0), get_original_indexing=True)
collink, col_reindex = cu.UPGMACluster(cu.calculateDistanceMatrix(log_pv_data,axis=1), get_original_indexing=True)


log_pv_data = log_pv_data[:, col_reindex][row_reindex, :]


row_labels = np.array(row_labels)[row_reindex].tolist()
col_labels = np.array(col_labels)[col_reindex].tolist()



hm = pdh.DendroHeatMap(log_pv_data,
                       left_dendrogram=rowlink,
                       top_dendrogram=collink,
                       row_labels=row_labels,
                       col_labels=col_labels,
                       window_height=8)

hm.show()
hm.export(os.path.join(data_dir, "hsf1_HS_nHS_tf_cobinding_heatmap.pdf"))