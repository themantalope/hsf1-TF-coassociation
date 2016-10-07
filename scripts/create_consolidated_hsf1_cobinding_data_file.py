import os
import pandas as pd
import datetime


data_dir = os.path.join("..", "data_files")
cctf_file = os.path.join(data_dir, "hsf1_HS_nonHS_cctf_feature_comparisons.xlsx")
encode_file = os.path.join(data_dir, "hsf1_HS_nonHS_encode_feature_comparisons.xlsx")

cctf = pd.read_excel(open(cctf_file))
encode = pd.read_excel(open(encode_file))

encode_groups = encode.groupby("TF").groups

outidxs = []

for g in encode_groups:
    gidxs = encode_groups[g]
    encode_g_df = encode.loc[gidxs, :]
    encode_g_means = encode_g_df.mean()
    cols = encode_g_means.index.values.tolist()
    print encode_g_means
    print cols

    g0idx = gidxs[0]
    encode.loc[g0idx, cols] = encode_g_means
    outidxs.append(g0idx)

outencode = encode.loc[outidxs, :]

outencode.set_index("TF", inplace=True)

consolidated = pd.concat([outencode, cctf])

print consolidated.columns

# print consolidated
outwriter = pd.ExcelWriter(os.path.join(data_dir, "hsf1_HS_nonHS_feature_comparisons_combined.v2.xlsx"))
consolidated.to_excel(outwriter)
outwriter.save()