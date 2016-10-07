import os
import matplotlib.pyplot as plt
import pandas as pd

data_dir = os.path.join("..", "data_files")
cctf_file = os.path.join(data_dir, "hsf1_features_cctf.xlsx")

cctf = pd.read_excel(open(cctf_file))

for c in cctf.columns:
    tflocdata = cctf[c]
    dists = tflocdata[tflocdata.index.isin(["Co-Occurance Enrichment", "My Odds Ratio", "Fisher Exact Odds Ratio", "Fisher Exact P-value"])]
    print dists