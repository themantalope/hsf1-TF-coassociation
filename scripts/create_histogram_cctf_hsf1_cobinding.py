import os
import matplotlib.pyplot as plt
import pandas as pd
import shutil

data_dir = os.path.join("..", "data_files")
cctf_file = os.path.join(data_dir, "hsf1_cctf_position_analysis.xlsx")

cctf = pd.read_excel(open(cctf_file), sheetname=None)
histo_figures = os.path.join(data_dir, "cctf_hsf1_histograms")

if not os.path.isdir(histo_figures):
    os.mkdir(histo_figures)
else:
    shutil.rmtree(histo_figures)
    os.mkdir(histo_figures)

for c in cctf.keys():
    curdf = cctf[c]
    dist_col = "HSF1 - {tf}".format(tf=c)
    tflocdata = cctf[c]
    dists = tflocdata[dist_col]
    dists = dists[~pd.isnull(dists)].values
    if len(dists) < 50: continue
    plt.hist(dists, bins=50)
    plt.title("HSF1 - {tf} Binding distances".format(tf=c))
    # plt.show()
    plt.savefig(os.path.join(histo_figures, "hsf1_{tf}_association_dist.pdf".format(tf=c)))
