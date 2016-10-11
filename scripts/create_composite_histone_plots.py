import matplotlib.pyplot as plt
import os
import pickle
import numpy as np
import shutil

data_dir = os.path.join("..","data_files")
histone_data_dir = os.path.join(data_dir, "histone_marks", "regions_output_files")

dfs = [os.path.join(histone_data_dir,f) for f in os.listdir(histone_data_dir) if "pickle" in f]

output_dir = os.path.join(data_dir, "histone_marks", "composite_plots")
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)
else:
    shutil.rmtree(output_dir)
    os.mkdir(output_dir)

for f in dfs:
    name = os.path.basename(f).replace("_marks.pickle", "")
    data = pickle.load(open(f))
    for k in data:
        n_regions = float(len(data[k]))
        signals = [x[1] for x in data[k]]
        avg_sig = np.nansum(np.array(signals), axis=0)
        avg_sig /= n_regions

        plt.plot([i-250 for i in range(len(avg_sig))],avg_sig)
        plt.title("{n} and {m}".format(n=name, m=k))
        plt.xlabel("Distance from ChIP-seq Peak")
        plt.ylabel("Read per billion basepair")
        plt.savefig(os.path.join(output_dir, "{n}_{m}_composite_plot.pdf".format(n=name, m=k)))
        plt.show()

