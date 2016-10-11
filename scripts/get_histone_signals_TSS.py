import os
import cPickle as pickle
import pybedtools as pbt
import pyBigWig as bw
import numpy as np
from tqdm import tqdm
import shutil
import pandas as pd

data_dir = os.path.join("..", "data_files")
hg19_tss_file = os.path.join(data_dir, "hg19_TSS.ranged.1.bed")
histone_marks_dir = os.path.join(data_dir, "histone_marks")
mark_dirs = [os.path.join(histone_marks_dir, d) for d in os.listdir(histone_marks_dir) if "_" not in d]

hg19_tss = pbt.BedTool(hg19_tss_file).sort()

tss_output = {}

for d in mark_dirs:
    name = os.path.basename(d)
    mark_file = [os.path.join(d,f) for f in os.listdir(d) if ".bw" in f]
    assert len(mark_file) ==1, "didn't find one file: {d}, {l}".format(d=d, l=mark_file)
    mark = bw.open(mark_file[0])
    total_signal = mark.header()["sumData"]
    tss_output[name] = []
    # curdf = pd.DataFrame(columns=[int(i-1000) for i in range(2000)])
    print "Getting data for: ", name
    for feature in tqdm(hg19_tss):
        fields = feature.fields
        tss = int(fields[1])
        fieldname = fields[3]
        chr = fields[0]
        start = tss - 2000
        end = tss + 2000
        try:
            vals = np.array(mark.values(chr, start, end))
            vals *= 1e9
            vals /= total_signal
            tss_output[name].append((fieldname, vals))
        except RuntimeError:
            continue

    data = np.array([x[1] for x in tss_output[name]])
    curdf = pd.DataFrame(data=data, index=[x[0] for x in tss_output[name]], columns=[i-2000 for i in range(4000)])
    tss_output[name] = curdf

if not os.path.isdir(os.path.join(histone_marks_dir, "regions_output_files", "hg19_tss_marks")):
    os.mkdir(os.path.join(histone_marks_dir, "regions_output_files", "hg19_tss_marks"))
else:
    shutil.rmtree(os.path.join(histone_marks_dir, "regions_output_files", "hg19_tss_marks"))
    os.mkdir(os.path.join(histone_marks_dir, "regions_output_files", "hg19_tss_marks"))


store = pd.HDFStore(os.path.join(histone_marks_dir, "regions_output_files", "hg19_tss_marks.hdf5"))

for k in tss_output:
    curdf = tss_output[k]
    print "all naners?: ", curdf.isnull().values.all()
    store[k] = curdf
    # pickle.dump(tss_output[k], open( os.path.join(histone_marks_dir, "regions_output_files", "hg19_tss_marks", "hg19_tss_{m}_marks.pickle".format(m=k)), "w"))

store.close()