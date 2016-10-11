import os
import pyBigWig as bw
import pandas as pd
import numpy as np
import pybedtools as pbt
import math
import shutil
from tqdm import tqdm
import cPickle as pickle

data_dir = os.path.join("..", "data_files")
hsf1_nonHS_regions_file = os.path.join(data_dir, "hsf1_nonHS_exclusive.bed")
hsf1_HS_regions_file = os.path.join(data_dir, "hsf1_HS_exclusive.bed")
hg19_tss_file = os.path.join(data_dir, "hg19_TSS.ranged.1.bed")
e2f1_regions_file = os.path.join(data_dir, "cc_tf_bedfiles", "E2F1", "E2F1.bed")
tp53_regions_file = os.path.join(data_dir, "cc_tf_bedfiles", "TP53", "TP53.bed")

histone_dir = os.path.join(data_dir, "histone_marks")

histone_files = []
for d in os.listdir(histone_dir):
    if "output" in d: continue
    d = os.path.join(histone_dir, d)
    if not os.path.isdir(d): continue
    else:
        dfiles = [os.path.join(d,f) for f in os.listdir(d) if ".bw" in f]
        assert  len(dfiles) == 1, "didn't find one bw file for {d}: {l}".format(d=d, l=dfiles)
        histone_files.append(dfiles[0])


#for now, use the center of the HSF1 regions as the peak
hsf1_HS_regions = pbt.BedTool(hsf1_HS_regions_file)
hsf1_HS_regions_with_peak_fields = []
for feature in hsf1_HS_regions:
    fields = feature.fields
    start = int(fields[1])
    end = int(fields[2])
    mid = int(math.floor(  (start + end)/2.0  ) - start)
    fields.append(mid)
    hsf1_HS_regions_with_peak_fields.append(fields)

hsf1_HS_regions = pbt.BedTool(hsf1_HS_regions_with_peak_fields).sort()

#repeat for the nonHS HSF1 regions
hsf1_nonHS_regions = pbt.BedTool(hsf1_nonHS_regions_file)
hsf1_nonHS_regions_with_peak_fields = []
for feature in hsf1_nonHS_regions:
    fields = feature.fields
    start = int(fields[1])
    end = int(fields[2])
    mid = int(math.floor( (start + end)/2.0 ) - start)
    fields.append(mid)
    hsf1_nonHS_regions_with_peak_fields.append(fields)

hsf1_nonHS_regions = pbt.BedTool(hsf1_nonHS_regions_with_peak_fields).sort()


hg19_tss = pbt.BedTool(hg19_tss_file).sort()
#for the other TFs, get only the peaks near TSS
e2f1_regions = pbt.BedTool(e2f1_regions_file).sort()
e2f1_close = e2f1_regions.closest(hg19_tss, D="a", k=2)
e2f1_close_fields = []
chridxs = [i for i in range(len(e2f1_close[0].fields)) if "chr" in e2f1_close[0].fields[i]]
e2f1_peak_idx = chridxs[1]-1
tss_start_idx = chridxs[1]+1
uniques = []
for feature in e2f1_close:
    fields = feature.fields
    if int(feature.fields[-1]) == -1: continue
    e2f1_peak_pos = int(feature.fields[1]) + int(feature.fields[e2f1_peak_idx])
    tss = int(feature.fields[tss_start_idx])

    dist = abs(e2f1_peak_pos - tss)
    if dist <= 2000:
        if ":".join(fields[0:3]) not in uniques:
            e2f1_close_fields.append(fields[0:e2f1_peak_idx+1])
            uniques.append(":".join(fields[0:3]))


e2f1_regions = pbt.BedTool(e2f1_close_fields).sort()




#repeat for tp53
tp53_regions = pbt.BedTool(tp53_regions_file).sort()
tp53_close = tp53_regions.closest(hg19_tss, D="a", k=2)
tp53_close_fields = []
chridxs = [i for i in range(len(tp53_close[0].fields)) if "chr" in tp53_close[0].fields[i]]
tp53_peak_idx = chridxs[1]-1
tss_start_idx = chridxs[1]+1
uniques = []
for feature in tp53_close:
    fields = feature.fields
    if int(feature.fields[-1]) == -1: continue

    tp53_peak_pos = int(feature.fields[1]) + int(feature.fields[tp53_peak_idx])
    tss = int(feature.fields[tss_start_idx])

    dist = abs(tp53_peak_pos - tss)
    if dist <= 2000:
        if ":".join(fields[0:3]) not in uniques:
            tp53_close_fields.append(fields[0:tp53_peak_idx + 1])
            uniques.append(":".join(fields[0:3]))

tp53_regions = pbt.BedTool(tp53_close_fields).sort()





tf_regions_bts = [("HSF1 HS", hsf1_HS_regions),
                  ("HSF1 nonHS", hsf1_nonHS_regions),
                  ("E2F1", e2f1_regions),
                  ("TP53", tp53_regions)]

outdf = {}
tf_output_dir = os.path.join(data_dir, "histone_marks", "regions_output_files")
if not os.path.isdir(tf_output_dir):
    os.mkdir(tf_output_dir)
else:
    shutil.rmtree(tf_output_dir)
    os.mkdir(tf_output_dir)


for r in tf_regions_bts:
    tf = r[0]
    tf_bed = r[1].remove_invalid().sort()
    tf_output = {}
    for mark in histone_files:
        name = os.path.basename(mark)
        name = name[0:name.find(".")]
        tf_output[name] = []
        histone_bw = bw.open(mark)
        total_signal = histone_bw.header()["sumData"]
        for feature in tqdm(tf_bed):
            peak_pos = int(feature.fields[-1]) + int(feature.fields[1])
            start = peak_pos - 250
            end = peak_pos + 250
            chr = feature.fields[0]
            try:
                vals = np.array(histone_bw.values(chr, start, end))
                vals *= 1e9
                vals /= total_signal
                tf_output[name].append(("-".join([chr, str(start), str(end)]), vals))
            except RuntimeError:
                continue

    pick_file = os.path.join(tf_output_dir, "{r}_marks.pickle".format(r=tf.replace(" ", "_")))
    pickle.dump(tf_output, open(pick_file, "wb"))

    # for k in tf_output:
    #
    #     pickle.dump(tf_output)
    #     cols = [i-250 for i in range(len(tf_output[k][0][1]))]
    #     curdf = pd.DataFrame(index=[x[0] for x in tf_output[k]], columns=cols)
    #     for x in tf_output[k]:
    #         curdf.loc[x[0], cols] = x[1]
    #
    #     curdf.to_pickle(os.path.join(tf_output_dir, "{r}_{n}_marks.pickle".format(r=tf.replace(" ", "_"), n=k)))
    #
    #
