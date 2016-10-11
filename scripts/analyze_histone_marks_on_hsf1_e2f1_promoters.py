import os
import pybedtools as pbt
try: import cPickle as pickle
except: import pickle
import math
import matplotlib.pyplot as plt
import shutil
import pandas as pd
import numpy as np

data_dir = os.path.join("..", "data_files")
hg19_tss_file = os.path.join(data_dir, "hg19_TSS.ranged.1.bed")
hsf1_HS_regions_file = os.path.join(data_dir, "hsf1_HS_exclusive.bed")
hsf1_nonHS_regions_file = os.path.join(data_dir, "hsf1_nonHS_exclusive.bed")
e2f1_regions_file = os.path.join(data_dir, "cc_tf_bedfiles", "E2F1", "E2F1.bed")
tp53_regions_file = os.path.join(data_dir, "cc_tf_bedfiles", "TP53", "TP53.bed")
histone_marks_regions_dir = os.path.join(data_dir, "histone_marks", "regions_output_files")
tss_histone_file = os.path.join(histone_marks_regions_dir, "hg19_tss_marks.pickle")

hg19_tss = pbt.BedTool(hg19_tss_file).sort()
hsf1_nonHS_regions = pbt.BedTool(hsf1_nonHS_regions_file).sort()
hsf1_HS_regions = pbt.BedTool(hsf1_HS_regions_file).sort()
e2f1_regions = pbt.BedTool(e2f1_regions_file).sort()
tp53_regions = pbt.BedTool(tp53_regions_file).sort()

outf = os.path.join(histone_marks_regions_dir, "tfs_TSS_names.pickle")
if not os.path.isfile(outf):
    print "creating: ", outf
    #get the close TSS in each regions file
    hsf1_nonHS_close = hsf1_nonHS_regions.closest(hg19_tss, D="a", k=2)
    hsf1_nonHS_close_names = []
    chridxs = [i for i in range(len(hsf1_nonHS_close[0].fields)) if "chr" in hsf1_nonHS_close[0].fields[i] and "-" not in hsf1_nonHS_close[0].fields[i]]
    tss_chr_idx = chridxs[1]
    tss_name_idx = tss_chr_idx+3
    tss_strand_idx = -2
    for feature in hsf1_nonHS_close:
        fields = feature.fields
        if int(fields[-1]) == -1: continue
        hsf1_mid = int( math.floor( (int(fields[1]) + int(fields[2]))/2.0 ) )
        tss = int(fields[tss_chr_idx+1])
        dist = abs(hsf1_mid - tss)
        if dist <= 2000:
            hsf1_nonHS_close_names.append((fields[tss_name_idx], fields[tss_strand_idx]))



    #get the close names for hsf1 HS
    hsf1_HS_close = hsf1_HS_regions.closest(hg19_tss, D="a", k=2)
    hsf1_HS_close_names = []
    chridxs = [i for i in range(len(hsf1_HS_close[0].fields)) if "chr" in hsf1_HS_close[0].fields[i] and "-" not in hsf1_HS_close[0].fields[i]]
    tss_chr_idx = chridxs[1]
    tss_name_idx = tss_chr_idx+3
    for feature in hsf1_HS_close:
        fields = feature.fields
        if int(fields[-1]) == -1: continue
        hsf1_mid = int(math.floor( (int(fields[1]) + int(fields[2]))/2.0))
        tss = int(fields[tss_chr_idx+1])
        dist = abs(hsf1_mid - tss)

        if dist <= 2000:
            hsf1_HS_close_names.append((fields[tss_name_idx], fields[tss_strand_idx]))


    #get close names for e2f1
    e2f1_close = e2f1_regions.closest(hg19_tss, D="a", k=2)
    e2f1_close_names = []
    chridxs = [i for i in range(len(e2f1_close[0].fields)) if "chr" in e2f1_close[0].fields[i] and "-" not in e2f1_close[0].fields[i]]
    tss_chr_idx = chridxs[1]
    tss_name_idx = tss_chr_idx+3
    for feature in e2f1_close:
        fields = feature.fields
        if int(fields[-1]) == -1: continue
        e2f1_peak = int(fields[-1]) + int(fields[1])
        tss = int(fields[tss_chr_idx+1])
        dist = abs(e2f1_peak - tss)

        if dist <= 2000:
            e2f1_close_names.append((feature[tss_name_idx], fields[tss_strand_idx]))


    #get names for tp53

    tp53_close = tp53_regions.closest(hg19_tss, D="a", k=2)
    tp53_close_names = []
    chridxs = [i for i in range(len(tp53_close[0].fields)) if "chr" in tp53_close[0].fields[i] and "-" not in tp53_close[0].fields[i]]
    tss_chr_idx = chridxs[1]
    tss_name_idx = tss_chr_idx+3
    print "TP53"
    print tp53_close[0].fields
    print tss_chr_idx
    print tss_name_idx
    raw_input()
    for feature in tp53_close:
        fields = feature.fields
        if int(fields[-1]) == -1:continue
        tp53_peak = int(fields[-1]) + int(fields[1])
        tss = int(fields[tss_chr_idx+1])

        dist = abs(tp53_peak - tss)
        if dist <= 2000:
            tp53_close_names.append((fields[tss_name_idx], fields[tss_strand_idx]))

    names = {"HSF1_nonHS":hsf1_nonHS_close_names, "HSF1_HS":hsf1_HS_close_names, "E2F1":e2f1_close_names, "TP53":tp53_close_names}
    pickle.dump(names, open(outf, "wb"))
else:
    print "loading: ", outf
    names = pickle.load(open(outf))

print "Loading tss histone data..."

# hg19_tss_marks_dir = os.path.join(histone_marks_regions_dir, "hg19_tss_marks")
# tssfs = [os.path.join(hg19_tss_marks_dir,f) for f in os.listdir(hg19_tss_marks_dir) if "pickle" in f]
# tss_histone_data = {}
# for f in tssfs:
#     mark = os.path.basename(f).replace("_marks.pickle", "").replace("hg19_tss_", "")
#     tss_histone_data[mark] = pickle.load(open(f))


store = pd.HDFStore(os.path.join(histone_marks_regions_dir, "hg19_tss_marks.hdf5"))

tss_histone_data = {}
for k in store.keys():
    mykey = k.replace("/", "")
    tss_histone_data[mykey] = store[k]

# tss_histone_data = pickle.load(open(tss_histone_file))
print "Loaded tss histone data."



outfigures_dir = os.path.join(histone_marks_regions_dir, "composite_plots_tss")
if not os.path.isdir(outfigures_dir):
    os.mkdir(outfigures_dir)
else:
    shutil.rmtree(outfigures_dir)
    os.mkdir(outfigures_dir)

# for k in names:
#     print k
#     # if k == "TP53":continue
#     names_list = names[k]
#     avg_sigs = []
#     for hk in tss_histone_data:
#         names_list_copy = list(names_list)
#         tss_histone_list = tss_histone_data[hk]
#         signals = []
#         while len(names_list_copy) > 0:
#             curname, direction = names_list_copy.pop()
#             if curname not in tss_histone_data[hk].index: continue
#             histone_data_row = tss_histone_data[hk].loc[curname, :]
#             if len(histone_data_row.shape) > 1:
#                 histone_data_row = histone_data_row.mean(axis=0)
#             vals = histone_data_row.values
#             if direction == "-": vals = vals[::-1]
#             signals.append(vals)
#
#         signals = np.array(signals)
#         nsigs = signals.shape[0]
#         avg_sig = np.nansum(signals, axis=0)
#         avg_sig /= nsigs
#         avg_sigs.append((k,avg_sig))
#
#     fig, ax = plt.subplots(nrows=1, ncols=1)
#     for asig in avg_sigs:
#         ax.plot([i-asig[1].size/2.0 for i in range(asig[1].size)], asig[1], label=asig[0])
#
#     ax.set_xlabel("Distance from TSS")
#     ax.set_ylabel("Reads per billion reads")
#     ax.set_title("{m} Average Signal at TSS".format(m=hk))
#     plt.show()
#     fig.savefig(os.path.join(outfigures_dir, "{m}_tss_composite_plot.pdf".format(m=hk)))
#     plt.show()


for hk in tss_histone_data:
    avg_sigs = []

    for k in names:
        signals = []
        names_list_copy = list(names[k])
        while len(names_list_copy) > 0 :
            curname, direction = names_list_copy.pop()
            if curname not in tss_histone_data[hk].index: continue
            histone_data_row = tss_histone_data[hk].loc[curname, :]
            if len(histone_data_row.shape) > 1:
                histone_data_row = histone_data_row.mean(axis=0)
            vals = histone_data_row.values
            if direction == "-": vals = vals[::-1]
            signals.append(vals)

        signals = np.array(signals)
        nsigs = signals.shape[0]
        avg_sig = np.nansum(signals, axis=0)
        avg_sig /= nsigs
        avg_sigs.append((k, avg_sig))

    fig, ax = plt.subplots(nrows=1, ncols=1)
    for asig in avg_sigs:
        ax.plot([i-asig[1].size/2.0 for i in range(asig[1].size)], asig[1], label=asig[0])

    ax.set_xlabel("Distance from TSS")
    ax.set_ylabel("Reads per billion reads")
    ax.set_title("{m} Average Signal at TSS".format(m=hk))
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    fig.savefig(os.path.join(outfigures_dir, "{m}_tss_composite_plot.pdf".format(m=hk)))
    # plt.show()


