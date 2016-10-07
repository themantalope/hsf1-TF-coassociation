import os
import chrom_utils
import pybedtools as pbt
import re
import pandas as pd

ext_patt = re.compile("(.bed)(?=$)")

data_dir = os.path.join("..", "data_files")
cc_tf_dir = os.path.join(data_dir, "cc_tf_bedfiles")
tss_bed_file = os.path.join(data_dir, "hg19_TSS.ranged.bed")

mendillo_HSF1_file = os.path.join(data_dir, "mendillo_merged", "MCF10_Regions_HSF1.peak.bed")

cc_tf_files = []
cc_tf_dirs = [os.path.join(cc_tf_dir,d) for d in os.listdir(cc_tf_dir) if os.path.isdir(os.path.join(cc_tf_dir,d))]

for d in cc_tf_dirs:
    ccfs = [os.path.join(d,f) for f in os.listdir(d) if os.path.isfile(os.path.join(d,f)) and len(re.findall(ext_patt, f)) != 0 and "peak" not in f]
    cc_tf_files += ccfs

all_files = cc_tf_files

all_file_tuple = []
for i, f1 in enumerate(all_files):
    if "RB" in f1 or "MYC" in f1: continue
    all_file_tuple.append((mendillo_HSF1_file, f1))



tss_bed = pbt.BedTool(tss_bed_file)

min_ovr_ratio = chrom_utils.calculate_TSS_co_binding_minimum_overlap_ratio_multiprocessing(tss_bed, all_file_tuple, useflank=None)
# min_ovr_ratio = []

# for t in all_file_tuple:
#     min_ovr_ratio.append(chrom_utils.calculate_TSS_co_binding_minimum_overlap_ratio(tss_bed, pbt.BedTool(t[0]), pbt.BedTool(t[1]), useflank=None))


final_out = []
for i in range(0, len(all_file_tuple)):
    t = all_file_tuple[i]
    ovr = min_ovr_ratio[i]
    final_out.append((os.path.basename(t[0]), os.path.basename(t[1]), ovr))

outdf = chrom_utils.convert_tuple_list_to_dataframe(final_out)
writer = pd.ExcelWriter(os.path.join(data_dir, "hsf1_vs_cc_tf.xlsx"))
outdf.to_excel(writer)
writer.save()
print "Done"