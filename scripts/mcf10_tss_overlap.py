import pybedtools as pbt
import os
import pickle
import chrom_utils
import pandas as pd

data_dir = os.path.join("..", "data_files")
tss_file = os.path.join(data_dir, "hg19_TSS.ranged.2500.bed")

mendillo_list_file = os.path.join(data_dir, "mendillo_peak_bed_files.pickle")
encode_list_file = os.path.join(data_dir, "encode_peak_bed_files.pickle")

mendillo_files = [f for f in pickle.load(open(mendillo_list_file)) if "MCF10" in f or "Mcf10" in f]
encode_files = [f for f in pickle.load(open(encode_list_file)) if "MCF10" in f or "Mcf10" in f]

all_files = mendillo_files + encode_files

tss_bed = pbt.BedTool(tss_file)

output = []
for f in all_files:
    fbt = pbt.BedTool(f)
    ratio = float(len(tss_bed.intersect(fbt).sort().merge())) / float(len(fbt))
    output.append(ratio)

final_out = []
for i, f in enumerate(all_files):
    final_out.append((chrom_utils.get_cell_type_from_filename(os.path.basename(f)) + "-" + chrom_utils.get_tf_type_from_filename(os.path.basename(f)) + "-" + chrom_utils.get_exp_type_from_filename(os.path.basename(f)), output[i]))

rows = [x[0] for x in final_out]
ratios = [x[1] for x in final_out]

# print final_out
# print rows
# print ratios

df = pd.DataFrame(data=ratios, index=rows, columns=["TSS Overlap Fraction"])

writer = pd.ExcelWriter(os.path.join(data_dir, "mcf10_tss_overlap.xlsx"))

df.to_excel(writer)
writer.save()