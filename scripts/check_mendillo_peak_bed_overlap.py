import pandas as pd
import os
import pickle
import pybedtools as pbt
import chrom_utils


data_dir = os.path.join("..", "data_files")
tss_bed_file = os.path.join(data_dir, "hg19_TSS.ranged.2500.bed")
mendillo_list_file = os.path.join(data_dir, "mendillo_peak_bed_files.pickle")

mendillo_list = pickle.load(open(mendillo_list_file, "r"))

tss_bed = pbt.BedTool(tss_bed_file)

output = []

for i, f1 in enumerate(mendillo_list):
    f1_bt = pbt.BedTool(f1)
    f1_cell = chrom_utils.get_cell_type_from_mendillo_filename(os.path.basename(f1))
    f1_tf = chrom_utils.get_tf_type_from_mendillo_filename(os.path.basename(f1))
    for j, f2 in enumerate(mendillo_list):
        if i <= j: continue
        else:
            f2_bt = pbt.BedTool(f2)
            f2_cell = chrom_utils.get_cell_type_from_mendillo_filename(os.path.basename(f2))
            f2_tf = chrom_utils.get_tf_type_from_mendillo_filename(os.path.basename(f2))

            output.append((f1_cell+"-"+f1_tf,
                           f2_cell+"-"+f2_tf,
                           chrom_utils.calculate_TSS_co_binding_minimum_overlap_ratio(tss_bed, f1_bt, f2_bt)))


rows = []
cols = []

for x in output:
    if x[0] not in rows:
        rows.append(x[0])
    if x[1] not in cols:
        cols.append(x[1])


df = pd.DataFrame(index=rows, columns=cols)

writer = pd.ExcelWriter(os.path.join(data_dir, "mendillo_peak_bed_method_check.xlsx"))

for x in output:
    df.set_value(x[0], x[1], x[2])

df.to_excel(writer)
writer.save()
