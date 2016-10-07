import pandas as pd
import os
import chrom_utils
import pybedtools as pbt


data_dir = os.path.join("..", "data_files")
mendillo_merged_dir = os.path.join(data_dir, "mendillo_merged")
tss_bed_file = os.path.join(data_dir, "hg19_TSS.ranged.2500.bed")

tss_bed = pbt.BedTool(tss_bed_file)


mendillo_files = [os.path.join(mendillo_merged_dir,f) for f in os.listdir(mendillo_merged_dir) if "bed" in f]

# print "\n".join(mendillo_files)

output = []

for i, f1 in enumerate(mendillo_files):
    f1_bed = pbt.BedTool(f1)
    f1_cell = chrom_utils.get_cell_type_from_mendillo_filename(os.path.basename(f1))
    f1_tf = chrom_utils.get_tf_type_from_mendillo_filename(os.path.basename(f1))
    for j, f2 in enumerate(mendillo_files):
        if j<= i: continue
        else:
            f2_bed = pbt.BedTool(f2)
            f2_cell = chrom_utils.get_cell_type_from_mendillo_filename(os.path.basename(f2))
            f2_tf = chrom_utils.get_tf_type_from_mendillo_filename(os.path.basename(f2))

            if f2_tf != f1_tf:
                raise ValueError("TF types do not match: ", f1_tf, f2_tf)
            output.append((f1_cell + "-" + f1_tf,
                           f2_cell + "-" + f2_tf,
                           chrom_utils.calculate_TSS_co_binding_minimum_overlap_ratio(tss_bed=tss_bed,a_bed=f1_bed,b_bed=f2_bed)))


rows = []
cols = []

for x in output:
    if x[0] not in rows:
        rows.append(x[0])
    if x[1] not in cols:
        cols.append(x[1])

df = pd.DataFrame(index=rows, columns=cols)

print df
print rows
print cols

for x in output:
    df.set_value(x[0], x[1], x[2])

writer = pd.ExcelWriter(os.path.join(data_dir, "mendillo_cobinding_merged.xlsx"))
df.to_excel(writer)
writer.save()
