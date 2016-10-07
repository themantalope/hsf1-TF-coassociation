import os
import pybedtools as pbt
import chrom_utils
import pandas as pd

data_dir = os.path.join("..", "data_files")
mendillo_unmerged_dir = os.path.join(data_dir, "mendillo")
mendillo_merged_dir = os.path.join(data_dir, "mendillo_merged")


replicate_files = [os.path.join(mendillo_unmerged_dir,f) for f in os.listdir(mendillo_unmerged_dir) if ("R1" in f or "R2" in f) and ("Colon" not in f and "Breast" not in f)]

rep_pairs = [{"1":f, "2":""} for f in replicate_files if "R1" in f]
cell_types = set()
for rp in rep_pairs:
    r1 = rp["1"]
    r2 = r1.replace("R1", "R2")
    rp["2"] = r2
    cell_types.add(chrom_utils.get_cell_type_from_filename(os.path.basename(r1)))

# print rep_pairs

df = pd.DataFrame(index=cell_types, columns=["# Peaks R1", "# Peaks R2", "# Within 2500 bp", "# Merged", "# Intersect"])

for rp in rep_pairs:
    r1_bed = pbt.BedTool(rp["1"]).sort()
    r2_bed = pbt.BedTool(rp["2"]).sort()

    for ct in cell_types:
        if ct in rp["1"]: cur_cell_type = ct;break
    closest = r1_bed.closest(r2_bed, d=True, t="first")
    good_peaks = []
    for feature in closest:
        if int(feature[-1]) < 2500:
            good_peaks.append(feature)

    out_merged_name = rp["1"].replace(".bed", ".merged.bed")
    out_merged_name = os.path.join(mendillo_merged_dir, os.path.basename(out_merged_name.replace("_R1", "")))

    if len(r1_bed) > len(r2_bed):
        m = r2_bed.merge(i=r1_bed, d=200, o="distinct",c=4).remove_invalid().sort().saveas(out_merged_name)
    else:
        m = r1_bed.merge(i=r2_bed, d=200, o="distinct", c=4).remove_invalid().sort().saveas(out_merged_name)
    intersection = r1_bed.intersect(r2_bed)
    df.set_value(index=cur_cell_type, col="# Peaks R1", value=len(r1_bed))
    df.set_value(index=cur_cell_type, col="# Peaks R2", value=len(r2_bed))
    df.set_value(index=cur_cell_type, col="# Within 2500 bp", value=len(good_peaks))
    df.set_value(index=cur_cell_type, col="# Merged", value=len(m))
    df.set_value(index=cur_cell_type, col="# Intersect", value=len(intersection))


writer = pd.ExcelWriter(os.path.join(data_dir, "mendillo_replicate_analysis.xlsx"))
df.to_excel(writer)
writer.save()