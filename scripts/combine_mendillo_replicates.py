import shutil
import os
import pybedtools as pbt
import re
import chrom_utils

data_dir = os.path.join("..", "data_files")
mendillo_dir = os.path.join(data_dir, "mendillo")
mendillo_merged_dir = os.path.join(data_dir, "mendillo_merged")



if not os.path.isdir(mendillo_merged_dir):
    os.mkdir(mendillo_merged_dir)
else:
    shutil.rmtree(mendillo_merged_dir)
    os.mkdir(mendillo_merged_dir)

rep_pattern = re.compile("\_R\d")

mendillo_files = [os.path.join(mendillo_dir, f) for f in os.listdir(mendillo_dir) if "bed" in f and re.search(rep_pattern, f)]
# print mendillo_files
#now group files by cell type

other_files = [os.path.join(mendillo_dir, f) for f in os.listdir(mendillo_dir) if "bed" in f and not re.search(rep_pattern, f)]

for i,f1 in enumerate(mendillo_files):
    if "Colon" in f1 or "Breast" in f1: continue
    else:
        f1_cell = chrom_utils.get_cell_type_from_mendillo_filename(os.path.basename(f1))
        f1_cell = re.sub(rep_pattern, "", f1_cell)
        f1_cell = f1_cell.replace("_", "")
        f2 = [f for f in mendillo_files if f1_cell in os.path.basename(f) and f != f1]
        if len(f2) != 1:
            raise ValueError("Didn't get one other file!: ", f2, f1, f1_cell)
        f2 = f2[0]
        f1_bed = pbt.BedTool(f1)
        f2_bed = pbt.BedTool(f2)

        merge_file = os.path.join(mendillo_merged_dir, "{c}_Regions_HSF1.bed".format(c=f1_cell))
        x=pbt.BedTool()
        f1_bed.cat(f2_bed, postmerge=False, force_truncate=False).saveas(merge_file)

        merge_bed = pbt.BedTool(merge_file)
        merge_bed.sort().saveas(merge_file)


for f in other_files:
    shutil.copy(f, os.path.join(mendillo_merged_dir, os.path.basename(f)))

print "Done merging"