import pybedtools as pbt
import os
import chrom_utils
import pandas as pd

data_dir = os.path.join("..", "data_files")
mendillo_HS_dir = os.path.join(data_dir, "mendillo_HS")
tss_file = os.path.join(data_dir, "hg19_TSS.ranged.2500.bed")

tss_bed = pbt.BedTool(tss_file)

hs_files = [os.path.join(mendillo_HS_dir, f) for f in os.listdir(mendillo_HS_dir) if "bed" in f]

output = {}

for f in hs_files:
    hs_bed = pbt.BedTool(f)
    cell = chrom_utils.get_cell_type_from_mendillo_filename(os.path.basename(f))
    tss_intersection = tss_bed.intersect(hs_bed)
    tss_names = [tss_intersection[n].name for n in range(len(tss_intersection))]
    total_intersection_len = float(len(tss_intersection))
    total_hs_binding_sites = float(len(hs_bed))
    promoter_binding_rate = total_intersection_len / total_hs_binding_sites
    output[cell] = {"genes":tss_names, "tss_intersection":total_intersection_len, "hs_binding_sites":total_hs_binding_sites, "promoter_binding_rate":promoter_binding_rate}

writer = pd.ExcelWriter(os.path.join(data_dir, "hsf1_HS_binding.xlsx"))
for c in output:
    df = pd.DataFrame(output[c])
    df.to_excel(writer,sheet_name=c)

writer.save()


