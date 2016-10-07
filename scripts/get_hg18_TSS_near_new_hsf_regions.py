import os
import pybedtools as pbt

data_dir = os.path.join("..", "data_files")
new_hsf_bed_file = os.path.join(data_dir, "hsf1_0.15_merged_regions_hse.bed")
hg18_TSS_bed_file = os.path.join(data_dir, "hg18_TSS.ranged.0.bed")

hg18_TSS_bed = pbt.BedTool(hg18_TSS_bed_file)
new_hsf_bed = pbt.BedTool(new_hsf_bed_file)

print len(new_hsf_bed)

# for feature in new_hsf_bed:
#     print abs(int(feature[1]) - int(feature[2]))

closest = hg18_TSS_bed.sort().closest(new_hsf_bed.sort(),D="a")
close_genes = []

for feature in closest:
    if abs(int(feature[-1])) <= 1500 and int(feature[-1]) != -1:
        # print feature.fields
        close_genes.append(feature.fields)

print len(close_genes)
print len(closest)



close_genes = pbt.BedTool(close_genes)
gene_set = set()

for feature in close_genes:
    gene_set.add(feature[1])
    # print feature.fields


out_tss_list = os.path.join(data_dir, "hg18_TSS_filtered_by_hsf_0.15_hse.bed")

close_genes.saveas(out_tss_list)