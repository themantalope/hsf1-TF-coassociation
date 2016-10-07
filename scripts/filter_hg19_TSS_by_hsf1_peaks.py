import pybedtools as pbt
import os

data_dir = os.path.join("..", "data_files")
hg19_tss_file = os.path.join(data_dir, "hg19_TSS.txt")
hsf1_bed_file = os.path.join(data_dir, "hsf1_0.17_merged_regions_hg19.bed")


#first need to make a hg19 tss bed file
hg19_tss_lines = open(hg19_tss_file, "r").readlines()

hg19_fields = []

for line in hg19_tss_lines:
    if "#" in line: continue
    else:
        parts = line.split("\t")
        chr = parts[1]
        start = stop = parts[3]
        strand = parts[2]
        name = parts[0]
        hg19_fields.append([chr, start, stop, name, strand, "0"])


hsf1_bed = pbt.BedTool(hsf1_bed_file)
hg19_tss = pbt.BedTool(hg19_fields).sort()
hg19_tss.saveas(os.path.join(data_dir, "hg19_TSS.ranged.0.bed"))

close = hg19_tss.sort().closest(hsf1_bed.sort(), D="a")
close_fields = []

for feature in close:
    if abs(int(feature[-1])) <= 1500 and int(feature[-1]) != -1:
        close_fields.append(feature.fields[0:5])

hg19_close = pbt.BedTool(close_fields).sort().saveas(os.path.join(data_dir, "hg19_TSS_filtered_hsf1_0.17_lifted.bed"))
