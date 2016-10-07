import os
import pybedtools as pbt

data_dir = os.path.join("..", "data_files")
hg19_tss_file = os.path.join(data_dir, "hg19_TSS.ranged.0.bed")

hg19_tss = pbt.BedTool(hg19_tss_file)

seen_features = []
uniques = []
for feature in hg19_tss:
    if "_" in feature.chrom: continue
    loc = ":".join(feature.fields[0:3])
    if loc not in seen_features:
        seen_features.append(loc)
        uniques.append(feature.fields)

print "number of unique features: ", len(uniques)

outfile = os.path.join(data_dir, "hg19_TSS.ranged.0.unique.bed")
unique_bed = pbt.BedTool(uniques)
unique_bed.sort().saveas(outfile)