import os
import pybedtools as pbt
import math

data_dir = os.path.join("..", "data_files")
hg19_tss_file = os.path.join(data_dir, "hg19_TSS.ranged.0.unique.bed")

hg19_tss = pbt.BedTool(hg19_tss_file)

newregions = []
for feature in hg19_tss:
    fields = feature.fields
    newstart = int(fields[1]) - 2000
    newstop = int(fields[2]) + 2000
    if newstart < 0: newstart = 1

    fields[1] = str(newstart)
    fields[2] = str(newstop)
    newregions.append(fields)

hg19_tss_4000 = pbt.BedTool(newregions).sort()

hg19_tss_merged = hg19_tss_4000.merge()

print len(hg19_tss)
print len(hg19_tss_4000)
print len(hg19_tss_merged)


bp_500_counts = 0
for feature in hg19_tss_merged:
    reg_length = abs(int(feature.fields[1]) - int(feature.fields[2]))
    bp_500_counts += int( math.floor( float(reg_length) / float(500) ) )

print "number of 500 bp regions: ", bp_500_counts