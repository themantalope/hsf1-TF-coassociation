import os
import chrom_utils
import pybedtools as pbt
import re
import shutil

data_dir = os.path.join("..", "data_files")
merged_dir = os.path.join(data_dir, "mendillo_merged")
tss_file = os.path.join(data_dir, "hg19_TSS.ranged.2500.bed")
hs_dir = os.path.join(data_dir, "mendillo_HS")
hs_files = [os.path.join(hs_dir,f) for f in os.listdir(hs_dir) if os.path.isfile(os.path.join(hs_dir, f)) and "+" in f]


old_fs = [os.path.join(hs_dir, f) for f in os.listdir(hs_dir) if "peak." in f]
for f in old_fs:
    os.remove(f)


tss_bed = pbt.BedTool(tss_file)

hs_cell_type_patt = re.compile("(?<=^)(.*?)(?=\+)")

#first compare the overall association between the two files

cell_types = set()
for f in hs_files:
    cell_types.add(re.findall(hs_cell_type_patt, os.path.basename(f))[0])

print "Cell types: ", list(cell_types)
print "\nComparing HS Replicates"
for t in cell_types:
    tfiles = set()
    for f in hs_files:
        if t in f: tfiles.add(f)

    tfiles = list(tfiles)
    if len(tfiles) != 2: raise ValueError("Didn't get the right number of files: {s}".format(s=tfiles))
    print "TSS overlap for {c}: {v}".format(c=t, v=chrom_utils.calculate_TSS_co_binding_minimum_overlap_ratio(tss_bed, pbt.BedTool(tfiles[0]), pbt.BedTool(tfiles[1])))
    print "Overall overlap for {c}: {v}".format(c=t, v=float(len(pbt.BedTool(tfiles[0]).intersect(pbt.BedTool(tfiles[1])))) / float(min(len(pbt.BedTool(tfiles[0])), len(pbt.BedTool(tfiles[1])))) )

    #next merge the two files together
    merged_name = os.path.join(hs_dir,"{t}_HEATSHOCK_Regions_HS_HSF1.bed".format(t=t))
    merged_bed = pbt.BedTool(tfiles[0]).cat(pbt.BedTool(tfiles[1]),postmerge=False, force_truncate=False)
    merged_bed.sort().saveas(merged_name)


#now we need to reset the peaks for each file

new_hs_files = [os.path.join(hs_dir,f) for f in os.listdir(hs_dir) if "_HEATSHOCK" in f]

for f in new_hs_files:
    chrom_utils.set_region_start_stop_based_on_peak(f)


#compare the amount of overlap to the non heatshock files

non_hs_files = [os.path.join(merged_dir,f) for f in os.listdir(merged_dir) if ".peak.bed" in f and chrom_utils.get_cell_type_from_mendillo_filename(f) in cell_types]

# print non_hs_files
print "\n"
print "Comparing HS to non-HS"
for f in non_hs_files:
    non_hs_bed = pbt.BedTool(f)
    cell = chrom_utils.get_cell_type_from_mendillo_filename(os.path.basename(f))
    hs_bed_file = [h for h in new_hs_files if cell in h][0]
    hs_bed = pbt.BedTool(hs_bed_file)

    print "Overall overlap for {c}: {v}".format(c=cell, v= float(len(hs_bed.intersect(non_hs_bed).merge())) / float(min(len(hs_bed), len(non_hs_bed))))
    print "TSS overlap for {c}: {v}".format(c=cell, v=chrom_utils.calculate_TSS_co_binding_minimum_overlap_ratio(tss_bed, hs_bed, non_hs_bed))

