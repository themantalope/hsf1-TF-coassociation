import pandas as pd
import os
import pybedtools as pbt



data_dir = os.path.join("..", "data_files")
hsf1_regions_file = os.path.join(data_dir, "hsf1_0.17_merged_regions_hg19.bed")
hsf1_nonHS_regions_file = os.path.join(data_dir, "hsf1_nonHS_exclusive.bed")
hsf1_HS_regions_file = os.path.join(data_dir, "hsf1_HS_exclusive.bed")
e2f1_regions_file = os.path.join(data_dir, "cc_tf_bedfiles", "E2F1", "E2F1.bed")
hg19_tss_file = os.path.join(data_dir, "hg19_TSS.ranged.0.unique.bed")

hsf1_regions = pbt.BedTool(hsf1_regions_file)
hsf1_nonHS_regions = pbt.BedTool(hsf1_nonHS_regions_file)
hsf1_HS_regions = pbt.BedTool(hsf1_HS_regions_file)
e2f1_regions = pbt.BedTool(e2f1_regions_file)
hg19_tss = pbt.BedTool(hg19_tss_file)

#get the number of e2f1 peaks near the TSS
e2f1_tss_close = hg19_tss.sort().closest(e2f1_regions.sort(), D="a")
uniques = []
count = 0
for feature in e2f1_tss_close:
    if int(feature.fields[-1]) == -1: continue

    # print feature.fields
    e2f1_start = int(feature.fields[-10])
    e2f1_peak = int(feature.fields[-2])

    e2f1_peak_pos = e2f1_start + e2f1_peak
    dist = min( abs(int(feature.fields[1]) - e2f1_peak_pos), abs(int(feature.fields[2]) - e2f1_peak_pos) )
    if dist <= 2000:
        if ":".join(feature.fields[0:3]) not in uniques:
            count += 1
            uniques.append(":".join(feature.fields[0:3]))

#compute the number of e2f1 peaks within 500 bp of HSF1 in the various contexts

e2f1_hsf1_regions_close = hsf1_regions.sort().closest(e2f1_regions.sort(), D="a")
uniques = []
hsf1_count = 0
for feature in e2f1_hsf1_regions_close:
    if int(feature.fields[-1]) == -1: continue
    e2f1_start = int(feature.fields[-10])
    e2f1_peak = int(feature.fields[-2])
    e2f1_peak_pos = e2f1_start + e2f1_peak

    hsf1_mid = (int(feature.fields[1]) + int(feature.fields[2])) / 2.0

    # print feature.fields
    # print e2f1_peak_pos
    # raw_input()
    dist = min(abs(hsf1_mid - e2f1_peak_pos), abs(hsf1_mid - e2f1_peak_pos))
    if dist <= 500:
        if ":".join(feature.fields[0:3]) not in uniques:
            hsf1_count += 1
            uniques.append(":".join(feature.fields[0:3]))


#repeat for nonHS
e2f1_hsf1_nonHS_regions_close = hsf1_nonHS_regions.sort().closest(e2f1_regions.sort(), D="a")
uniques = []
hsf1_nonHS_count = 0
for feature in e2f1_hsf1_nonHS_regions_close:
    if int(feature.fields[-1]) == -1: continue
    e2f1_start = int(feature.fields[-10])
    e2f1_peak = int(feature.fields[-2])
    e2f1_peak_pos = e2f1_start + e2f1_peak

    hsf1_mid = (int(feature.fields[1]) + int(feature.fields[2])) / 2.0

    # print feature.fields
    # print e2f1_peak_pos
    # raw_input()
    dist = min(abs(hsf1_mid - e2f1_peak_pos), abs(hsf1_mid - e2f1_peak_pos))
    if dist <= 500:
        if ":".join(feature.fields[0:3]) not in uniques:
            hsf1_nonHS_count += 1
            uniques.append(":".join(feature.fields[0:3]))


#repeat for HS
e2f1_hsf1_HS_regions_close = hsf1_HS_regions.sort().closest(e2f1_regions.sort(), D="a")
uniques = []
hsf1_HS_count = 0
for feature in e2f1_hsf1_HS_regions_close:
    if int(feature.fields[-1]) == -1: continue
    e2f1_start = int(feature.fields[-10])
    e2f1_peak = int(feature.fields[-2])
    e2f1_peak_pos = e2f1_start + e2f1_peak

    hsf1_mid = (int(feature.fields[1]) + int(feature.fields[2])) / 2.0

    # print feature.fields
    # print e2f1_peak_pos
    # raw_input()
    dist = min(abs(hsf1_mid - e2f1_peak_pos), abs(hsf1_mid - e2f1_peak_pos))
    if dist <= 500:
        if ":".join(feature.fields[0:3]) not in uniques:
            hsf1_HS_count += 1
            uniques.append(":".join(feature.fields[0:3]))

print "E2F1 TSS count: ", count
print "HSF1 E2F1 count: ", hsf1_count
print "HSF1 nonHS E2F1 count: ", hsf1_nonHS_count
print "HSF1 HS E2F1 count: ", hsf1_HS_count

print "HSF1 nonHS frac: ", float(hsf1_nonHS_count) / float(len(hsf1_nonHS_regions))
print "HSF1 HS frac: ", float(hsf1_HS_count) / float(len(hsf1_HS_regions))

