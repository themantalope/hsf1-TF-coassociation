import pandas as pd
import pybedtools as pbt
import os
import math

data_dir = os.path.join("..", "data_files")
cctf_dir = os.path.join(data_dir, "cc_tf_bedfiles")
e2f1_regions_file = os.path.join(cctf_dir, "E2F1", "E2F1.bed")
bmyb_regions_file = os.path.join(cctf_dir, "BMYB", "BMYB.bed")
foxm1_regions_file = os.path.join(cctf_dir, "FOXM1", "FOXM1.bed")
lin9_regions_file = os.path.join(cctf_dir, "LIN9", "LIN9.bed")
lin54_regions_file = os.path.join(cctf_dir, "LIN54", "LIN54.bed")
p130_regions_file = os.path.join(cctf_dir, "P130", "P130.bed")
hg19_tss_file = os.path.join(data_dir, "hg19_TSS.ranged.0.unique.bed")

e2f1_regions = pbt.BedTool(e2f1_regions_file)
bmyb_regions = pbt.BedTool(bmyb_regions_file)
foxm1_regions = pbt.BedTool(foxm1_regions_file)
lin9_regions = pbt.BedTool(lin9_regions_file)
lin54_regions = pbt.BedTool(lin54_regions_file)
p130_regions = pbt.BedTool(p130_regions_file)
hg19_tss = pbt.BedTool(hg19_tss_file)

#get the total number of unique TSS

uniques = []
for feature in hg19_tss:
    if ":".join(feature.fields[0:3]) not in uniques:
        uniques.append(":".join(feature.fields[0:3]))

num_unique_tss = len(uniques)

#get the number of e2f1 peaks near the TSS
e2f1_tss_close = hg19_tss.sort().closest(e2f1_regions.sort(), D="a")
uniques = []
e2f1_tss_count = 0
e2f1_tss_fields = []
for feature in e2f1_tss_close:
    if int(feature.fields[-1]) == -1: continue

    # print feature.fields
    # print feature.fields[6:16]
    # raw_input()
    e2f1_start = int(feature.fields[-10])
    e2f1_peak = int(feature.fields[-2])

    e2f1_peak_pos = e2f1_start + e2f1_peak
    dist = min( abs(int(feature.fields[1]) - e2f1_peak_pos), abs(int(feature.fields[2]) - e2f1_peak_pos) )
    if dist <= 2000:
        if ":".join(feature.fields[0:3]) not in uniques:
            e2f1_tss_count += 1
            uniques.append(":".join(feature.fields[0:3]))
            e2f1_tss_fields.append(feature.fields[6:16])


e2f1_tss_regions = pbt.BedTool(e2f1_tss_fields).sort()

#get the number near bmyb
uniques = []
bmyb_counts = 0
e2f1_bmyb_close = e2f1_tss_regions.closest(bmyb_regions.sort(), D="a")
print "BYMB"
print e2f1_bmyb_close[0].fields
print

for feature in e2f1_bmyb_close:
    e2f1_peak = int(feature.fields[9])
    e2f1_start = int(feature.fields[1])
    e2f1_peak_pos = e2f1_peak + e2f1_start

    bmyb_peak = int(feature.fields[-2])
    bmyb_start = int(feature.fields[-10])
    bmyb_peak_pos = bmyb_peak + bmyb_start

    dist = abs(e2f1_peak_pos - bmyb_peak_pos)
    if dist <= 500:
        if ":".join(feature.fields[0:3]) not in uniques:
            bmyb_counts += 1
            uniques.append(":".join(feature.fields[0:3]))


#we also need to get a count of the peaks near the tss for each TF

bmyb_hg19_tss_close = hg19_tss.sort().closest(bmyb_regions.sort(), D="a")
uniques = []
bmyb_tss_count = 0
for feature in bmyb_hg19_tss_close:
    if int(feature.fields[-1]) == -1: continue
    # print feature.fields
    # raw_input()
    bmyb_peak = int(feature.fields[-2])
    bmyb_start = int(feature.fields[-10])
    bmyb_peak_pos = bmyb_peak + bmyb_start

    dist = min( abs(bmyb_peak_pos - int(feature.fields[1])), abs(bmyb_peak_pos - int(feature.fields[2])) )
    if dist <= 2000:
        if ":".join(feature.fields[0:3]) not in uniques:
            bmyb_tss_count += 1
            uniques.append(":".join(feature.fields[0:3]))


#get the number near FOXM1

uniques = []
foxm1_counts = 0
e2f1_foxm1_close = e2f1_tss_regions.closest(foxm1_regions.sort(), D="a")
print "FOXM1"
print e2f1_foxm1_close[0].fields
print

for feature in e2f1_foxm1_close:
    e2f1_peak = int(feature.fields[9])
    e2f1_start = int(feature.fields[1])
    e2f1_peak_pos = e2f1_peak + e2f1_start

    foxm1_peak = int(feature.fields[-2])
    foxm1_start = int(feature.fields[-10])
    foxm1_peak_pos = foxm1_peak + foxm1_start

    dist = abs(e2f1_peak_pos - foxm1_peak_pos)
    if dist <= 500:
        if ":".join(feature.fields[0:3]) not in uniques:
            foxm1_counts += 1
            uniques.append(":".join(feature.fields[0:3]))


foxm1_hg19_tss_close = hg19_tss.sort().closest(foxm1_regions.sort(), D="a")
uniques = []
foxm1_tss_counts = 0
for feature in foxm1_hg19_tss_close:
    if int(feature.fields[-1]) == -1: continue
    foxm1_peak = int(feature.fields[-2])
    foxm1_start = int(feature.fields[-10])
    foxm1_peak_pos = foxm1_peak + foxm1_start

    dist = min( abs(foxm1_peak_pos - int(feature.fields[1])), abs(foxm1_peak_pos - int(feature.fields[2])) )
    if dist <= 2000:
        if ":".join(feature.fields[0:3]) not in uniques:
            foxm1_tss_counts += 1
            uniques.append(":".join(feature.fields[0:3]))


#get the number near LIN9

uniques = []
lin9_counts = 0
e2f1_lin9_close = e2f1_tss_regions.closest(lin9_regions.sort(), D="a")
print "LIN9"
print e2f1_lin9_close[0].fields
print

for feature in e2f1_lin9_close:
    e2f1_peak = int(feature.fields[9])
    e2f1_start = int(feature.fields[1])
    e2f1_peak_pos = e2f1_peak + e2f1_start

    # lin9_peak = int(feature.fields[-2])
    lin9_start = int(feature.fields[11])
    lin9_stop = int(feature.fields[12])
    lin9_peak_pos = int( math.floor( (lin9_start + lin9_stop) / 2.0) )

    dist = abs(lin9_peak_pos - e2f1_peak_pos)
    if dist <= 500:
        if ":".join(feature.fields[0:3]) not in uniques:
            lin9_counts +=1
            uniques.append(":".join(feature.fields[0:3]))


lin9_hg19_tss_close = hg19_tss.sort().closest(lin9_regions.sort(), D="a")
lin9_tss_count = 0
uniques = []

for feature in lin9_hg19_tss_close:
    if int(feature.fields[-1]) == -1: continue
    # print feature.fields
    lin9_start = int(feature.fields[7])
    lin9_stop = int(feature.fields[8])
    lin9_peak_pos = int(math.floor((lin9_start + lin9_stop) / 2.0))

    dist = min( abs(lin9_peak_pos - int(feature.fields[1])), abs(lin9_peak_pos - int(feature.fields[2])) )

    if dist <= 2000:
        if ":".join(feature.fields[0:3]) not in uniques:
            lin9_tss_count += 1
            uniques.append(":".join(feature.fields[0:3]))



#get the number near LIN54
uniques = []
lin54_counts = 0
e2f1_lin54_close = e2f1_tss_regions.closest(lin54_regions.sort(), D="a")
print "LIN54"
print e2f1_lin54_close[0].fields
print

for feature in e2f1_lin54_close:
    e2f1_peak = int(feature.fields[9])
    e2f1_start = int(feature.fields[1])
    e2f1_peak_pos = e2f1_peak + e2f1_start

    lin54_start = int(feature.fields[11])
    lin54_stop = int(feature.fields[12])
    lin54_peak_pos = int( math.floor( (lin54_start + lin54_stop) / 2.0) )

    dist = abs(e2f1_peak_pos - lin54_peak_pos)
    if dist <= 500:
        if ":".join(feature.fields[0:3]) not in uniques:
            lin54_counts += 1
            uniques.append(":".join(feature.fields[0:3]))


lin54_hg19_tss_close = hg19_tss.sort().closest(lin54_regions.sort(), D="a")
uniques = []
lin54_tss_count = 0
for feature in lin54_hg19_tss_close:
    if int(feature.fields[-1]) == -1: continue
    lin54_start = int(feature.fields[7])
    lin54_stop = int(feature.fields[8])
    lin54_peak_pos = int(math.floor((lin54_start + lin54_stop) / 2.0))

    dist = min(abs(lin54_peak_pos - int(feature.fields[1])), abs(lin54_peak_pos - int(feature.fields[2])))

    if dist <= 2000:
        if ":".join(feature.fields[0:3]) not in uniques:
            lin54_tss_count += 1
            uniques.append(":".join(feature.fields[0:3]))


#repeat for p130

uniques = []
p130_counts = 0
e2f1_p130_close = e2f1_tss_regions.closest(p130_regions.sort(), D="a")
print "P130"
print e2f1_p130_close[0].fields
print

for feature in e2f1_p130_close:
    e2f1_peak = int(feature.fields[9])
    e2f1_start = int(feature.fields[1])
    e2f1_peak_pos = e2f1_peak + e2f1_start

    p130_start = int(feature.fields[11])
    p130_stop = int(feature.fields[12])
    p130_peak_pos = int( math.floor( (p130_start + p130_stop) / 2.0) )

    dist = abs(e2f1_peak_pos - p130_peak_pos)
    if dist <= 500:
        if ":".join(feature.fields[0:3]) not in uniques:
            p130_counts += 1
            uniques.append(":".join(feature.fields[0:3]))


p130_hg19_tss_close = hg19_tss.sort().closest(p130_regions.sort(), D="a")
uniques = []
p130_tss_count = 0
for feature in p130_hg19_tss_close:
    if int(feature.fields[-1]) == -1: continue
    p130_start = int(feature.fields[7])
    p130_stop = int(feature.fields[8])
    p130_peak_pos = int(math.floor((p130_start + p130_stop) / 2.0))

    dist = min(abs(p130_peak_pos - int(feature.fields[1])), abs(p130_peak_pos - int(feature.fields[2])))

    if dist <= 2000:
        if ":".join(feature.fields[0:3]) not in uniques:
            p130_tss_count += 1
            uniques.append(":".join(feature.fields[0:3]))



#save the data

outdf = pd.DataFrame(index=["BMYB", "FOXM1", "LIN9", "LIN54"], columns=["N Sites Near E2F1", "N Sites Near TSS", "N E2F1 Sites", "N TSS Sites", "Ratio"])

outdf.set_value("BMYB", "N Sites Near E2F1", bmyb_counts)
outdf.set_value("BMYB", "N Sites Near TSS", bmyb_tss_count)
outdf.set_value("BMYB", "N E2F1 Sites", e2f1_tss_count)
outdf.set_value("BMYB", "N TSS Sites", num_unique_tss)
outdf.set_value("BMYB", "Ratio", (float(bmyb_counts)/e2f1_tss_count)/(float(bmyb_tss_count)/num_unique_tss) )

outdf.set_value("FOXM1", "N Sites Near E2F1", foxm1_counts)
outdf.set_value("FOXM1", "N Sites Near TSS", foxm1_tss_counts)
outdf.set_value("FOXM1", "N E2F1 Sites", e2f1_tss_count)
outdf.set_value("FOXM1", "N TSS Sites", num_unique_tss)
outdf.set_value("FOXM1", "Ratio", (float(foxm1_counts)/e2f1_tss_count)/(float(foxm1_tss_counts)/num_unique_tss) )

outdf.set_value("LIN9", "N Sites Near E2F1", lin9_counts)
outdf.set_value("LIN9", "N Sites Near TSS", lin9_tss_count)
outdf.set_value("LIN9", "N E2F1 Sites", e2f1_tss_count)
outdf.set_value("LIN9", "N TSS Sites", num_unique_tss)
outdf.set_value("LIN9", "Ratio", (float(lin9_counts)/e2f1_tss_count)/(float(lin9_tss_count)/num_unique_tss) )

outdf.set_value("LIN54", "N Sites Near E2F1", lin54_counts)
outdf.set_value("LIN54", "N Sites Near TSS", lin54_tss_count)
outdf.set_value("LIN54", "N E2F1 Sites", e2f1_tss_count)
outdf.set_value("LIN54", "N TSS Sites", num_unique_tss)
outdf.set_value("LIN54", "Ratio", (float(lin54_counts)/e2f1_tss_count)/(float(lin54_tss_count)/num_unique_tss) )

outdf.set_value("P130", "N Sites Near E2F1", p130_counts)
outdf.set_value("P130", "N Sites Near TSS", p130_tss_count)
outdf.set_value("P130", "N E2F1 Sites", e2f1_tss_count)
outdf.set_value("P130", "N TSS Sites", num_unique_tss)
outdf.set_value("P130", "Ratio", (float(p130_counts)/e2f1_tss_count)/(float(p130_tss_count)/num_unique_tss) )

writer = pd.ExcelWriter(os.path.join(data_dir, "e2f1_other_cctf_overlap.xlsx"))
outdf.to_excel(writer)
writer.save()

