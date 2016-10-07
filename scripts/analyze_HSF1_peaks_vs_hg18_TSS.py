import os
import pybedtools as pbt
import chrom_utils
import pandas as pd
import shutil
import datetime

now = datetime.datetime.now().isoformat()

data_dir = os.path.join("..", "data_files")
mendillo_cell_type_files = os.path.join(data_dir, "mendillo_good_cell_lines.txt")
mendillo_dir = os.path.join(data_dir, "mendillo_merged")
hg18_tss_file = os.path.join(data_dir, "hg18_TSS.ranged.0.bed")



mendillo_hsf1_bed_files = [os.path.join(mendillo_dir, f) for f in os.listdir(mendillo_dir) if "merged" not in f]

cell_types = []
for f in mendillo_hsf1_bed_files:
    cell_types.append(chrom_utils.get_cell_type_from_filename(os.path.basename(f)))

for i, _ in enumerate(cell_types):
    cell_types[i] = cell_types[i].strip()


# print mendillo_hsf1_bed_files
hg18_tss = pbt.BedTool(hg18_tss_file)


#we need to check the hg18tss file to make sure that we don't count redundant TSSs

good_hg18_tss = []

hg18_tss.sort()


hg18_tss_fields = [feature.fields for feature in hg18_tss]

i = 0
while i < len(hg18_tss_fields) - 1:
    cur_tss = hg18_tss_fields[i]
    next_tss = hg18_tss_fields[i+1]

    if "_" in cur_tss[0]:
        i+=1
        continue
    elif "NR" in cur_tss[3]:
        i+=1
        continue

    dist = abs(int(cur_tss[2]) - int(next_tss[1]))

    if dist == 0:
        good_hg18_tss.append(cur_tss)
        while dist == 0:
            i+=1
            dist = abs(int(hg18_tss_fields[i][2]) - int(hg18_tss_fields[i+1][1]))
    elif dist <= 2000:
        good_hg18_tss.append(cur_tss)
        i += 1
    else:
        good_hg18_tss.append(cur_tss)

    i+=1
    # print i


hg18_tss_unique = pbt.BedTool(good_hg18_tss)
# print len(hg18_tss_unique)
# print [feature.fields for feature in hg18_tss_unique[0:10]]
temp = os.path.join(data_dir, "temp.bed")
hg18_tss_unique.saveas(temp)
# raw_input()

hg18_tss = hg18_tss_unique

outfile = os.path.join(data_dir, "hg18_TSS.ranged.0.uniquecoding.bed")
hg18_tss.saveas(outfile)
hg18_tss.sort()

if not False:

    intermediate_dir = os.path.join(data_dir, "intermediate_files")

    if not os.path.isdir(intermediate_dir):
        os.mkdir(intermediate_dir)
    else:
        shutil.rmtree(intermediate_dir)
        os.mkdir(intermediate_dir)

    genes_with_peak = set()

    for f in mendillo_hsf1_bed_files:
        cell = chrom_utils.get_cell_type_from_filename(os.path.basename(f))
        hsfbed = pbt.BedTool(f).sort()
        #reset the peak region to be +/- 500 bp from the peak

        # bed_fields = [feature.fields for feature in hsfbed]
        # # print "before: ", bed_fields[0]
        # for fields in bed_fields:
        #     peak = int(fields[-1])
        #     peak_pos = peak + int(fields[1])
        #     fields[1] = str(peak_pos - 200)
        #     fields[2] = str(peak_pos + 200)
        #     fields[-1] = str(200)
        # # print "after: ", bed_fields[0]
        #
        #
        # hsfbed = pbt.BedTool(bed_fields).sort()
        # hsfbed = hsfbed.merge(i=hsfbed.fn, d=200).sort()
        # raw_input()

        #find the closest features across the two files
        closest = hg18_tss.sort().closest(hsfbed, D="a")
        chr_idxs = [i for i in range(len(closest[0].fields)) if "chr" in closest[0].fields[i]]
        hsf1_start_idx = chr_idxs[1]+1
        under_2000_kb = []
        for feature in closest:
            if int(feature.fields[-1]) == -1: continue
            hsf1_peak = int(feature.fields[-2])
            hsf1_peak_pos = int(feature.fields[hsf1_start_idx]) + hsf1_peak

            tss_start = int(feature.fields[1])
            dist = abs( hsf1_peak_pos - tss_start )


            if dist <= 2000 and int(feature[-1]) != -1:
                under_2000_kb.append(feature)
                genes_with_peak.add(feature.name)

        under_2000_kb_bed = pbt.BedTool(under_2000_kb)

        intermediate_name = os.path.join(intermediate_dir, "{c}_HSF1_hg18_TSS_under2kb.bed".format(c=cell))
        under_2000_kb_bed.saveas(intermediate_name)


    intermediate_files = [os.path.join(intermediate_dir,f) for f in os.listdir(intermediate_dir) if os.path.isfile(os.path.join(intermediate_dir,f))]
    # print intermediate_files
    present_df = pd.DataFrame(index=list(genes_with_peak), columns=cell_types)
    distance_df = pd.DataFrame(index=list(genes_with_peak), columns=cell_types)

    for f in intermediate_files:
        # print f
        cell_type = None
        for ct in cell_types:
            if ct in f: cell_type = ct; break

        if cell_type is None:
            raise ValueError("Couldn't find cell type: {f}".format(f=f))


        hsf1_peaks_TSS = pbt.BedTool(f)
        genes = set()
        for feature in hsf1_peaks_TSS:
            # print feature
            # print feature.name
            # print feature[-1]
            # raw_input()
            genes.add(feature.name)
            distance_df.set_value(index=feature.name, col=cell_type, value=int(feature[-1]))
            present_df.set_value(index=feature.name, col=cell_type, value=1.0)

        # print genes
    writer = pd.ExcelWriter(os.path.join(data_dir, "hsf1_hg18_TSS_peaks_2000kb.xlsx"))
    present_df["fraction"] = present_df.apply(lambda x: x[:].sum()/float(len(x)), axis=1)
    distance_df.to_excel(writer, sheet_name="distance")
    present_df.to_excel(writer, sheet_name="present")
    writer.save()

else:
    data = pd.read_excel(os.path.join(data_dir, "hsf1_hg18_TSS_peaks_2000kb.xlsx"), sheetname=None)
    present_df = data["present"]
    distance_df = data["distance"]
    intermediate_dir = os.path.join(data_dir, "intermediate_files")
    intermediate_files = [os.path.join(intermediate_dir,f) for f in os.listdir(intermediate_dir) if os.path.isfile(os.path.join(intermediate_dir,f))]


#now for each level of threshold we need to create a tss file that only contains TSSs that are found within some % of
#the datasets

cutoffs = list(set(present_df["fraction"].values.tolist()))
# print cutoffs

merged_hsf1_regions_dir = os.path.join(data_dir, "hsf1_merged_region_files")
if not os.path.isdir(merged_hsf1_regions_dir):
    os.mkdir(merged_hsf1_regions_dir)

for cutoff in cutoffs:
    genes_in_cutoff = present_df[present_df["fraction"] >= cutoff].index.tolist()
    genes_not_in_cutoff = present_df[~(present_df["fraction"]>=cutoff)].index.tolist()
    features_for_cutoff = []
    features_for_not_cutoff = []
    for f in intermediate_files:
        hsf1_peaks_TSS = pbt.BedTool(f)
        for feature in hsf1_peaks_TSS:
            if feature.name in genes_in_cutoff:
                features_for_cutoff.append(feature[6:10])
            if feature.name in genes_not_in_cutoff:
                features_for_not_cutoff.append(feature[6:10])

    cutoff_features = pbt.BedTool(features_for_cutoff).sort().remove_invalid()
    not_cutoff_features = pbt.BedTool(features_for_not_cutoff).sort().remove_invalid()

    cutoff_merged_features = cutoff_features.merge(i=cutoff_features.fn, d=200).remove_invalid().sort()
    not_cutoff_merged_features = not_cutoff_features.merge(i=not_cutoff_features.fn, d=200).remove_invalid().sort()


    print "cutoff: ", cutoff
    print "total features across datasets: ", len(features_for_cutoff)
    print "total after merging: ", len(cutoff_merged_features)
    merged_name = os.path.join(merged_hsf1_regions_dir, "hsf1_{c}_merged_regions.bed".format(c=cutoff))
    not_cutoff_merged_name = os.path.join(merged_hsf1_regions_dir, "hsf1_{c}_merged_regions.negatives.bed".format(c=cutoff))

    cutoff_merged_features = [feature.fields for feature in cutoff_merged_features]
    for i, feature in enumerate(cutoff_merged_features):
        feature.append("feature{n}|{c}-{start}-{stop}".format(n=i, c=feature[0], start=feature[1],stop=feature[2]))

    cutoff_merged_features = pbt.BedTool(cutoff_merged_features)

    not_cutoff_merged_features = [feature.fields for feature in not_cutoff_merged_features]
    for i, feature in enumerate(not_cutoff_merged_features):
        feature.append("feature{n}|{c}-start-{stop}".format(n=i,c=feature[0],start=feature[1],stop=feature[2]))

    not_cutoff_merged_features = pbt.BedTool(not_cutoff_merged_features)

    cutoff_merged_features.saveas(merged_name)
    not_cutoff_merged_features.saveas(not_cutoff_merged_name)

metafile = os.path.join(merged_hsf1_regions_dir, "log.txt")
with open(metafile, "w") as f:
    f.write("Data generated at : {n}".format(n=now))

f.close()