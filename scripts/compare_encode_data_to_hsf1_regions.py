import pandas as pd
import os
import pybedtools as pbt
import pickle
import chrom_utils
import mygene
import shutil
import math
import re

data_dir = os.path.join("..", "data_files")
encode_files_list_file = os.path.join(data_dir, "encode_files.pickle")

# hsf1_tss_hg19_file = os.path.join(data_dir, "hg19_TSS_filtered_hsf1_0.17_lifted.bed")
hg19_TSS_file = os.path.join(data_dir, "hg19_TSS.ranged.1.bed")
hsf1_regions_hg19_file = os.path.join(data_dir, "hsf1_0.17_merged_regions_hg19.refiltered.bed")

encode_files_list = pickle.load(open(encode_files_list_file))
encode_files_list.sort()

encode_files_data = [(os.path.basename(f),
                      f,
                      chrom_utils.get_cell_type_from_filename(os.path.basename(f)) + "-" +
                      chrom_utils.get_tf_type_from_filename(os.path.basename(f)) + "-" +
                      chrom_utils.get_exp_type_from_filename(os.path.basename(f))) for f in encode_files_list]

hg19_TSS = pbt.BedTool(hg19_TSS_file)
hsf1_regions = pbt.BedTool(hsf1_regions_hg19_file)
hsf1_regions.sort()

hsf1_hg19_TSS_fields = []
hsf1_hg19_tss_close = hg19_TSS.sort().closest(hsf1_regions, D="a")
chr_idxs = [i for i in range(len(hsf1_hg19_tss_close[0].fields)) if re.match(chrom_utils.chr_start_pattern, hsf1_hg19_tss_close[0].fields[i])]
hsf1_start_idx = chr_idxs[1]+1
for feature in hsf1_hg19_tss_close:
    if feature.fields[-1] == -1: continue
    hsf1_start = int(feature.fields[hsf1_start_idx])
    hsf1_end = int(feature.fields[hsf1_start_idx+1])
    hsf1_peak_pos = int(math.floor( (hsf1_start + hsf1_end)/2.0 ))

    tss_start = int(feature.fields[1])
    dist = hsf1_peak_pos - tss_start
    if abs(dist) <= 2000:
        hsf1_hg19_TSS_fields.append(feature.fields)

hsf1_tss_hg19 = pbt.BedTool(hsf1_hg19_TSS_fields).sort()


# assert len([x[0] for x in new_encode_files_data]) == len(set([x[0] for x in new_encode_files_data]))

close_df = pd.DataFrame(index=set([feature.name for feature in hsf1_tss_hg19]), columns=[x[0] for x in encode_files_data])

for fd in encode_files_data:
    bedf = fd[1]
    encode_bed = pbt.BedTool(bedf)
    close = hsf1_tss_hg19.sort().closest(encode_bed.sort(), D="a")
    for feature in close:
        if abs(int(feature[-1])) <= 1500 and int(feature[-1]) != -1:
            # print feature.name
            # print feature.name in close_df.index.values.tolist()
            close_df.set_value(index=feature.name, col=fd[0], value=int(feature[-1]))


writer = pd.ExcelWriter(os.path.join(data_dir, "encode_hsf1_tss_proximity.xlsx"))
close_df.to_excel(writer)
writer.save()

print "Analyzing HSF1 sites vs modencode"

#analyze the proximity to hsf1 peaks

intermediate_dir = os.path.join(data_dir, "hsf1_encode_region_comparisons")
if not os.path.isdir(intermediate_dir):
    os.mkdir(intermediate_dir)
else:
    shutil.rmtree(intermediate_dir)
    os.mkdir(intermediate_dir)



for fd in encode_files_data:
    bedf = fd[1]
    encode_bed = pbt.BedTool(bedf)
    close_regions = hsf1_regions.closest(encode_bed.sort(), D="a")
    hits = []
    for feature in close_regions:
        # print feature
        # raw_input()
        encode_peak = int(feature[-2])
        encode_start = int(feature[5])
        encode_peak_pos = encode_peak + encode_start

        hsf1_start = int(feature[1])
        hsf1_end = int(feature[2])
        hsf1_peak_pos = int(math.floor((hsf1_start - hsf1_end) / 2.0))

        dist = hsf1_peak_pos - encode_peak_pos

        if abs(dist) <= 500:
            fields = feature.fields
            fields[-1] = dist
            hits.append(fields)

        # if abs(int(feature[-1])) <= 500 and int(feature[-1]) != -1:
        #     hits.append(feature.fields)

    outfile = os.path.join(intermediate_dir, "hsf1_{e}_le_500.bed".format(e=fd[2]))
    pbt.BedTool(hits).sort().saveas(outfile)


comparison_files = [os.path.join(intermediate_dir, f) for f in os.listdir(intermediate_dir) if "bed" in f]

intermediate_regions_tss_dir = os.path.join(data_dir, "hsf1_encode_region_comparisons_tss")

if not os.path.isdir(intermediate_regions_tss_dir): os.mkdir(intermediate_regions_tss_dir)
else:
    shutil.rmtree(intermediate_regions_tss_dir)
    os.mkdir(intermediate_regions_tss_dir)

mg = mygene.MyGeneInfo()
query = mg.querymany([feature.name for feature in hsf1_tss_hg19], scopes=["refseq"], fields="symbol", returnall=True)
# print query["out"]
# raw_input()
syms = []
for q in query["out"]:
    syms.append((q["query"],q["symbol"]))



df = pd.DataFrame(index=[feature.name for feature in hsf1_tss_hg19], columns=["Symbol"]+ [os.path.basename(f).replace("_le_500.bed", "") for f in comparison_files])

for s in syms:
    df.set_value(index=s[0], col="Symbol", value=s[1])

print "Getting TSS comparisons"
for f in comparison_files:
    intermediate_bed = pbt.BedTool(f)
    col_name = os.path.basename(f).replace("_le_500.bed", "")
    close = hsf1_tss_hg19.sort().closest(intermediate_bed, D="a")
    outfeatures = []
    outfile = os.path.join(intermediate_regions_tss_dir, "hsf1_{cn}_tss.bed".format(cn=col_name))
    for feature in close:
        if abs(int(feature[-1])) <= 2000 and int(feature[-1]) != -1:
            outfeatures.append(feature.fields)

    outbed = pbt.BedTool(outfeatures).sort()
    outbed.saveas(outfile)

    for feature in outbed:
        df.set_value(index=feature.name, col=col_name, value=int(feature[-1]))



writer = pd.ExcelWriter(os.path.join(data_dir, "hsf1_encode_TSS_comparisons.xlsx"))
df.to_excel(writer)
writer.save()