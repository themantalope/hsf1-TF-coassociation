import pandas as pd
import os
import pybedtools as pbt
import mygene
import math

data_dir = os.path.join("..", "data_files")
hg19_tss_file = os.path.join(data_dir, "hg19_TSS.ranged.1.bed")
hsf1_nonHS_file = os.path.join(data_dir , "hsf1_nonHS_exclusive.bed")
gw_annotations_file = os.path.join(data_dir, "grant_whitfield_probe_annotations.xlsx")
gw_data_file = os.path.join(data_dir, "grant_whitfield_cell_cycle_expression.xlsx")

hg19_tss = pbt.BedTool(hg19_tss_file)
hsf1_nonHS = pbt.BedTool(hsf1_nonHS_file)

close = hsf1_nonHS.sort().closest(hg19_tss.sort(), D="a")
close_fields = []
for feature in close:
    hsf1_peak = int( math.floor((int(feature.fields[1]) + int(feature.fields[2]))/2.0 ) )
    tss = int(feature.fields[5])

    dist = abs(tss - hsf1_peak)

    if dist <= 2000:
        close_fields.append(feature.fields)


names = [f[7] for f in close_fields]

# mg = mygene.MyGeneInfo()
#
# query = mg.querymany(names, scopes=["refseq"], fields="ensembl.gene, ensembl.transcript", returnall=True)
# # print query
# ensgs = []
# ensts = []
#
# for q in query["out"]:
#     if "ensembl" not in q: continue
#
#     if isinstance(q["ensembl"], list):
#         ensgs += [x["gene"] for x in q["ensembl"]]
#     else:
#         ensgs.append(q["ensembl"]["gene"])
#
#     if isinstance(q["ensembl"], list):
#         ensts += [x for t in q["ensembl"] for x in t["transcript"]]
#     else:
#         ensts += [x for x in q["ensembl"]["transcript"]]
#
# #see how many of these transcripts we can pull out from the annotations file
#
print "Loading data"
gw_annotations = pd.read_excel(open(gw_annotations_file))
gw_data = pd.read_excel(open(gw_data_file), sheetname=None)
gw_data = gw_data["original"]
print "Data loaded"
#
# print ensts+names
#
# gw_annotations_cut = gw_annotations[gw_annotations["ENSEMBL_ID"].isin(ensts+names)]
#
#
# print len(gw_data[gw_data["CLID"].isin(gw_annotations_cut["GW CC Phase ID"])])

#let's try making a quick annotation bed file
if not os.path.exists(os.path.join(data_dir, "gwatemp.bed")):
    bedfields = []
    for row in gw_annotations.iterrows():
        pos = row[1]["CHROMOSOMAL_LOCATION"]
        if pd.isnull(pos) or "unmapped" in pos: continue
        pos = pos.replace(":", "-")
        chr = pos.split("-")[0]
        start = int(pos.split("-")[1])
        end = int(pos.split("-")[2])
        if start == end: continue
        if start > end: continue
        name = row[1]["GW CC Phase ID"]
        start = str( min(int(start), int(end)) )
        end  = str( max(int(start), int(end)) )
        bedfields.append([chr, start, end, name])

    temp_bed = pbt.BedTool(bedfields)
    temp_bed.remove_invalid().sort().saveas(os.path.join(data_dir, "gwatemp.bed"))
else:
    temp_bed = pbt.BedTool(os.path.join(data_dir, "gwatemp.bed"))
    print len(temp_bed)
    temp_bed.remove_invalid()


close = hsf1_nonHS.sort().closest(temp_bed.sort().remove_invalid(), D="a")
close_probes = []
print close[0].fields
print len(close)
for feature in close:
    hsf1_mid = int(math.floor(int(feature.fields[1]) + int(feature.fields[2])/2.0))
    start = int(feature.fields[5])
    end = int(feature.fields[6])

    dist = abs( min( hsf1_mid - start, hsf1_mid - end ) )
    if dist <= 2000:
        close_probes.append(feature.fields)

print len(close_probes)
