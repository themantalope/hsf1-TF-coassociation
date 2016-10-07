import pandas as pd
import json
import os
import pybedtools as pbt
import mygene

mg = mygene.MyGeneInfo()

data_dir = os.path.join("..", "data_files")
oscillations_file = os.path.join(data_dir, "dominguez_wang_periodic_expression_modified.xlsx")
hsf_regions_file = os.path.join(data_dir, "hsf1_0.17_merged_regions_hg19.refiltered.bed")
tss_file = os.path.join(data_dir, "hg19_TSS.ranged.1.bed")

tss = pbt.BedTool(tss_file)
hsf1 = pbt.BedTool(hsf_regions_file)
close = tss.sort().closest(hsf1.sort(), D="a")
good_features = []
for feature in close:
    if int(feature.fields[-1]) == -1:continue
    if abs(int(feature.fields[-1])) <= 2000:
        good_features.append(feature.fields)

len(good_features)
goodbed = pbt.BedTool(good_features)
names = [feature.name for feature in goodbed]
query = mg.querymany(names, scopes=["refseq"], fields="ensembl.gene", returnall=True)
print query
ensgs = []
for q in query["out"]:
    if "ensembl" not in q: continue
    if isinstance(q["ensembl"], list):
        ensgs += [g["gene"] for g in q["ensembl"]]
    else:
        ensgs.append(q["ensembl"]["gene"])

print ensgs
# raw_input()


oscillations = pd.read_excel(open(oscillations_file))

good_oscs = oscillations[oscillations["ENSG"].isin(ensgs)]


writer = pd.ExcelWriter(os.path.join(data_dir, "hsf1_DW_reg_oscillatory_genes.xlsx"))
good_oscs.to_excel(writer)


writer.save()