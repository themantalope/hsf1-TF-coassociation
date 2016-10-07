import pandas as pd
import os
import pybedtools as pbt
import matplotlib.pyplot as plt

data_dir = os.path.join("..", "data_files")
grant_whitfield_file = os.path.join(data_dir, "grant_whitfield_cell_cycle_expression.xlsx")
hsf1_target_genes_file = os.path.join(data_dir,"hsf1_0.17_merged_regions_hg19.refiltered.bed")
probe_annotations_file = os.path.join(data_dir, "grant_whitfield_probe_annotations.xlsx")
hg18_TSS_bed_file = os.path.join(data_dir, "hg19_TSS.ranged.1.bed")

grant_whitfield_expression = pd.read_excel(open(grant_whitfield_file), sheetname=None)


# print grant_whitfield_expression

# hsf1_target_genes = pd.read_excel(open(hsf1_target_genes_file), sheetname=None)
probe_ids = pd.read_excel(open(probe_annotations_file))


# gene_list_sheet = hsf1_target_genes["present"]
# gene_list_sheet = gene_list_sheet[gene_list_sheet["fraction"] >= 0.1]

# refseq_genes = list(set(gene_list_sheet.index.tolist()))
hsf1_target_genes = pbt.BedTool(hsf1_target_genes_file)
hg18_TSS = pbt.BedTool(hg18_TSS_bed_file)
close = hg18_TSS.sort().closest(hsf1_target_genes.sort(), D="a").sort()
cout = []
for feature in close:
    if abs(int(feature[-1])) <= 1500 and int(feature[-1]) != -1:
        cout.append(feature.fields)

close_bed = pbt.BedTool(cout)

outfile = os.path.join(data_dir, "hg18_TSS_filtered_by_hsf_0.15.bed")
pbt.BedTool([[feature[0], feature[1], str(int(feature[2]) + 1), feature[3]] for feature in close_bed]).sort().saveas(outfile)

refseq_genes = []
for feature in close_bed:
    refseq_genes.append(feature[3])



print "number of refseq ids: ", len(refseq_genes)

probes_in_set = probe_ids[probe_ids["REFSEQ"].isin(refseq_genes)]

print "number of found probes: ", len(probes_in_set)
probe_names = probes_in_set["NAME"].values.tolist()

mapping = grant_whitfield_expression["ID_Mapping"]

hsf1_mapped = mapping[mapping["ProbeID"].isin(probe_names)]

# print hsf1_mapped["GeneSymbol"].values.tolist()

# mg = mygene.MyGeneInfo()
#
# query = mg.querymany(refseq_genes, scopes=["refseq"], fields="entrezgene", returnall=True)
# entrez = set()
#
# for q in query["out"]:
#
#     entrez.add(q["entrezgene"])

# print len(entrez)

# print entrez

# agilent_id = grant_whitfield_expression["ID_Mapping"]
# agilent_rows = agilent_id[agilent_id["EntrezGeneID"].isin([str(e) for e in entrez])]
#
# print "agilent rows found: ", len(agilent_rows)

hsf1_cc1_data = grant_whitfield_expression["CC1"]
hsf1_cc1_data = hsf1_cc1_data[hsf1_cc1_data["Hours"].isin(hsf1_mapped["GW_ProbeID"].values)]

original_rows = grant_whitfield_expression["original"]
original_rows = original_rows[original_rows["CLID"].isin(hsf1_mapped["GW_ProbeID"].values)]
# print original_rows["Phase"].values


writer = pd.ExcelWriter(os.path.join(data_dir, "quick_look_hsf1_targets.xlsx"))
original_rows.to_excel(writer)
writer.save()
# for row in hsf1_cc1_data.iterrows():
#     rd = row[1]
#     # print rd
#     # raw_input()
#     plt.plot(range(0,len(rd.iloc[1:])), rd.iloc[1:].values)
#     plt.show()