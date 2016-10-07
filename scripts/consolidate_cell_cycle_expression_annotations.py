import os
import pandas as pd
import mygene as mg
import pybedtools as pbt
import pickle

mg= mg.MyGeneInfo()

data_dir = os.path.join("..", "data_files")
gw_expression_file = os.path.join(data_dir, "grant_whitfield_cell_cycle_expression.xlsx")
gw_annotations_file = os.path.join(data_dir, "grant_whitfield_probe_annotations.xlsx")
dw_expression_file = os.path.join(data_dir, "dominguez_wang_periodic_expression_modified.xlsx")
hg19_tss_file = os.path.join(data_dir, "hg19_TSS.ranged.1.bed")

refseq_ids = [feature.name for feature in pbt.BedTool(hg19_tss_file) if "_" not in feature[0]]

#first do the gw expression

gw_expression_dict = pd.read_excel(open(gw_expression_file), sheetname=None)
gw_expression = gw_expression_dict["original"]
gw_annotations = gw_expression_dict["ID_Mapping"]

rs_annotations = gw_annotations[gw_annotations["RefSeqAccession"].isin(refseq_ids)]
rs_expression = gw_expression[gw_expression["CLID"].isin(rs_annotations["GW_ProbeID"].values)]

gw_rs_phases = []

for row in rs_expression.iterrows():
    idx = row[0]
    rowdata = row[1]
    gwid = rowdata["CLID"]
    phase = rowdata["Phase"]
    rs_id = rs_annotations.loc[rs_annotations["GW_ProbeID"]==gwid, "RefSeqAccession"]
    gw_rs_phases.append((rs_id.iloc[0], phase))


#get the dw data

pfile = os.path.join(data_dir, "temp_refseq_query.pickle")
if os.path.isfile(pfile):
    query = pickle.load(open(pfile))
else:
    query = mg.querymany(refseq_ids, scopes=["refseq"], fields="ensembl.gene", returnall=True)
    pickle.dump(query, open(pfile, "wb"))



# print query["out"]
# raw_input()

ensembl_genes = []

for q in query["out"]:
    if "ensembl" not in q: continue
    if isinstance(q["ensembl"], list):
        ensembl_genes.append((q["query"], [x["gene"] for x in q["ensembl"]]))
    else:
        ensembl_genes.append((q["query"], [q["ensembl"]["gene"]]))

all_ensembl = [item for sublist in ensembl_genes for item in [x for x in sublist[1]] ]
dw_expression = pd.read_excel(open(dw_expression_file))

dw_expression = dw_expression[dw_expression["ENSG"].isin(all_ensembl)]


# print ensembl_genes

dw_rs_phase = []
for row in dw_expression.iterrows():
    idx = row[0]
    rowdata = row[1]
    ensembl_id = rowdata["ENSG"]
    phase = rowdata["Cell Cycle Label"]
    rs_ids = [x[0] for x in ensembl_genes if ensembl_id in x[1]]
    # dw_rs_phase.append(rs_id, phase)
    for rsi in rs_ids:
        dw_rs_phase.append((rsi, phase))

outdf = pd.DataFrame(index=refseq_ids, columns=["GW Phase", "DW Phase", "Genomic Position"])

for gwr in gw_rs_phases:
    outdf.set_value(gwr[0],"GW Phase", value=gwr[1])

for dwr in dw_rs_phase:
    outdf.set_value(dwr[0], "DW Phase", value=dwr[1])

for feature in pbt.BedTool(hg19_tss_file):
    if feature.name in outdf.index.values:
        outdf.set_value(feature.name, "Genomic Position", value=":".join(feature.fields[0:3]))

outname = os.path.join(data_dir, "cell_cycle_phase_ids.xlsx")
writer = pd.ExcelWriter(outname)
outdf.to_excel(writer)
writer.save()