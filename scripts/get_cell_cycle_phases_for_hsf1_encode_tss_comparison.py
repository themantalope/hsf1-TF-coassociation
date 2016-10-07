import pandas as pd
import os
import mygene

mg = mygene.MyGeneInfo()

data_dir = os.path.join("..", "data_files")

dw_expression_file = os.path.join(data_dir, "dominguez_wang_periodic_expression_modified.xlsx")
gw_expression_file = os.path.join(data_dir, "grant_whitfield_cell_cycle_expression.xlsx")
hsf1_encode_TSS_file = os.path.join(data_dir, "hsf1_encode_TSS_comparisons.xlsx")

print "Loading data"
gw_expression_dict = pd.read_excel(open(gw_expression_file), sheetname=None)
dw_expression = pd.read_excel(open(dw_expression_file))
hsf1_encode_TSS = pd.read_excel(open(hsf1_encode_TSS_file))

gw_id_mapping = gw_expression_dict["ID_Mapping"]
gw_expression = gw_expression_dict["original"]

tss_ids = hsf1_encode_TSS.index.tolist()

gw_ids = gw_id_mapping[gw_id_mapping["RefSeqAccession"].isin(tss_ids)]
gw_ids = gw_ids[["RefSeqAccession","GW_ProbeID"]]

gw_expression = gw_expression[gw_expression["CLID"].isin(gw_ids["GW_ProbeID"])]

for row in gw_ids.iterrows():
    idx = row[0]
    rowdata = row[1]
    refseq = rowdata["RefSeqAccession"]
    gw_id = rowdata["GW_ProbeID"]
    # print rowdata
    try:
        phase = gw_expression.loc[gw_expression["CLID"] == gw_id, "Phase"]
        # print phase
        # print phase.iloc[0]
        hsf1_encode_TSS.set_value(refseq, "GW Phase", phase.iloc[0])
    except KeyError:
        continue


#get data for DW
dw_query = mg.querymany(tss_ids, scopes=["refseq"], fields="ensembl.gene", returnall=True)

# print dw_query
# raw_input()

ensembl_genes = []
for q in dw_query["out"]:
    # print q
    if "ensembl" not in q: continue
    if isinstance(q["ensembl"], list):
        ensembl_genes.append(( q["query"], [x["gene"] for x in q["ensembl"]] ))
    else:
        ensembl_genes.append((q["query"], [q["ensembl"]["gene"]]))

all_ensembl_genes = []
for eg in ensembl_genes:
    all_ensembl_genes += eg[1]

dw_expression = dw_expression[dw_expression["ENSG"].isin(all_ensembl_genes)]

for row in dw_expression.iterrows():
    idx = row[0]
    rowdata = row[1]
    phase = rowdata["Cell Cycle Label"]
    ensg = rowdata["ENSG"]
    for eg in ensembl_genes:
        if ensg in eg[1]:
            hsf1_encode_TSS.set_value(eg[0], "DW Phase", phase)
            break

writer = pd.ExcelWriter(hsf1_encode_TSS_file)
hsf1_encode_TSS.to_excel(writer)
writer.save()
