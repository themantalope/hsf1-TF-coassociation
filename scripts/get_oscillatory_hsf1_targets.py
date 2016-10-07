import pandas as pd
import os
import pybedtools as pbt
import mygene

data_dir = os.path.join("..", "data_files")
hg19_TSS_file = os.path.join(data_dir, "hg19_TSS.ranged.1.bed")
cell_cycle_phase_ids_file = os.path.join(data_dir, "cell_cycle_phase_ids.xlsx")
hsf1_regions_file = os.path.join(data_dir, "hsf1_0.17_merged_regions_hg19.refiltered.bed")

gw_annotations_file = os.path.join(data_dir, "grant_whitfield_probe_annotations.xlsx")
gw_data_file = os.path.join(data_dir, "grant_whitfield_cell_cycle_expression.xlsx")
gw_annotations = pd.read_excel(open(gw_annotations_file))
gw_data = pd.read_excel(open(gw_data_file), sheetname=None)
gw_data_orig = gw_data["original"]


gw_annotations["GW CC Phase ID"] = gw_annotations.apply(lambda x: "AGI_HUM1_OLIGO_" + x["NAME"] if not pd.isnull(x["NAME"]) else "", axis=1)
quickwriter = pd.ExcelWriter(gw_annotations_file)
gw_annotations.to_excel(quickwriter)
quickwriter.save()

hg19_TSS = pbt.BedTool(hg19_TSS_file)
hsf1_regions = pbt.BedTool(hsf1_regions_file)
cell_cycle_phase_ids = pd.read_excel(open(cell_cycle_phase_ids_file))

close = hg19_TSS.sort().closest(hsf1_regions.sort(), D="a")

close_features = []

for feature in close:
    if abs(int(feature[-1])) <= 2000 and int(feature[-1]) != -1:
        close_features.append(feature.fields)

close_hsf1 = pbt.BedTool(close_features).sort()

osc_df = pd.DataFrame(index=[feature.name for feature in close_hsf1], columns=["HSF1 Feature","Symbol", "Cell Cycle Phase", "Study"])

mg = mygene.MyGeneInfo()

query = mg.querymany([feature.name for feature in close_hsf1], scopes=["refseq"], fields="symbol", returnall=True)


gw_annotations_rows = gw_annotations[gw_annotations["REFSEQ"].isin([feature.name for feature in close_hsf1])]
gw_annos_gw = gw_annotations_rows["GW CC Phase ID"].values.tolist()

gwhits = gw_data_orig[gw_data_orig["CLID"].isin(gw_annos_gw)]


print "number of rows of gw annotations found: ", gw_annotations_rows.shape
print "number found in study: ", gwhits.shape
raw_input()

id_pairs = []
for q in query["out"]:
    if "symbol" not in q: continue
    id_pairs.append((q["query"], q["symbol"]))


for feature in close_hsf1:
    tss_name = feature.name
    feature_loc = ":".join(feature.fields[0:3])
    #search for the refseq id first
    tss_row = cell_cycle_phase_ids.loc[tss_name, :]
    # print tss_row.shape
    if len(tss_row.shape) > 1:
        tss_row = tss_row.iloc[0, :]
    # print tss_row
    loc_rows = cell_cycle_phase_ids.loc[cell_cycle_phase_ids["Genomic Position"] == feature_loc, :]
    loc_rows = loc_rows[~(pd.isnull(loc_rows[["GW Phase", "DW Phase"]])) ]
    # print loc_rows
    # raw_input()
    phase = None
    study = None

    if not pd.isnull(tss_row["DW Phase"]):
        phase = tss_row["DW Phase"]
        study = "DW"

    if not pd.isnull(tss_row["GW Phase"]):
        phase = tss_row["GW Phase"]
        study = "GW"

    if not pd.isnull(loc_rows[["GW Phase", "DW Phase"]]).values.all() and phase is None:
        gw = loc_rows[loc_rows["GW Phase"].notnull()]["GW Phase"]
        dw = loc_rows[loc_rows["DW Phase"].notnull()]["DW Phase"]

        # print gw
        # print dw
        if not dw.empty:
            # print gw
            phase = dw.iloc[0]
            study = "DW"

        if not gw.empty:
            # print dw
            # print dw.iloc[0]
            phase = gw.iloc[0]
            study = "GW"

    try:
        sym, = (x[1] for x in id_pairs if x[0] == tss_name)
    except ValueError:
        sym = None

    if sym not in osc_df["Symbol"].values.tolist():
        osc_df.set_value(index=tss_name, col="Symbol", value=sym)
        osc_df.set_value(index=tss_name, col="HSF1 Feature", value=feature.fields[9])
        osc_df.set_value(index=tss_name, col="Cell Cycle Phase", value=phase)
        osc_df.set_value(index=tss_name, col="Study", value=study)

writer = pd.ExcelWriter(os.path.join(data_dir, "cell_cycle_osc_genes_hsf1_promoters.xlsx"))
osc_df.to_excel(writer)
writer.save()






