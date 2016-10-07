import os
import shutil
import pybedtools as pbt
import pandas as pd
import mygene as mg

mg = mg.MyGeneInfo()

data_dir = os.path.join("..", "data_files")
hg19_tss_file = os.path.join(data_dir, "hg19_TSS.ranged.1.bed")
hsf1_cctf_region_files_dir = os.path.join(data_dir, "hsf1_cctf_region_comparisons")
cc_phase_ids_file = os.path.join(data_dir, "cell_cycle_osc_genes_hsf1_promoters.xlsx")

cobinding_files = [os.path.join(hsf1_cctf_region_files_dir,f) for f in os.listdir(hsf1_cctf_region_files_dir) if ".bed" in f]

hg19_tss = pbt.BedTool(hg19_tss_file)

intermediate_files_dir = os.path.join(data_dir, "hsf1_cobinding_setup_intermediate_files")
if not os.path.isdir(intermediate_files_dir):
    os.mkdir(intermediate_files_dir)
else:
    shutil.rmtree(intermediate_files_dir)
    os.mkdir(intermediate_files_dir)


cc_phase_ids = pd.read_excel(open(cc_phase_ids_file))
#maybe not the most elegant way, but let's go feature by feature in each file

print "Analyzing files"
output = {}
for cbf in cobinding_files:
    tf_type = os.path.basename(cbf).replace("hsf1_", "").replace("_le_500.bed", "")
    print "Analyzing data for: ", tf_type
    cb_bed = pbt.BedTool(cbf)
    close_5 = hg19_tss.sort().closest(cb_bed.sort(), D="a", k=5)
    good_fields = []
    uniques = []
    for feature in close_5:
        if abs(int(feature[-1])) <= 2000 and int(feature[-1]) != -1:
            if ":".join(feature.fields[0:3]) not in uniques:
                good_fields.append(feature.fields)
                uniques.append(":".join(feature.fields[0:3]))

    good_bed = pbt.BedTool(good_fields).sort()
    # print len(good_bed)
    # print good_bed[0].fields
    # raw_input()
    tss_names = [feature.name for feature in good_bed]
    curdf = pd.DataFrame( columns=["TSS Position","Symbol", "HSF1 Feature Name", "HSF1 - TSS", "{ct} - TSS".format(ct=tf_type), "HSF1 - {ct}".format(ct=tf_type), "Direction", "Phase"])

    query = mg.querymany(tss_names, scopes=["refseq"], fields="symbol", returnall=True)
    name_pairs = [(q["query"], q["symbol"]) for q in query["out"]]

    for feature in good_bed:
        tss_pos_name = ":".join(feature.fields[0:3])
        tss_name = feature.name
        hsf1_feature_name = feature.fields[9]
        direction = feature.fields[4]
        hsf1_midpoint = (int(feature.fields[7]) + int(feature.fields[8]))/2.0
        tss_pos = int(feature.fields[1])
        tf_peak = int(feature.fields[11]) + int(feature.fields[-3])
        hsf1_tss_dist = hsf1_midpoint - tss_pos if direction == "+" else tss_pos - hsf1_midpoint
        tf_tss_dist = tf_peak - tss_pos if direction == "+" else tss_pos - tf_peak
        hsf1_tf_dist = hsf1_midpoint - tf_peak if direction == "+" else tf_peak - hsf1_midpoint

        try:
            sym, = (x[1] for x in name_pairs if x[0] == tss_name)
        except ValueError:
            sym = None

        ccphase = cc_phase_ids[cc_phase_ids["HSF1 Feature"] == hsf1_feature_name]

        if ccphase.shape[0] == 1:
            ccphase = ccphase.iloc[0,:]
            phase = ccphase["Cell Cycle Phase"]
        elif ccphase.shape[0] < 1:
            phase = None
            # print tss_name
        else:
            ccphase = ccphase.iloc[0,:]
            phase = ccphase["Cell Cycle Phase"]

        if not sym in curdf["Symbol"].values.tolist():
            curdf.set_value(index=tss_name, col="TSS Position", value=tss_pos_name)
            curdf.set_value(index=tss_name, col="Symbol", value=sym)
            curdf.set_value(index=tss_name, col="HSF1 Feature Name", value=hsf1_feature_name)
            curdf.set_value(index=tss_name, col="HSF1 - TSS", value=hsf1_tss_dist)
            curdf.set_value(index=tss_name, col="{ct} - TSS".format(ct=tf_type), value=tf_tss_dist)
            curdf.set_value(index=tss_name, col="HSF1 - {ct}".format(ct=tf_type), value=hsf1_tf_dist)
            curdf.set_value(index=tss_name, col="Direction", value=direction)
            curdf.set_value(index=tss_name, col="Phase", value=phase)


    #save the data
    outname = os.path.join(intermediate_files_dir, "hsf1_{tf}_tss.bed".format(tf=tf_type))
    good_bed.saveas(outname)
    output[tf_type] = curdf


writer = pd.ExcelWriter(os.path.join(data_dir, "hsf1_cctf_position_analysis.xlsx"))
for k in output.keys():
    output[k].to_excel(writer, sheet_name=k)

writer.save()
