import os
import pandas as pd
import chrom_utils
import pickle
import pybedtools as pbt
import progressbar
import math
# import pprint

data_dir = os.path.join("..", "data_files")
hsf1_encode_comparison_dir = os.path.join(data_dir, "hsf1_encode_region_comparisons")
hsf1_regions_file = os.path.join(data_dir, "hsf1_0.17_merged_regions_hg19.refiltered.bed")
hg19_TSS_file = os.path.join(data_dir, "hg19_TSS.ranged.1.bed")

hg19_TSS = pbt.BedTool(hg19_TSS_file)

encode_files_list_file = os.path.join(data_dir, "encode_files.pickle")
encode_files_list = pickle.load(open(encode_files_list_file))

# pprint.pprint(encode_files_list)
# raw_input()

#we need to see how frequently a peak from the encode data is found near a TSS

encode_files_data = [(chrom_utils.get_cell_type_from_filename(os.path.basename(f)),
                      chrom_utils.get_tf_type_from_filename(os.path.basename(f)),
                      chrom_utils.get_exp_type_from_filename(os.path.basename(f)),
                      f) for f in encode_files_list]

tss_count_df = pd.DataFrame(index=["-".join(x[0:3]) for x in encode_files_data], columns=["TSS Count", "Total Binding Sites", "Genomic Fraction"])

pbar = progressbar.ProgressBar(maxval=len(encode_files_data), widgets=["Getting TSS counts for ENCODE data", progressbar.Bar(), progressbar.Percentage()]).start()

for i, efd in enumerate(encode_files_data):
    efile = efd[3]
    bt = pbt.BedTool(efile)
    #reset the regions in bt so that it compares peak position to TSS regions
    newbtfeatures = []
    for feature in bt:
        peak_pos = int(feature[-1]) + int(feature[1])
        fields = feature.fields
        fields[1] = str(peak_pos)
        fields[2] = str(peak_pos)
        newbtfeatures.append(fields)

    bt = pbt.BedTool(newbtfeatures)
    close = bt.sort().closest(hg19_TSS.sort(), D="a")
    tss_count = 0
    total_peaks = len(bt)
    features_to_check = []
    for feature in close:
        if abs(int(feature[-1])) <= 2000 and int(feature[-1]) != -1:
            pos = ":".join(feature.fields[0:3])
            if pos not in features_to_check:
                tss_count += 1
                features_to_check.append(pos)

    idx_name = "-".join(efd[0:3])
    tss_count_df.set_value(idx_name, "TSS Count", tss_count)
    tss_count_df.set_value(idx_name, "Total Binding Sites", total_peaks)
    tss_count_df.set_value(idx_name, "Genomic Fraction", float(tss_count) / float(len(hg19_TSS)))
    pbar.update(i+1)

pbar.finish()


twriter = pd.ExcelWriter(os.path.join(data_dir, "temp.xlsx"))
tss_count_df.to_excel(twriter)
twriter.save()

# now let's compare those features against the HSF1 features

hsf1_encode_comparison_files = [os.path.join(hsf1_encode_comparison_dir, f) for f in os.listdir(hsf1_encode_comparison_dir) if "bed" in f]

# print hsf1_encode_comparison_files
hsf1_regions = pbt.BedTool(hsf1_regions_file)

idx = [feature.name for feature in hsf1_regions]
# print idx

cols = [os.path.basename(f).replace("hsf1_", "").replace("_le_500.bed", "") for f in hsf1_encode_comparison_files]
# print cols

# pprint.pprint(idx)

dist_df = pd.DataFrame(index=idx + ["Co-Occurance Enrichment"], columns=cols)

pbar = progressbar.ProgressBar(maxval=len(hsf1_encode_comparison_files), widgets=["Comparing HSF1 - ENCODE TF data", progressbar.Bar(), progressbar.Percentage()]).start()


for i, f in enumerate(hsf1_encode_comparison_files):
    col_name = os.path.basename(f).replace("hsf1_", "").replace("_le_500.bed", "")
    bt = pbt.BedTool(f)
    for feature in bt:
        dist_df.set_value(feature.name.replace(":", "-"), col_name, int(feature[-1]))

    #get the TSS count for that dataset
    tss_count = tss_count_df.get_value(col_name, "TSS Count")
    tss_fraction = tss_count_df.get_value(col_name, "Genomic Fraction")
    # print tss_count
    co_occ_freq = float(len(bt)) / float(tss_count)
    hsf1_frac = float(len(bt)) / float(len(dist_df.index))
    enrich = hsf1_frac / tss_fraction

    dist_df.set_value("Co-Occurance Enrichment", col=col_name, value=enrich)



    pbar.update(i+1)

pbar.finish()


hsf1_mod_regions = []
for feature in hsf1_regions:
    peak = str( int(math.floor((int(feature.fields[1]) + int(feature.fields[2])) /2.0) - int(feature.fields[1])) )
    hsf1_mod_regions.append(feature.fields + [peak])

hsf1_regions_with_peak = pbt.BedTool(hsf1_mod_regions)


pbar = progressbar.ProgressBar(maxval=len(encode_files_data), widgets=["Computing fisher exact tests", progressbar.Bar(), progressbar.Percentage()]).start()

for i, f in enumerate(encode_files_data):
    encode_file = f[3]
    col_name = "-".join(f[0:3])
    oddsr, pval = chrom_utils.compute_fisher_exact_TF_cobinding_from_bedfiles(hsf1_regions_with_peak,
                                                                              encode_file,
                                                                              hg19_TSS_file,
                                                                              tss_length=len(hg19_TSS),
                                                                              tf_distance_thresh=500)

    dist_df.set_value("Fisher Exact Test Odds Ratio", col=col_name, value=oddsr)
    dist_df.set_value("Fisher Exact P-value", col=col_name,value=pval)
    col_parts = col_name.split("-")
    cell = col_parts[0]
    tf = col_parts[1]
    exp = col_parts[2]
    dist_df.set_value("Cell Line", col=col_name, value=cell)
    dist_df.set_value("Transcription Factor", col=col_name, value=tf)
    dist_df.set_value("Experiment", col=col_name, value=exp)
    pbar.update(i+1)

pbar.finish()

print "Saving data"

writer = pd.ExcelWriter(os.path.join(data_dir, "hsf1_features_encode.xlsx"))
dist_df.to_excel(writer, sheet_name="HSF1 TF Distance")
tss_count_df.to_excel(writer, sheet_name="ENCODE TF TSS")
writer.save()
