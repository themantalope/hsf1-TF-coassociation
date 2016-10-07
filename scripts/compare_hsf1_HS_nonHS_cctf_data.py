import pandas as pd
import os
import pybedtools as pbt
import chrom_utils
import pickle
import math
import progressbar

data_dir = os.path.join("..", "data_files")
hsf1_regions_file = os.path.join(data_dir, "hsf1_0.17_merged_regions_hg19.refiltered.bed")
hsf1_HS_regions_file = os.path.join(data_dir, "hsf1_HS_exclusive.bed")
hsf1_nonHS_regions_file = os.path.join(data_dir, "hsf1_nonHS_exclusive.bed")
hsf1_HS_intersection_file = os.path.join(data_dir, "hsf1_HS_nonHS_intersection_hg19.bed")
hg19_tss_file = os.path.join(data_dir, "hg19_TSS.ranged.1.bed")
cc_tf_dir = os.path.join(data_dir, "cc_tf_bedfiles")


cc_tfs_dirs = [os.path.join(cc_tf_dir, d) for d in os.listdir(cc_tf_dir) if os.path.isdir(os.path.join(cc_tf_dir,d)) and "RB" not in d]


#make temporary peaks for the hsf1 files

hsf1_regions = pbt.BedTool(hsf1_regions_file)
hsf1_peak_regions = []
for feature in hsf1_regions:
    peak = str(int(math.floor( (int(feature.fields[1]) + int(feature.fields[2]))/2.0 ) - int(feature.fields[1])))
    hsf1_peak_regions.append(feature.fields + [peak])

hsf1_regions = pbt.BedTool(hsf1_peak_regions).sort()

#hsf1 HS
hsf1_HS_regions = pbt.BedTool(hsf1_HS_regions_file)
hsf1_HS_peak_regions = []
for feature in hsf1_HS_regions:
    peak = str(int(math.floor( (int(feature.fields[1]) + int(feature.fields[2]))/2.0 ) - int(feature.fields[1])))
    hsf1_HS_peak_regions.append(feature.fields + [peak])

hsf1_HS_regions = pbt.BedTool(hsf1_HS_peak_regions).sort()

#hsf1 nonHS
hsf1_nonHS_regions = pbt.BedTool(hsf1_nonHS_regions_file)
hsf1_nonHS_peak_regions = []
for feature in hsf1_nonHS_regions:
    peak = str(int(math.floor( (int(feature.fields[1]) + int(feature.fields[2]))/2.0 ) - int(feature.fields[1])))
    hsf1_nonHS_peak_regions.append(feature.fields + [peak])

hsf1_nonHS_regions = pbt.BedTool(hsf1_nonHS_peak_regions).sort()

#hsf1 HS intersection
hsf1_HS_intersection_regions = pbt.BedTool(hsf1_HS_intersection_file)
hsf1_HS_intersection_peak_regions = []
for feature in hsf1_HS_intersection_regions:
    peak = str(int(math.floor( (int(feature.fields[1]) + int(feature.fields[2]))/2.0 ) - int(feature.fields[1])))
    hsf1_HS_intersection_peak_regions.append(feature.fields + [peak])

hsf1_HS_intersection_regions = pbt.BedTool(hsf1_HS_intersection_peak_regions).sort()

#do the analysis loop

outdf = pd.DataFrame(columns=["HSF1 Regions OR","HSF1 Regions Pv", "HSF1 HS Regions OR","HSF1 HS Regions Pv", "HSF1 nonHS Regions OR","HSF1 nonHS Regions Pv", "HSF1 HS Intersection OR","HSF1 HS Intersection Pv"])

hg19_tss = pbt.BedTool(hg19_tss_file).sort()

pbar = progressbar.ProgressBar(maxval=len(cc_tfs_dirs), widgets=["Calculating Fisher Exact Tests", progressbar.Bar(), progressbar.Percentage()]).start()


for i, d in enumerate(cc_tfs_dirs):
    cc_tf = os.path.basename(d)
    tf_bed_file = os.path.join(d, cc_tf + ".bed")

    row_name = cc_tf

    hsf1_oddsr, hsf1_pv = chrom_utils.compute_fisher_exact_TF_cobinding_from_bedfiles(hsf1_regions,
                                                                                      tf_bed_file,
                                                                                      tss_bedfile=hg19_tss_file,
                                                                                      tss_length=len(hg19_tss),
                                                                                      tf_distance_thresh=500,
                                                                                      tss_distance_thresh=2000)

    hsf1_HS_oddsr, hsf1_HS_pv = chrom_utils.compute_fisher_exact_TF_cobinding_from_bedfiles(hsf1_HS_regions,
                                                                                            tf_bed_file,
                                                                                            tss_bedfile=hg19_tss,
                                                                                            tss_length=len(hg19_tss),
                                                                                            tf_distance_thresh=500)

    hsf1_nonHS_oddsr, hsf1_nonHS_pv = chrom_utils.compute_fisher_exact_TF_cobinding_from_bedfiles(hsf1_nonHS_regions,
                                                                                                  tf_bed_file,
                                                                                                  tss_bedfile=hg19_tss,
                                                                                                  tss_length=len(
                                                                                                      hg19_tss),
                                                                                                  tf_distance_thresh=500)

    hsf1_HS_inters_oddsr, hsf1_HS_inters_pv = chrom_utils.compute_fisher_exact_TF_cobinding_from_bedfiles(
        hsf1_HS_intersection_regions,
        tf_bed_file,
        tss_bedfile=hg19_tss,
        tss_length=len(hg19_tss),
        tf_distance_thresh=500)

    outdf.set_value(row_name, "HSF1 Regions OR", value=hsf1_oddsr)
    outdf.set_value(row_name, "HSF1 Regions Pv", value=hsf1_pv)

    outdf.set_value(row_name, "HSF1 HS Regions OR", value=hsf1_HS_oddsr)
    outdf.set_value(row_name, "HSF1 HS Regions Pv", value=hsf1_HS_pv)

    outdf.set_value(row_name, "HSF1 nonHS Regions OR", value=hsf1_nonHS_oddsr)
    outdf.set_value(row_name, "HSF1 nonHS Regions Pv", value=hsf1_nonHS_pv)

    outdf.set_value(row_name, "HSF1 HS Intersection OR", value=hsf1_HS_inters_oddsr)
    outdf.set_value(row_name, "HSF1 HS Intersection Pv", value=hsf1_HS_inters_pv)



    pbar.update(i + 1)
pbar.finish()

writer = pd.ExcelWriter(os.path.join(data_dir, "hsf1_HS_nonHS_cctf_feature_comparisons.xlsx"))
outdf.to_excel(writer)
writer.save()