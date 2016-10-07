import os
import pandas as pd
import pybedtools as pbt
import shutil
import progressbar
import chrom_utils
import math
import re

data_dir = os.path.join("..", "data_files")
cc_tf_dir = os.path.join(data_dir, "cc_tf_bedfiles")
hg19_TSS_file = os.path.join(data_dir, "hg19_TSS.ranged.1.bed")
hsf1_regions_file = os.path.join(data_dir, "hsf1_0.17_merged_regions_hg19.refiltered.bed")

hsf1_regions = pbt.BedTool(hsf1_regions_file)
hg19_TSS = pbt.BedTool(hg19_TSS_file)


cc_tfs_dirs = [os.path.join(cc_tf_dir, d) for d in os.listdir(cc_tf_dir) if os.path.isdir(os.path.join(cc_tf_dir,d)) and "RB" not in d]

tf_tss_df = pd.DataFrame(index=[os.path.basename(d) for d in cc_tfs_dirs], columns=["TSS Count", "Total Binding Sites", "Genomic Fraction"])
dist_df = pd.DataFrame(index=[feature.name for feature in hsf1_regions] + [ "Co-Occurance Enrichment"], columns=[os.path.basename(d) for d in cc_tfs_dirs])

#make a modified hsf1 bedfile for the Fisher exact test analysis that makes the last column the peak, assuming
#that the peak is in the middle of the region

hsf1_mod_regions = []
for feature in hsf1_regions:
    peak = str( int(math.floor((int(feature.fields[1]) + int(feature.fields[2])) /2.0) - int(feature.fields[1])) )
    hsf1_mod_regions.append(feature.fields + [peak])

hsf1_regions_with_peak = pbt.BedTool(hsf1_mod_regions)
print len(hsf1_regions_with_peak)

intermediate_regions_dir = os.path.join(data_dir, "hsf1_cctf_region_comparisons")
if not os.path.isdir(intermediate_regions_dir): os.mkdir(intermediate_regions_dir)
else:
    shutil.rmtree(intermediate_regions_dir)
    os.mkdir(intermediate_regions_dir)

intermediate_regions_tss_dir = os.path.join(data_dir, "hsf1_cctf_region_comparisons_tss")
if not os.path.isdir(intermediate_regions_tss_dir): os.mkdir(intermediate_regions_tss_dir)
else:
    shutil.rmtree(intermediate_regions_tss_dir)
    os.mkdir(intermediate_regions_tss_dir)

pbar = progressbar.ProgressBar(maxval=len(cc_tfs_dirs), widgets=["Processing CC TF files", progressbar.Bar(), progressbar.Percentage()]).start()

for i, d in enumerate(cc_tfs_dirs):
    cc_tf = os.path.basename(d)
    tf_bed_file = os.path.join(d, cc_tf + ".bed")
    # print os.path.isfile(tf_bed_file)
    cctf_bed = pbt.BedTool(tf_bed_file)
    #make a temp bed where the start and stop are the peak
    # peakfields = []
    # for feature in cctf_bed:
    #     peak = int(feature[-1])
    #     newstart = int(feature[1]) + peak
    #     newstop = int(feature[2]) + peak
    #     newfields = feature.fields
    #     newfields[1] = newstart
    #     newfields[2] = newstop
    #     peakfields.append(newfields)
    #
    # cctf_peak_bed = pbt.BedTool(peakfields)
    #see how many places the tf is near the TSS
    tss_close = cctf_bed.sort().closest(hg19_TSS.sort(), D="a")
    ex_fields = tss_close[0].fields
    chridxs = [j for j in range(len(ex_fields)) if re.match(re.compile("(?<=^)(chr\w+)(?=$)"), ex_fields[j])]
    tss_start_idx = chridxs[1]+1
    tf_peak_idx = chridxs[1]-1
    print len(tss_close)
    # print tss_close[0].fields
    # raw_input()
    tss_count = 0
    total = float(len(cctf_bed))
    features_to_check = []
    for feature in tss_close:
        if int(feature.fields[-1]) == -1: continue
        cctf_peak = int(tf_peak_idx)
        cctf_peak_pos = int(feature.fields[1]) + cctf_peak

        tss_start = int(feature.fields[tss_start_idx])
        tss_end = int(feature.fields[tss_start_idx+1])

        dist = min( abs(cctf_peak_pos - tss_start), abs(cctf_peak_pos - tss_end) )
        # print dist
        if dist <= 2000:
            pos = ":".join(feature.fields[tss_start_idx:tss_start_idx+3])
            if pos not in features_to_check:
                tss_count += 1
                features_to_check.append(pos)



    if tss_count >= total:
        print "something is wrong, more tss that features in bed file"
        print "saving temp bed"
        print "file: ", tf_bed_file
        print "temp file: ", os.path.join(data_dir, "temp.bed")
        tempbt = pbt.BedTool(features_to_check).sort()
        tempbt.saveas(os.path.join(data_dir, "temp.bed"))
        print "tss count: ", tss_count
        print "num from bedtool: ", len(tempbt)
        print "len of original bed file: ", len(cctf_bed)
        raw_input()


    genomic_frac = float(tss_count) / float(len(hg19_TSS))
    tf_tss_df.set_value(cc_tf, "TSS Count", float(tss_count))
    tf_tss_df.set_value(cc_tf, "Total Binding Sites", total)
    tf_tss_df.set_value(cc_tf, "Genomic Fraction", float(tss_count) / float(len(hg19_TSS)))

    #compare where it is binding relative to hsf1
    # original_features = [feature for feature in cctf_bed]
    hsf1_close = hsf1_regions.sort().closest(cctf_bed.sort(), D="a")
    hits = []
    uniques = []
    for feature in hsf1_close:
        cctf_peak = int(feature[-2])
        cctf_start = int(feature[5])
        cctf_peak_pos = cctf_peak + cctf_start

        hsf1_start = int(feature[1])
        hsf1_end = int(feature[2])
        hsf1_peak = int(math.floor( (hsf1_start + hsf1_end) / 2.0 ))

        dist = min(hsf1_peak - cctf_peak_pos, hsf1_peak - cctf_peak_pos)

        if abs(dist) <= 500:
            p = ":".join(feature.fields[0:3])
            if p not in uniques:
                fields = feature.fields
                fields[-1] = dist
                hits.append(fields)
                uniques.append(p)

    outname = os.path.join(intermediate_regions_dir, "hsf1_{cctf}_le_500.bed".format(cctf=cc_tf))
    hitsbt = pbt.BedTool(hits)
    hitsbt.sort().saveas(outname)

    for feature in hitsbt:
        dist_df.set_value(feature.name, cc_tf, value=float(feature[-1]))

    # print len(dist_df.index)
    # print len(hsf1_regions)
    # print len(hsf1_regions_with_peak)
    hg19_tss_len = len(hg19_TSS)

    co_occ_freq = float(len(hitsbt)) / float(tss_count)
    hsf1_frac = float(len(hitsbt)) / float(len(hsf1_regions))
    tf_frac = (float(tss_count) / float(hg19_tss_len))

    # print hsf1_frac
    # print len(hitsbt)
    # print len(hsf1_regions)
    # print tf_frac == genomic_frac

    enrich = hsf1_frac / tf_frac

    my_oddsr = ( float(len(hitsbt)) / (float(len(hsf1_regions)) - float(len(hitsbt))) ) / ( float(tss_count) / ( float(hg19_tss_len) - float(tss_count) ) )

    dist_df.set_value("Co-Occurance Enrichment", cc_tf, value=enrich)
    dist_df.set_value("My Odds Ratio Calc", cc_tf,value=my_oddsr )

    #compute the fisher exact test for these two transcription factors
    #use the built in chrom_utils function
    #should use the precomputed tss_length

    # print tf_bed_file
    oddsr, pv = chrom_utils.compute_fisher_exact_TF_cobinding_from_bedfiles(hsf1_regions_with_peak,
                                                                            tf_bed_file,
                                                                            hg19_TSS_file,
                                                                            tss_length=hg19_tss_len,
                                                                            tss_distance_thresh=2000,
                                                                            tf_distance_thresh=500)
    dist_df.set_value("Fisher Exact Odds Ratio", cc_tf, value=oddsr)
    dist_df.set_value("Fisher Exact P-value", cc_tf, value=pv)

    # print "HSF1 + TF: ", len(hitsbt)
    # print "HSF1 Total: ", len(hsf1_regions)
    # print "HSF1 peak total: ", len(hsf1_regions_with_peak)
    # print "TF TSS count: ", tss_count
    # print "Total TSS: ", hg19_tss_len

    pbar.update(i+1)

pbar.finish()

writer = pd.ExcelWriter(os.path.join(data_dir, "hsf1_features_cctf.xlsx"))
dist_df.to_excel(writer, sheet_name="HSF1 TF Distance")
tf_tss_df.to_excel(writer, sheet_name="CCTF TSS")
writer.save()
