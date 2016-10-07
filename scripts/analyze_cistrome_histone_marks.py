import os
import pybedtools as pbt
import pandas as pd
import chrom_utils
import mygene

data_dir = os.path.join("..", "data_files")
cistrome_dir = os.path.join(data_dir, "cistrome_files")
hg19_TSS_bed_file = os.path.join(data_dir, "hg19_TSS.ranged.0.unique.bed")
hsf1_regions_bed_file = os.path.join(data_dir, "hsf1_0.17_merged_regions_hg19.bed")
hsf1_basal_regions_file = os.path.join(data_dir, "hsf1_nonHS_exclusive.bed")
hsf1_HS_regions_file = os.path.join(data_dir, "hsf1_HS_exclusive.bed")
cell_cycle_phase_file = os.path.join(data_dir, "cell_cycle_phase_ids.xlsx")

cell_cycle_phase = pd.read_excel(open(cell_cycle_phase_file))
mg = mygene.MyGeneInfo()

cell_cycle_bed_fields = []
for row in cell_cycle_phase.iterrows():
    idx = row[0]
    rowdata = row[1]
    newfields = rowdata["Genomic Position"].split(":")
    name = "."
    if not pd.isnull(rowdata["DW Phase"]): name = rowdata["DW Phase"]
    if not pd.isnull(rowdata["GW Phase"]): name = rowdata["GW Phase"]
    newfields.append(idx+":"+name)
    # print newfields
    # raw_input()
    cell_cycle_bed_fields.append(newfields)

cc_bed = pbt.BedTool(cell_cycle_bed_fields)


hg19_TSS = pbt.BedTool(hg19_TSS_bed_file)
hsf1_regions = pbt.BedTool(hsf1_regions_bed_file)
hsf1_basal_regions = pbt.BedTool(hsf1_basal_regions_file)
hsf1_HS_regions = pbt.BedTool(hsf1_HS_regions_file)


histone_folders = [os.path.join(cistrome_dir,d) for d in os.listdir(cistrome_dir) if os.path.isdir(cistrome_dir) and "H" in d]

outdf = pd.DataFrame(index=[os.path.basename(d) for d in histone_folders], columns=["HSF1", "HSF1 Frac", "HSF1 - basal", "HSF1 - basal Frac", "HSF1 - HS", "HSF1 - HS Frac", "TSS", "TSS Frac", "HSF1/TSS Frac Ratio", "HSF1 - basal/TSS Frac Ratio", "HSF1 - HS/TSS Frac Ratio"])

outdfs = {"comparisons":outdf}


#make an "extended tss" temp bedfile

extended = []
for feature in hg19_TSS:
    fields = feature.fields
    newstart = int(feature[1]) - 1000
    newstop = int(feature[1]) + 1000
    if newstart < 0: newstart = 1
    if newstart > newstop :
        print "start: ", newstart, "stop: ", newstop
        raw_input()
    fields[1] = str(newstart)
    fields[2] = str(newstop)
    extended.append(fields)

extended_TSS = pbt.BedTool(extended)

for mark_dir in histone_folders:
    mark_file = os.path.join(mark_dir,os.path.basename(mark_dir) + ".bed")
    mark = os.path.basename(mark_dir)
    mark_bed = pbt.BedTool(mark_file)

    #count the number of features where there is at least one peak overlapping

    # mark_tss_close = extended_TSS.sort().closest(mark_bed.sort(), D="a")
    # tss_peak_count = 0
    # for feature in mark_tss_close:
    #     if abs(int(feature[-1])) <= 2000 and int(feature[-1]) != -1:
    #         tss_peak_count += 1
    #
    # tss_peak_count = float(tss_peak_count)

    tss_intersection = extended_TSS.sort().intersect(mark_bed)
    # tss_unique_intersection = []
    # uniques = []
    # for feature in tss_intersection:
    #     if ":".join(feature.fields[0:3]) not in uniques:
    #         tss_unique_intersection.append(feature.fields)
    #         uniques.append(":".join(feature.fields[0:3]))
    #
    #
    # tss_intersection = pbt.BedTool(tss_unique_intersection)
    tss_peak_count = float(len(tss_intersection))

    outdf.set_value(mark, "TSS", tss_peak_count)
    tss_frac = tss_peak_count / chrom_utils.calculate_total_basepair_length(extended_TSS)
    outdf.set_value(mark, "TSS Frac", tss_frac)

    #get the intersection with HSF1
    hsf1_intersection = hsf1_regions.intersect(mark_bed)
    hsf1_intersection_count = float(len(hsf1_intersection))
    outdf.set_value(mark, "HSF1", hsf1_intersection_count)
    hsf1_frac = hsf1_intersection_count / chrom_utils.calculate_total_basepair_length(hsf1_regions)
    outdf.set_value(mark, "HSF1 Frac", hsf1_frac)
    outdf.set_value(mark, "HSF1/TSS Frac Ratio", hsf1_frac / tss_frac)

    #get intersection with HSF1 basal
    hsf1_basal_intersection = hsf1_basal_regions.intersect(mark_bed)
    hsf1_basal_count = float(len(hsf1_basal_intersection))
    hsf1_basal_frac = hsf1_basal_count / chrom_utils.calculate_total_basepair_length(hsf1_basal_regions)
    outdf.set_value(mark, "HSF1 - basal", hsf1_basal_count)
    outdf.set_value(mark, "HSF1 - basal Frac", hsf1_basal_frac)
    outdf.set_value(mark, "HSF1 - basal/TSS Frac Ratio", hsf1_basal_frac / tss_frac)

    #get intersection with HSF1 HS
    hsf1_HS_intersection = hsf1_HS_regions.intersect(mark_bed)
    hsf1_HS_count = float(len(hsf1_HS_intersection))
    hsf1_HS_frac = hsf1_HS_count / chrom_utils.calculate_total_basepair_length(hsf1_HS_regions)
    outdf.set_value(mark, "HSF1 - HS", hsf1_HS_count)
    outdf.set_value(mark, "HSF1 - HS Frac", hsf1_HS_frac)
    outdf.set_value(mark, "HSF1 - HS/TSS Frac Ratio", hsf1_HS_frac / tss_frac)

    #now let's find which genes are close to the HSF1 intersection data


    hsf1_intersection_close = hsf1_intersection.closest(cc_bed.sort(), D="a")
    uniq_locs = []
    good_locs = []
    for feature in hsf1_intersection_close:
        if ":".join(feature.fields[4:6]) in uniq_locs: continue
        phase = feature.fields[-2]
        phase = phase.split(":")[1]
        if phase == ".": continue
        else:
            good_locs.append(feature.fields)

    mark_genes_df = pd.DataFrame(columns=["Symbol", "Phase"])
    genes_to_search = []
    phase_pairs = []
    for gl in good_locs:
        phase = gl[-2].split(":")[1]
        gene = gl[-2].split(":")[0]
        # mark_genes_df.set_value(gene, "Phase", phase)
        genes_to_search.append(gene)
        phase_pairs.append((gene, phase))

    query = mg.querymany(genes_to_search, scopes=["refseq"], fields="symbol", returnall=True)
    for q in query["out"]:
        if "symbol" not in q: continue
        else:
            if q["symbol"] in mark_genes_df["Symbol"].values.tolist(): continue
            else:
                mark_genes_df.set_value(index=q["query"], col="Symbol", value=q["symbol"])
                cphase = [x[1] for x in phase_pairs if x[0] == q["query"]][0]
                mark_genes_df.set_value(index=q["query"], col="Phase", value=cphase)

    outdfs[mark] = mark_genes_df



writer = pd.ExcelWriter(os.path.join(data_dir, "hsf1_histone_marks.xlsx"))
for k in outdfs.keys():
    outdfs[k].to_excel(writer, sheet_name=k)

writer.save()






