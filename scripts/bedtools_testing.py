import os
import pybedtools as pbt
import chrom_utils


wd = '/Users/antalek/Documents/School/Research/aging_regulation/data_files/mendillo_merged'
os.chdir(wd)

a_bed = pbt.BedTool("NCI1703_Regions_HSF1.bed")
b_bed = pbt.BedTool("NCIH441_Regions_HSF1.bed")
tss_bed = pbt.BedTool("../hg19_TSS.ranged.2500.bed")

chrom_utils.set_region_start_stop_based_on_peak("NCI1703_Regions_HSF1.bed")

print chrom_utils.calculate_TSS_co_binding_minimum_overlap_ratio(tss_bed=tss_bed, a_bed=a_bed,b_bed=b_bed)

g = pbt.genome_registry.hg19

print len(a_bed.slop(b=5000,g=g ).remove_invalid().intersect(b_bed))