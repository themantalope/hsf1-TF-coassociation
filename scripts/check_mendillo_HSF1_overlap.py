import os
import pybedtools as pbt
import pprint
import chrom_utils

mendillo_dir = os.path.join("..", "data_files", "mendillo_sorted")

mendillo_files = [os.path.join(mendillo_dir,f) for f in os.listdir(mendillo_dir) if ".bed" in f]

tss_bed_file = os.path.join("..","data_files", "hg19_TSS.ranged.2500.bed")

for i, f1 in enumerate(mendillo_files):
    f1_bt = pbt.BedTool(f1)
    print "{f} overlap with human tss: ".format(f=f1), len(pbt.BedTool(tss_bed_file).intersect(f1_bt))
    print "len of {f}: ".format(f=f1), len(f1_bt)
    print "fraction of hsf1 tss overlap: ", float(len(pbt.BedTool(tss_bed_file).intersect(f1_bt))) / float(len(f1_bt))
    # for j, f2 in enumerate(mendillo_files):
    #     if j <= i: continue
    #     else:
    #         f2_bt = pbt.BedTool(f2)
    #         print "{c1} and {c2}".format(c1=chrom_utils.get_cell_type_from_mendillo_filename(f1), c2=chrom_utils.get_cell_type_from_mendillo_filename(f2))
    #         print "number of overlapping sites: ", len(f1_bt.intersect(f2_bt))
    #         print "number of sites non overlapping: ", len(f1_bt) + len(f2_bt) - len(f1_bt.intersect(f2_bt))
    #         print "overlap ratio: ", float(len(f1_bt.intersect(f2_bt))) / float(len(f1_bt) + len(f2_bt) - len(f1_bt.intersect(f2_bt)))
