import os
import chrom_utils
import re

ext_patt = re.compile("(.bed)(?=$)")

data_dir = os.path.join("..", "data_files")
cc_tf_dir = os.path.join(data_dir, "cc_tf_bedfiles" )

cc_tf_files = []
cc_tf_dirs = [os.path.join(cc_tf_dir,d) for d in os.listdir(cc_tf_dir) if os.path.isdir(os.path.join(cc_tf_dir,d))]

for d in cc_tf_dirs:
    ccfs = [os.path.join(d,f) for f in os.listdir(d) if os.path.isfile(os.path.join(d,f)) and len(re.findall(ext_patt, f)) != 0 and "peak" not in f]
    cc_tf_files += ccfs


for f in cc_tf_files:
    print f
    chrom_utils.set_region_start_stop_based_on_peak(f)