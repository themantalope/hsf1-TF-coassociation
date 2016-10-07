import pybedtools as pbt
import os
import pandas as pd
import chrom_utils
import progressbar
import pickle

data_dir = os.path.join("..", "data_files")
tss_bed_file = os.path.join(data_dir, "hg19_TSS.ranged.2500.bed")
mendillo_dir = os.path.join(data_dir, "mendillo_merged")
encode_dir = os.path.join(data_dir, "encode", "narrow_peaks", "human")
encode_cell_types_file = os.path.join(data_dir, "human_encode_epithelial_cell_types.txt")

encode_files_file = os.path.join(data_dir, "encode_peak_bed_files.pickle")
mendillo_files_file = os.path.join(data_dir, "mendillo_peak_bed_files.pickle")

encode_files = pickle.load(open(encode_files_file))
mendillo_files = pickle.load(open(mendillo_files_file))



mendillo_file_data = [(f, chrom_utils.get_cell_type_from_mendillo_filename(os.path.basename(f)), chrom_utils.get_tf_type_from_mendillo_filename(os.path.basename(f))) for f in mendillo_files]

encode_file_data = [(f, chrom_utils.get_cell_type_from_encode_filename(os.path.basename(f)), chrom_utils.get_tf_type_from_encode_filename(os.path.basename(f))) for f in encode_files]


#
# print mendillo_file_data
# raw_input()
#
#
# print encode_file_data
# raw_input()



all_file_data = mendillo_file_data + encode_file_data


tss_bed = pbt.BedTool(tss_bed_file)

output = []

for i , f1 in enumerate(all_file_data):
    f1_bed = pbt.BedTool(f1[0])
    print "Comparing against: ", f1[1] + "-" + f1[2]
    pbar = progressbar.ProgressBar(maxval=len(all_file_data), widgets=["Getting co-associations", progressbar.Bar()]).start()
    for j, f2 in enumerate(all_file_data):
        if j <= i: continue
        else:
            f2_bed = pbt.BedTool(f2[0])

            output.append((f1[1] + "-" + f1[2],
                           f2[1] + "-" + f2[2],
                           chrom_utils.calculate_TSS_co_binding_minimum_overlap_ratio(tss_bed=tss_bed,a_peak_bed=f1_bed,b_peak_bed=f2_bed)))
            pbar.update(j+1)
            # print j, f2[0]

    pbar.finish()

rows = []
cols = []

for x in output:
    if x[0] not in rows:
        rows.append(x[0])

    if x[1] not in cols:
        cols.append(x[1])

df = pd.DataFrame(index=rows, columns=cols)

for x in output:
    df.set_value(x[0], x[1], x[2])

writer = pd.ExcelWriter(os.path.join(data_dir, "mendillo_encode_comparison.xlsx"))
df.to_excel(writer)
writer.save()

