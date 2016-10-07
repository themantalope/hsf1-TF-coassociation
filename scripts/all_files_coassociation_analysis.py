import os
import chrom_utils
import timeit
import cPickle as pickle
import pandas as pd

data_dir = os.path.join("..", "data_files")
tss_bed_file = os.path.join(data_dir, "hg19_TSS.ranged.2500.bed")
mendillo_list_file = os.path.join(data_dir, "mendillo_peak_bed_files.pickle")
encode_list_file = os.path.join(data_dir, "encode_peak_bed_files.pickle")

encode_files = pickle.load(open(encode_list_file))
mendillo_files = pickle.load(open(mendillo_list_file))



all_file_data = mendillo_files + encode_files
all_files_tuple = []
for i, t1 in enumerate(all_file_data):
    for j, t2 in enumerate(all_file_data):
        if j <= i:continue
        else:

            all_files_tuple.append((t1, t2))

print "Running multiprocessing"
start = timeit.default_timer()
min_overlap_ratio = chrom_utils.calculate_TSS_co_binding_minimum_overlap_ratio_multiprocessing(tss_bed_file, all_files_tuple)
end = timeit.default_timer()
print "Multiprocessing run time: ", (end - start)

final_output = []

assert len(min_overlap_ratio) == len(all_files_tuple)

for i in range(len(all_files_tuple)):
    t = all_files_tuple[i]
    final_output.append((chrom_utils.get_cell_type_from_filename(os.path.basename(t[0])) + "-" + chrom_utils.get_tf_type_from_filename(os.path.basename(t[0])) + "-" + chrom_utils.get_exp_type_from_filename(os.path.basename(t[0])),
                         chrom_utils.get_cell_type_from_filename(os.path.basename(t[1])) + "-" + chrom_utils.get_tf_type_from_filename(os.path.basename(t[1])) + "-" + chrom_utils.get_exp_type_from_filename(os.path.basename(t[1])) ,
                        min_overlap_ratio[i]))


mcf10_df = chrom_utils.convert_tuple_list_to_dataframe(final_output)


writer = pd.ExcelWriter(os.path.join(data_dir,"all_files_TSS_coassociation.xlsx"))
mcf10_df.to_excel(writer)
writer.save()