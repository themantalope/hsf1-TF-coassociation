import pickle
import os
import chrom_utils
import pybedtools as pbt

data_dir = os.path.join("..", "data_files")
epi_cell_types_file = os.path.join(data_dir, "human_encode_epithelial_cell_types.txt")
mendillo_dir = os.path.join(data_dir, "mendillo")
encode_data_dir = os.path.join(data_dir, "encode", "narrow_peaks", "human")
tss_bed_file = os.path.join(data_dir, "hg19_TSS.ranged.2500.bed")


#get a list of the acceptable cell types
f = open(epi_cell_types_file, "r")
epi_cell_types = f.readlines()
f.close()

epi_cell_types = [e.strip() for e in epi_cell_types]

# print os.listdir(encode_data_dir)
# print epi_cell_types

encode_files = [os.path.join(encode_data_dir,f) for f in os.listdir(encode_data_dir) if os.path.isfile(os.path.join(encode_data_dir,f))  and "Tfbs" in f]

encode_files = [f for f in encode_files if "Histone" not in f]

encode_files = [f for f in encode_files if chrom_utils.get_cell_type_from_encode_filename(f) in epi_cell_types]

encode_file_data = [(f, chrom_utils.get_cell_type_from_encode_filename(f), chrom_utils.get_tf_type_from_encode_filename(f)) for f in encode_files]

mendillo_files = [os.path.join(mendillo_dir,f) for f in os.listdir(mendillo_dir) if os.path.isfile(os.path.join(mendillo_dir,f)) and ".bed" in f]

mendillo_file_data = [(f, chrom_utils.get_cell_type_from_mendillo_filename(os.path.basename(f)), chrom_utils.get_tf_type_from_mendillo_filename(os.path.basename(f))) for f in mendillo_files]


all_file_data = mendillo_file_data + encode_file_data

#now we can start computing the jaccard index between file types.

tss_bed = pbt.BedTool(tss_bed_file)

output = []
errors = []


for i, f1 in enumerate(all_file_data):
    f1_bt = pbt.BedTool(f1[0])
    for j, f2 in enumerate(all_file_data):
        if j <= i:
            continue
        elif "HSF" not in f1[2] or "HSF" not in f2[2]:
            continue
        else:
            print f1[0]
            print f2[0]
            try:
                # print chrom_utils.compute_tss_jaccard_intersection(tss_bed=tss_bed, a_bed=pbt.BedTool(f1[0]),b_bed=pbt.BedTool(f2[0]))
                # raw_input()
                output.append((f1[1]+"-"+f1[2],
                               f2[1] + "-" + f2[2],
                               chrom_utils.calculate_TSS_co_binding_minimum_overlap_ratio(tss_bed=tss_bed, a_bed=f1_bt, b_bed=f2[0])))
            except BaseException, e:
                print "hit error :/"
                errors.append((f1[0], f2[0], str(e)))
                print str(e)
                raw_input()


#dump it to a pickle for now

pickle.dump(output, open(os.path.join(data_dir, "initial_output.pickle"), "w"))
pickle.dump(errors, open(os.path.join(data_dir, "error_files_from_initial_test.pickle"), "w"))