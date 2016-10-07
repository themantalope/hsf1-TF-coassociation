import os
import chrom_utils
import pickle


data_dir = os.path.join("..", "data_files")
mendillo_dir = os.path.join(data_dir, "mendillo_merged")
encode_dir = os.path.join(data_dir, "encode", "narrow_peaks", "human")

encode_epi_cell_file = os.path.join(data_dir, "human_encode_epithelial_cell_types.txt")
f = open(encode_epi_cell_file)
encode_cell_types = f.readlines()
f.close()

encode_cell_types= [e.strip() for e in encode_cell_types]

encode_files = [os.path.join(encode_dir, f) for f in os.listdir(encode_dir) if "Tfbs" in f and ".bed" not in f]
encode_files = [f for f in encode_files if chrom_utils.get_cell_type_from_encode_filename(os.path.basename(f)) in encode_cell_types]

pickle.dump(encode_files, open(os.path.join(data_dir,"encode_files.pickle"),"w"))

mendillo_files = [os.path.join(mendillo_dir, f) for f in os.listdir(mendillo_dir) if "bed" in f]


pickle.dump(mendillo_files, open(os.path.join(data_dir, "mendillo_files.pickle"), "w"))

all_files = mendillo_files + encode_files

for f in all_files:
    chrom_utils.set_region_start_stop_based_on_peak(f)


encode_peak_bed_files = [os.path.join(encode_dir,f) for f in os.listdir(encode_dir) if "peak.bed" in f]
mendillo_peak_bed_files = [os.path.join(mendillo_dir, f) for f in os.listdir(mendillo_dir) if "peak.bed" in f]

pickle.dump(encode_peak_bed_files, open(os.path.join(data_dir, "encode_peak_bed_files.pickle"), "w"))
pickle.dump(mendillo_peak_bed_files, open(os.path.join(data_dir, "mendillo_peak_bed_files.pickle"), "w"))


print "Done"
