import pybedtools as pbt
import os
import pprint


fly_peak_dir = os.path.join("..","data_files", "encode", "narrow_peaks", "fly")
fly_tss_file = os.path.join("..","data_files", "d_melanogaster_all_refined_TSS.ranged.bed")
fly_tss_bed = pbt.BedTool(fly_tss_file)

fly_peaks_files = [os.path.join(fly_peak_dir,f) for f in os.listdir(fly_peak_dir) if "S2-" in f]
pprint.pprint(fly_peaks_files)


output = []

for i, f1 in enumerate(fly_peaks_files):
    for j, f2 in enumerate(fly_peaks_files):
        if j <= i:
            continue
        else:
            i_bed = pbt.BedTool(f1)
            j_bed = pbt.BedTool(f2)
            i_j_int = i_bed.intersect(j_bed)
            i_j_tss = i_j_int.intersect(fly_tss_bed)
            output.append(i_j_tss)


print output
output[0].saveas(os.path.join("..", "data_files", "test_intersection.bed"))