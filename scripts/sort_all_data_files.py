import pybedtools as pbt
import os
import progressbar

data_dir = os.path.join("..", "data_files")
mendillo_dir = os.path.join(data_dir, "mendillo")
encode_dir = os.path.join(data_dir, "encode", "narrow_peaks", "human")
endcode_sorted_dir = os.path.join(data_dir, "encode", "narrow_peaks","human_sorted")
mendillo_sorted_dir = os.path.join(data_dir, "mendillo_sorted")

if not os.path.isdir(endcode_sorted_dir):
    os.mkdir(endcode_sorted_dir)

if not os.path.isdir(mendillo_sorted_dir):
    os.mkdir(mendillo_sorted_dir)


encode_files = [os.path.join(encode_dir,f) for f in os.listdir(encode_dir) if "region" in f]
mendillo_files = [os.path.join(mendillo_dir,f) for f in os.listdir(mendillo_dir) if "bed" in f]
# print encode_files
# print mendillo_files



pbar = progressbar.ProgressBar(maxval=len(encode_files), widgets=["Sorting encode files", progressbar.Bar()]).start()

for idx, f in enumerate(encode_files):
    t = pbt.BedTool(f)
    end_fn = os.path.basename(t.fn)
    t.sort(chrThenSizeA=True).saveas(os.path.join(endcode_sorted_dir, end_fn))
    pbar.update(idx+1)

pbar.finish()

pbar = progressbar.ProgressBar(maxval=len(mendillo_files), widgets=["Sorting mendillo files", progressbar.Bar()]).start()


for idx, f in enumerate(mendillo_files):
    t = pbt.BedTool(f)
    end_fn = os.path.basename(t.fn)
    t.sort(chrThenSizeA=True).saveas(os.path.join(mendillo_sorted_dir, end_fn))
    pbar.update(idx+1)

pbar.finish()



