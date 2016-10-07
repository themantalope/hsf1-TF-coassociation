import os
import pybedtools as pbt
import math

data_dir = os.path.join("..", "data_files")
mendillo_hs_dir = os.path.join(data_dir, "mendillo_HS")
hg18_TSS_file = os.path.join(data_dir, "hg18_TSS.ranged.0.bed")
mendillo_files = [os.path.join(mendillo_hs_dir, f) for f in os.listdir(mendillo_hs_dir) if "peak" not in f and "-R" not in f]


all_merged = pbt.BedTool(mendillo_files[0]).sort()
all_merged = all_merged.cat(*[pbt.BedTool(mendillo_files[1]), pbt.BedTool(mendillo_files[2])], postmerge=True, d=200)
# print len(all_merged)

hg18_TSS = pbt.BedTool(hg18_TSS_file).sort()

close = hg18_TSS.closest(all_merged.sort(), D="a")
good_features = []

unique_pos = []
for feature in close:
    if abs(int(feature[-1])) <= 1500 and int(feature[-1]) != -1:
        if ":".join(feature.fields[0:3]) not in unique_pos:
            good_features.append(feature.fields[6:9])
            unique_pos.append(":".join(feature.fields[0:3]))

short_count = 0
for i, feature in enumerate(good_features):
    start = int(feature[1])
    stop = int(feature[2])
    mid = int(math.floor((start + stop)/2.0))
    if abs(start - stop) < 400:
        newstart = str(mid - 200)
        newstop = str(mid + 200)
        feature[1] = newstart
        feature[2] = newstop
        short_count += 1
    feature.append("feature{n}".format(n=i)+"|"+"-".join(feature)+"-1")

print short_count

# print len(good_features)
outfile = os.path.join(data_dir, "hsf1_HS_merged_regions.bed")
outbed = pbt.BedTool(good_features).sort()
outbed.saveas(outfile)
print len(outbed)