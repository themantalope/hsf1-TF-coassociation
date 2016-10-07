import os
import pybedtools as pbt
import datetime

data_dir = os.path.join("..", "data_files")
hg19_tss_file = os.path.join(data_dir, "hg19_TSS.ranged.1.bed")
hsf1_regions_file = os.path.join(data_dir, "hsf1_0.17_hg19.bed")

hg19_tss = pbt.BedTool(hg19_tss_file)
hsf1_regions = pbt.BedTool(hsf1_regions_file)

#also need to remerge the file after lifting

hsf1_regions_merged = hsf1_regions.sort().merge()



close = hsf1_regions_merged.sort().closest(hg19_tss.sort(), D="a")
# print close[0].fields; raw_input()
good_hsf1_regions = []
uniques = []
print len(close)
for feature in close:
    if int(feature.fields[-1]) == -1:continue
    if abs(int(feature.fields[-1])) <= 2000:
        if ":".join(feature.fields[4:7]) not in uniques:
            good_hsf1_regions.append(feature.fields[0:3])
            uniques.append(":".join(feature.fields[4:7]))






# closehsf1 = hsf1_regions.sort().closest(hg19_tss.sort(), D="a")
# test_regions = []
# for feature in closehsf1:
#     if int(feature.fields[-1]) == -1: continue
#     if abs(int(feature.fields[-1])) <= 2000 : test_regions.append(feature.fields)
#
#
# if len(test_regions) != len(good_hsf1_regions):
#     print test_regions
#     temp = pbt.BedTool(test_regions).sort().saveas(os.path.join(data_dir, "temp{dt}.bed".format(dt=datetime.datetime.now().isoformat())))
#     raw_input()
#     print good_hsf1_regions


# print "len of test: ", len(test_regions)
# print "len of good: ", len(good_hsf1_regions)


# str1 = "\n".join(["\t".join(feature.fields) for feature in hsf1_regions])
# str2 = "\n".join(["\t".join(gf) for gf in good_hsf1_regions])

# print "start and result are exactly the same: ", str1 == str2




print "number of good features after filter: ", len(good_hsf1_regions)
print "number of features total: ", len(hsf1_regions)
print "difference: ", len(hsf1_regions) - len(good_hsf1_regions)




for i, feature in enumerate(good_hsf1_regions):
    feature.append("feature{n}|{c}-{s}-{e}-1".format(n=i, c=feature[0],s=feature[1], e=feature[2]))

# print good_hsf1_regions[0]; raw_input()
# print good_hsf1_regions

outfile = os.path.join(data_dir, "hsf1_0.17_merged_regions_hg19.refiltered.bed")
pbt.BedTool(good_hsf1_regions).sort().saveas(outfile)