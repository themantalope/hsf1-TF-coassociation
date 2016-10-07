import pybedtools as pbt
import os
import pandas as pd

data_file = os.path.join("..", "data_files", "hg19_TSS.txt")

data = pd.read_table(data_file)


start_groups = data.groupby("txStart")
# print len(start_groups)
# print start_groups.groups.values()[0]
# raw_input()

features = []

g = pbt.genome_registry

gap = 1


all_startgroups = start_groups.groups

for startgroup in all_startgroups:
    idx = all_startgroups[startgroup][0]
    rowdata = data.loc[idx,:]

    if "NR" in rowdata["#name"]: continue
    if "_" in rowdata["chrom"]: continue

    newstart = rowdata["txStart"] - gap
    newstop = rowdata["txStart"] + gap

    if newstart < g.hg19[rowdata["chrom"]][0]:
        newstart = g.hg19[rowdata["chrom"]][0]

    if newstop > g.hg19[rowdata["chrom"]][1]:
        newstop = g.hg19[rowdata["chrom"]][1]


    features.append((rowdata["chrom"], newstart, newstop, rowdata["#name"], ".", rowdata["strand"]))

# for row in data.iterrows():
#     idx = row[0]
#     rowdata = row[1]
#
#     if "NR" in rowdata["#name"]: continue
#     if "_" in rowdata["chrom"]: continue
#
#     newstart = rowdata["txStart"] - gap
#     newstop = rowdata["txStart"] + gap
#
#     if newstart < g.hg19[rowdata["chrom"]][0]:
#         newstart = g.hg19[rowdata["chrom"]][0]
#
#     if newstop > g.hg19[rowdata["chrom"]][1]:
#         newstop = g.hg19[rowdata["chrom"]][1]
#
#
#     features.append((rowdata["chrom"], newstart, newstop, rowdata["#name"], ".", rowdata["strand"]))


#next we need to check if the features are redundant or within 2000 bp of each other

good_features = []

i = 0
while i < len(features) - 1:
    curfeatures = features[i]
    nextfeatures = features[i+1]

    dist = abs( int(curfeatures[2]) - int(nextfeatures[1]) )


    if dist <= 2000:
        good_features.append(curfeatures)
        while dist <= 2000:
            if i == len(features) - 1: break
            dist = abs( int(curfeatures[2]) - int(features[i][1]))
            i+=1
    else:
        good_features.append(curfeatures)

    i+=1



a = pbt.BedTool(good_features).sort()
a.saveas(os.path.join("..", "data_files", "hg19_TSS.ranged.{g}.bed".format(g=gap)))
print "Done"