import pybedtools as pbt
import os
import pandas as pd

data_file = os.path.join("..", "data_files", "hg18_TSS.txt")

data = pd.read_table(data_file)

features = []

g = pbt.genome_registry

gap = 0

for row in data.iterrows():
    idx = row[0]
    rowdata = row[1]
    newstart = rowdata["txStart"] - gap
    newstop = rowdata["txStart"] + gap

    if newstart < g.hg18[rowdata["chrom"]][0]:
        newstart = g.hg18[rowdata["chrom"]][0]

    if newstop > g.hg18[rowdata["chrom"]][1]:
        newstop = g.hg18[rowdata["chrom"]][1]


    features.append((rowdata["chrom"], newstart, newstop, rowdata["#name"], ".", rowdata["strand"]))


a = pbt.BedTool(features)
a.saveas(os.path.join("..", "data_files", "hg18_TSS.ranged.{g}.bed".format(g=gap)))
print "Done"