import pandas as pd
import os
import mygene
import pybedtools as pbt


data_dir = os.path.join("..", "data_files")
tss_bed_file = os.path.join(data_dir, "hg19_TSS.ranged.2500.bed")
mcf10_bed_file = os.path.join(data_dir, "mendillo_merged", "MCF10_Regions_HSF1.bed")

tss_bed = pbt.BedTool(tss_bed_file)
mcf10_bed = pbt.BedTool(mcf10_bed_file)

tss_mcf10_i = tss_bed.intersect(mcf10_bed).sort()

print "number of intersections: ", len(tss_mcf10_i)
print "number of features in mcf10: ", len(mcf10_bed)

genes = {f.name for f in tss_mcf10_i}
print genes

mg = mygene.MyGeneInfo()

query = mg.getgenes(genes, scopes=["refseq"], fields="symbol", species="human")


df = pd.DataFrame(columns=["RefSeq", "Entrez", "Symbol"])

query = [dict(x) for x in set(sum([l.items() for l in query]))]

for i, q, in enumerate(query):
    df.set_value(i, "RefSeq", q["query"])
    df.set_value(i, "Entrez", q["_id"])
    df.set_value(i, "Symbol", q["symbol"])


writer = pd.ExcelWriter(os.path.join(data_dir, "mcf10_hsf1_tss_gene_intersection.xlsx"))
df.to_excel(writer)
writer.save()


