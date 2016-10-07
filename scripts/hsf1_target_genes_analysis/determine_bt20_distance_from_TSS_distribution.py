import pandas as pd
import os
import pybedtools as pbt
import mygene
import matplotlib.pyplot as plt

data_dir = os.path.join("..", "..", "data_files")
hg18_tss_file = os.path.join(data_dir, "hg18_TSS.ranged.2500.bed")
hg18_exact_tss_file = os.path.join(data_dir, "hg18_TSS.ranged.0.bed")
bt20_data_file = os.path.join(data_dir, "mendillo_merged", "BT20_Regions_HSF1.bed")

hg18_bed = pbt.BedTool(hg18_tss_file)
hg18_bed_exact = pbt.BedTool(hg18_exact_tss_file)
bt20_bed = pbt.BedTool(bt20_data_file)

#we may also want to make a bed file for bt20 based solely on the peaks



print "overall intersection: ", len(hg18_bed.intersect(bt20_bed))
print "number of annotations in bt20: ", len(bt20_bed)

close_genes = hg18_bed_exact.sort().closest(bt20_bed, D="a",t="first")
# print close_genes[0]
# print close_genes[0][-1]
# print len(close_genes)
under_2500 = []

for gene in close_genes:
    if abs(int(gene[-1])) <= 2500:
        under_2500.append(gene.name)


print "number of genes under 2500: ", len(under_2500)
print under_2500[0]

out_file = os.path.join(data_dir, "bt20_hg18_within_2500.xlsx")
if not os.path.isfile(out_file):
    mg = mygene.MyGeneInfo()

    query = mg.querymany(under_2500, scopes=["refseq"], fields="symbol, ensembl.gene", returnall=True)

    unfound = query["missing"]

    df = pd.DataFrame(index=under_2500, columns=["Ensembl", "Symbol"])

    for find in query["out"]:
        if "notfound" in find: continue
        if isinstance(find, list):
            find = find[0]

        try:
            if isinstance(find["ensembl"], list):
                find["ensembl"] = find["ensembl"][0]

            df.set_value(find["query"], "Ensembl", find["ensembl"]["gene"])
        except KeyError:
            df.set_value(find["query"], "Ensembl", None)

        try:
            if isinstance(find["symbol"], list):
                find["symbol"] = find["symbol"][0]

            df.set_value(find["query"], "Symbol", find["symbol"])
        except:
            df.set_value(find["query"], "Symbol", None)


    writer = pd.ExcelWriter(out_file)
    df.to_excel(writer)
    writer.save()

under_2500_bed = close_genes.filter(lambda x: abs(int(x[-1])) <= 2500).saveas()


print under_2500_bed[0]
print len(under_2500_bed), " <- should be approx 3317"

distances = []
peak_tss_distances = []
for gene in under_2500_bed:
    if int(gene[-1]) == -1:continue
    # print gene
    distances.append(int(gene[-1]))
    peak_tss_distances.append((int(gene[7]) + int(gene[15]) ) - int(gene[1]))

# print distances
print "good distances: ", len(peak_tss_distances)
fig = plt.figure()
plt.hist(peak_tss_distances, bins=100)
plt.title("Histogram of HSF1 peak to TSS distances in BT20 cells (hg18 genome)")
plt.xlabel("Distance from TSS")
plt.ylabel("Frequency")
fig.savefig("hsf1_peak_tss_distance_bt20.pdf")
plt.show()
