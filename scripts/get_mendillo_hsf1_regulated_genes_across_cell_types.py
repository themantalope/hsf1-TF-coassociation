import pandas as pd
import os


data_dir = os.path.join("..", "data_files")
mendillo_file = os.path.join(data_dir, "mendillo_hsf1_binding.xlsx")

mendillo = pd.read_excel(open(mendillo_file), sheetname=None)

gene_sheets = [k for k in mendillo.keys() if ("GENES" in k or "Gene" in k) and "HEATSHOCK" not in k]

# print gene_sheets


all_genes = set()


for s in gene_sheets:
    all_genes.update(mendillo[s]["ID1"])

outdf = pd.DataFrame(index=list(all_genes), columns=gene_sheets)

for s in gene_sheets:
    curdf = mendillo[s]
    genes = curdf["ID1"].values.tolist()
    for g in genes:
        outdf.set_value(index=g,col=s,value=1.0)


outdf["sum"] = outdf.sum(axis=1)
outdf["frac"] = outdf["sum"] / float(len(gene_sheets))

writer = pd.ExcelWriter(os.path.join(data_dir, "mendillo_hsf1_binding_genes.xlsx"))

outdf.to_excel(writer); writer.save()