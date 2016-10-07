import mygene
import pandas as pd
import os
import pickle
import json
# import pprint

data_dir = os.path.join("..", "data_files")
mendillo_HSF1_gene_file = os.path.join(data_dir, "mendillo_hsf1_binding_genes.xlsx")

mendillo_HSF1_gene = pd.read_excel(open(mendillo_HSF1_gene_file))

good_genes = mendillo_HSF1_gene

gene_list = good_genes.index.tolist()

mg = mygene.MyGeneInfo()

search = mg.querymany(gene_list, scopes=["refseq"], fields="ensembl.gene", species="human", returnall=True)

# pprint.pprint(search)
# raw_input()


bad_sym = search["missing"]

mapping = {}

for find in search["out"]:
    if "notfound" in find:continue
    if isinstance(find, list):
        find = find[0]

    if "ensembl" not in find:
        bad_sym.append(find["query"])
        continue

    if isinstance(find["ensembl"], list):
        find["ensembl"] = find["ensembl"][0]

    mapping[find["query"]] = find["ensembl"]["gene"]

f = open(os.path.join(data_dir, "unfound_hsf1_ids.txt"), "w")
f.write("\n".join(bad_sym))
f.close()

pickle.dump(mapping, open(os.path.join(data_dir, "hsf1_mendillo_genes_map.pickle"), "w"))
json.dump(mapping, open(os.path.join(data_dir, "hsf1_mendillo_genes_map.json"), "w"))