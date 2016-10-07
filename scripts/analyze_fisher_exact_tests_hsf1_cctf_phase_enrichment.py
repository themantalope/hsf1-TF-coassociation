import pandas as pd
import os
import scipy.stats as st
import pprint

data_dir = os.path.join("..", "data_files")
hsf1_cctf_position_file = os.path.join(data_dir, "hsf1_cctf_position_analysis.xlsx")

hsf1_cctf_position = pd.read_excel(open(hsf1_cctf_position_file), sheetname=None)

e2f1 = hsf1_cctf_position["E2F1"]
foxm1 = hsf1_cctf_position["FOXM1"]

# overlapping = set(e2f1.index.tolist()).intersection(set(foxm1.index.tolist()))
overlapping = set()


e2f1 = e2f1[~(e2f1.index.isin(overlapping))]
foxm1 = foxm1[~(foxm1.index.isin(overlapping))]

e2f1_total_genes = len(e2f1[~(e2f1["Phase"].isin(["G1/S", "S", "G1"]))].values.tolist())
e2f1_g1s_genes = len(e2f1[e2f1["Phase"].isin(["G1/S", "S", "G1"])].values.tolist())

foxm1_total_genes = len(foxm1[~(foxm1["Phase"].isin(["G1/S", "S", "G1"]))].values.tolist())
foxm1_g1s_genes = len(foxm1[foxm1["Phase"].isin(["G1/S", "S", "G1"])].values.tolist())

table = [[e2f1_g1s_genes, foxm1_g1s_genes], [e2f1_total_genes, foxm1_total_genes]]

pprint.pprint(table)

oddsr, pv = st.fisher_exact(table)


print "overlapping genes: ", len(set(e2f1.index.tolist()).intersection(set(foxm1.index.tolist())))
print "odds ratio: ", oddsr
print "pv: ", pv