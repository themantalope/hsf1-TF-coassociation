import os
import pybedtools as pbt
import requests
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
from sklearn import metrics
import chrom_utils

hse_patt = re.compile("(?P<hse>GAA[A-Z]{2}TTC[A-Z]{2}GAA)")

data_dir = os.path.join("..", "data_files")
hsf1_merged_region_dir = os.path.join(data_dir, "hsf1_merged_region_files")

base_url = "http://rest.ensembl.org/sequence/region/human/"

hsf1_merged_region_fasta_dir = os.path.join(data_dir, "hsf1_merged_regions_fasta")
if not os.path.isdir(hsf1_merged_region_fasta_dir):
    os.mkdir(hsf1_merged_region_fasta_dir)

hsf1_merged_regions_positive_files = [os.path.join(hsf1_merged_region_dir, f) for f in os.listdir(hsf1_merged_region_dir) if "negatives" not in f and "bed" in f]
hsf1_merged_regions_negative_files = [os.path.join(hsf1_merged_region_dir, f) for f in os.listdir(hsf1_merged_region_dir) if "negatives" in f and "bed" in f]

#make fasta file for each cutoff

pack_size = 45

prevfile = None
prevlines = []
for f in hsf1_merged_regions_positive_files:
    outfasta = os.path.join(hsf1_merged_region_fasta_dir, os.path.basename(f).replace(".bed", ".fasta"))
    print "getting data for: {f}".format(f=os.path.basename(outfasta))
    if os.path.exists(outfasta):continue
    else:
        hsf1_pos_bed = pbt.BedTool(f)

        if prevfile is not None:
            prevlines = open(prevfile, "r").readlines()

        with open(outfasta, "w") as of:
            package = []
            features = []
            for idx, feature in enumerate(hsf1_pos_bed):

                found = False
                for i in range(0, len(prevlines), 2):
                    if feature.name in prevlines[i]:
                        # print i
                        of.write(prevlines[i])
                        of.write(prevlines[i + 1])
                        found = True
                        break

                if found:
                    continue

                region = feature[0]+":"+feature[1]+".."+feature[2]+":1"
                package.append(region)
                features.append(feature.name)
                if len(package) == pack_size or idx == len(hsf1_pos_bed)-1:
                    req = requests.post(url=base_url+"?coord_system_version=NCBI36", data=json.dumps({"regions":package}), headers={"content-type":"application/json"})
                    seqs = req.json()
                    for i, seq in enumerate(seqs):
                        if "N" in seq["seq"]:continue
                        of.write(">"+features[i]+"\n")
                        of.write(seq["seq"]+"\n")
                    package = []
                    features = []
        prevfile = outfasta

prevfile = None
prevlines = []
for f in hsf1_merged_regions_negative_files:
    outfasta = os.path.join(hsf1_merged_region_fasta_dir, os.path.basename(f).replace(".bed", ".fasta"))
    print "getting data for: {f}".format(f=os.path.basename(outfasta))
    if os.path.exists(outfasta):continue
    else:
        hsf1_pos_bed = pbt.BedTool(f)

        if prevfile is not None:
            prevlines = open(prevfile, "r").readlines()
            print len(prevlines)

        with open(outfasta, "w") as of:
            package = []
            features = []
            for idx, feature in enumerate(hsf1_pos_bed):

                found = False
                for i in range(0, len(prevlines), 2):
                    if feature.name in prevlines[i]:
                        # print i
                        of.write(prevlines[i])
                        of.write(prevlines[i+1])
                        found = True
                        break

                if found:
                    continue

                region = feature[0]+":"+feature[1]+".."+feature[2]+":1"
                package.append(region)
                features.append(feature.name)
                if len(package) == pack_size or idx == len(hsf1_pos_bed)-1:
                    req = requests.post(url=base_url+"?coord_system_version=NCBI36", data=json.dumps({"regions":package}), headers={"content-type":"application/json"})
                    seqs = req.json()
                    for i, seq in enumerate(seqs):
                        if "N" in seq["seq"]: continue
                        of.write(">"+features[i]+"\n")
                        of.write(seq["seq"]+"\n")
                    package = []
                    features = []
        prevfile = outfasta

#ok we have the fasta files, now for each cutoff find the number of TP, FP, FN, TN





pos_fastas = [os.path.join(hsf1_merged_region_fasta_dir, f) for f in os.listdir(hsf1_merged_region_fasta_dir) if "negatives" not in f]
file_pairs = [(f, f.replace(".fasta", ".negatives.fasta")) for f in pos_fastas]

cutoff_patt = re.compile("(?<=\_)(?P<ct>\d+\.\d+)(?=\_)")


cts = []
for fp in file_pairs:
    cts.append(re.findall(cutoff_patt,fp[0])[0])

roc_df = pd.DataFrame(index=cts, columns=["TP", "FP", "FN", "TN", "sensitivity", "specificity", "1-specificity"])

for fp in file_pairs:
    pos = fp[0]
    neg = fp[1]
    cutoff = re.findall(cutoff_patt,pos)[0]
    pos_features = set([line.replace(">", "") for line in open(pos, "r").readlines() if ">" in line])
    neg_features = set([line.replace(">", "") for line in open(neg, "r").readlines() if ">" in line])
    #
    # n_TP = 0
    # for seq in pos_seqs:
    #     if re.findall(hse_patt, seq): n_TP += 1
    # n_FP = len(pos_seqs) - n_TP
    #
    # n_TN = 0
    # for seq in neg_seqs:
    #     if not re.findall(hse_patt, seq): n_TN += 1
    # n_FN = len(neg_seqs) - n_TN

    if len(pos_features) > 0:
        pos_fimo_df = chrom_utils.fimo_search_for_motif(pos, thresh=1e-4)
        n_TP = len(set(pos_fimo_df["sequence name"].values.tolist()))
        n_FP = len(pos_features) - n_TP
    else:
        n_TP = n_FP = 0

    if len(neg_features) > 0:
        neg_fimo_df = chrom_utils.fimo_search_for_motif(neg, thresh=1e-4)
        n_FN = len(set(neg_fimo_df["sequence name"].values.tolist()))
        n_TN = len(neg_features) - n_FN
    else:
        n_TN = n_FN = 0


    print "cutoff: ", cutoff
    print "TP: ", n_TP
    print "FP: ", n_FP
    print "FN: ", n_FN
    print "TN: ", n_TN
    print "sens: ", float(n_TP) / (float(n_TP) + float(n_FN))
    print "spe: ", float(n_TN) / (float(n_TN) + float(n_FP))
    # raw_input()


    roc_df.set_value(index=cutoff, col="TP", value=n_TP)
    roc_df.set_value(index=cutoff, col="FP", value=n_FP)
    roc_df.set_value(index=cutoff, col="FN", value=n_FN)
    roc_df.set_value(index=cutoff, col="TN", value=n_TN)
    roc_df.set_value(index=cutoff, col="sensitivity", value=float(n_TP) / (float(n_TP) + float(n_FN)))
    roc_df.set_value(index=cutoff, col="specificity", value=float(n_TN) / (float(n_TN) + float(n_FP)))
    roc_df.set_value(index=cutoff, col="1-specificity", value=1.0 - float(n_TN) / (float(n_TN) + float(n_FP)))


roc_outfile = os.path.join(data_dir, "hsf1_regions_hse_ROC_analysis.xlsx")
writer = pd.ExcelWriter(roc_outfile)
roc_df.to_excel(writer)
writer.save()

#now make a plot of the ROC

print
print

fig = plt.figure()
ax = plt.gca()
x = roc_df["1-specificity"].values
print "x: ", x
y = roc_df["sensitivity"].values
print "y: ", y

auc = metrics.auc(x=x,y=y,reorder=True)
print "auc: ", auc

labels = [str(c) for c in roc_df.index.tolist()]
print labels
# raw_input()
plt.plot(x, y)
plt.title("ROC for HSF1 Binding sites")
plt.xlabel("1-specificity (False Positive Rate)")
plt.ylabel("sensitivity")
for i, l in enumerate(labels):
    plt.annotate("{l:.2f}".format(l=float(l)), (x[i], y[i]))

plt.annotate("AUC: {a:.2f}".format(a=auc), (0.8, 0.4))
fig.savefig(os.path.join(data_dir, "hsf1_binding_ROC_relaxed_HSE.pdf"))

plt.show()