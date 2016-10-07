import os
import pybedtools as pbt
import requests
import json

data_dir = os.path.join("..", "data_files")
e2f1_hsf1_bed_file = os.path.join(data_dir, "hsf1_cctf_region_comparisons", "hsf1_E2F1_le_500.bed")
base_url = "http://rest.ensembl.org/sequence/region/human/"
outfasta = os.path.join(data_dir, "e2f1_near_hsf1_regions.fasta")

e2f1_hsf1_bed = pbt.BedTool(e2f1_hsf1_bed_file)
e2f1_regions = []
for i, feature in enumerate(e2f1_hsf1_bed):
    e2f1_regions.append(feature.fields[4:7] + ["e2f1_feature{n}".format(n=i)+ "|"  +":".join(feature.fields[4:7])])

e2f1_bed = pbt.BedTool(e2f1_regions).sort()

# print e2f1_bed[0].fields

package = []
pack_names = []
pack_length = 50


with open(outfasta, "w") as of:
    for i, feature in enumerate(e2f1_bed):
        region = feature[0]+":"+feature[1]+".."+feature[2]+":1"
        package.append(region)
        pack_names.append(feature.name)
        if len(package) == pack_length or i == len(e2f1_bed)-1:
            lepost = requests.post(base_url+"?coord_system_version=GRCh37", data=json.dumps({"regions":package}), headers={"content-type":"application/json"})
            seqs = lepost.json()
            for j, seq in enumerate(seqs):
                of.write(">"+pack_names[j]+"\n")
                of.write(seq["seq"]+"\n")
            package = []
            pack_names = []

