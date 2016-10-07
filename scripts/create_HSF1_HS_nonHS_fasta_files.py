import os
import pybedtools as pbt
import requests
import json
import tqdm

data_dir = os.path.join("..", "data_files")
hsf1_HS_exclusive_file = os.path.join(data_dir, "hsf1_HS_exclusive.bed")
hsf1_nonHS_exclusive_file = os.path.join(data_dir, "hsf1_nonHS_exclusive.bed")

hsf1_HS_exclusive = pbt.BedTool(hsf1_HS_exclusive_file)
hsf1_nonHS_exclusive = pbt.BedTool(hsf1_nonHS_exclusive_file)

out_HS_file = os.path.join(data_dir, "hsf1_HS_exclusive_regions.fasta")
out_nonHS_file = os.path.join(data_dir, "hsf1_nonHS_exclusive_regions.fasta")

base_url = "http://rest.ensembl.org/sequence/region/human/"

package = []
package_len = 50
package_names = []

with open(out_HS_file, "w") as of:
    for i, feature in enumerate(tqdm.tqdm(hsf1_HS_exclusive)):
        region = feature.fields[0]+":"+feature.fields[1]+".."+feature.fields[2]+":1"
        package.append(region)
        package_names.append(feature.name)

        if len(package) == package_len or i == len(hsf1_HS_exclusive) - 1:
            post = requests.post(base_url+"?coord_system_version=GRCh37", data=json.dumps({"regions":package}), headers={"content-type":"application/json"})
            seqs = post.json()
            for j, seq in enumerate(seqs):
                of.write(">"+package_names[j]+"\n")
                of.write(seq["seq"]+"\n")
            package = []
            package_names = []


with open(out_nonHS_file, "w") as of:
    for i, feature in enumerate(tqdm.tqdm(hsf1_nonHS_exclusive)):
        region = feature.fields[0]+":"+feature.fields[1]+".."+feature.fields[2]+":1"
        package.append(region)
        package_names.append(feature.name)

        if len(package) == package_len or i == len(hsf1_nonHS_exclusive) - 1:
            post = requests.post(base_url+"?coord_system_version=GRCh37", data=json.dumps({"regions":package}), headers={"content-type":"application/json"})
            seqs = post.json()
            for j, seq in enumerate(seqs):
                of.write(">"+package_names[j]+"\n")
                of.write(seq["seq"]+"\n")
            package = []
            package_names = []