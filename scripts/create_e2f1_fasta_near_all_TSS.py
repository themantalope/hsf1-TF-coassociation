import os
import pybedtools as pbt
import requests
import json
import tqdm
import sys

data_dir = os.path.join("..", "data_files")
e2f1_bed_file = os.path.join(data_dir, "cc_tf_bedfiles", "E2F1", "E2F1.bed")
hg19_TSS_bed_file = os.path.join(data_dir, "hg19_TSS.ranged.0.unique.bed")

base_url = "http://rest.ensembl.org/sequence/region/human/"
outfasta = os.path.join(data_dir, "e2f1_near_tss_regions_400.fasta")

hg19_TSS = pbt.BedTool(hg19_TSS_bed_file).sort()
e2f1_bed = pbt.BedTool(e2f1_bed_file).sort()

good_features = []
uniques = []

regions_as_peaks = []

for feature in e2f1_bed:

    #make a temp bedfile

    peak = int(feature[-1])
    peak_pos = int(feature[1]) + peak
    newstart = newend = peak_pos

    fields = feature.fields
    fields[1] = fields[2] = str(newstart)
    fields[3] = ":".join(feature.fields)
    regions_as_peaks.append(fields)

temp_bed = pbt.BedTool(regions_as_peaks).sort()

close = hg19_TSS.sort().closest(temp_bed.sort(), D="a")

# print close[0].fields
# raw_input()


for i, feature in enumerate(close):
    dist = abs(int(feature[-1]))
    if dist <= 2000 and int(feature[-1]) != -1:
        if feature.fields[9] not in uniques:
            good_features.append(feature.fields[6:10])
            uniques.append(feature.fields[9])



good_bed = pbt.BedTool(good_features).sort()

# print good_bed[0].fields
# raw_input()

package = []
package_len = 50
package_names = []

with open(outfasta, "w") as of:
    for i, feature in enumerate(tqdm.tqdm(good_bed)):
        name = feature.name
        name_parts = name.split(":")
        peak_pos = int(name_parts[-1]) + int(name_parts[1])
        reg_start = str(peak_pos - 200)
        reg_end = str(peak_pos + 200)
        region = name_parts[0]+":"+reg_start+".."+reg_end+ ":1"
        # print region
        # raw_input()
        package.append(region)
        package_names.append(feature.name)
        if len(package) == package_len or i == len(good_bed) - 1:
            post = requests.post(base_url+"?coord_system_version=GRCh37", data=json.dumps({"regions":package}), headers={"content-type":"application/json"})
            seqs = post.json()
            # print seqs
            for j, seq in enumerate(seqs):
                # print seq
                # tqdm.tqdm.write(seq, file=sys.stdout)
                of.write(">"+package_names[j]+"\n")
                of.write(seq["seq"]+"\n")
            package = []
            package_names = []



