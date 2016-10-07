import os
import pybedtools as pbt
import chrom_utils
import math
import requests
import json
import tempfile

data_dir = os.path.join("..", "data_files")
hsf1_bed_file = os.path.join(data_dir, "hsf1_merged_region_files", "hsf1_0.166666666667_merged_regions.bed")
hsf1_fasta_file = os.path.join(data_dir, "hsf1_merged_regions_fasta", "hsf1_0.166666666667_merged_regions.fasta")

hsf1_bed = pbt.BedTool(hsf1_bed_file)

narrow_features = []

hse_df = chrom_utils.fimo_search_for_motif(hsf1_fasta_file)

# print hse_df

hse_feature_names = hse_df["sequence name"].values.tolist()

good_initial_regions = []

for feature in hsf1_bed:
    if abs( int(feature[1]) - int(feature[2]) ) < 200 and feature.name not in hse_feature_names:
        narrow_features.append(feature.fields)
    else:
        good_initial_regions.append(feature.fields)

print "number narrow: ", len(narrow_features)
print "number in hsf1 bed: ", len(hsf1_bed)
print "number already with hse: ", len(list(set(hse_feature_names)))
print "number of good_initial: ", len(good_initial_regions)

#make a new list of feature fields extending with width of each field
for i,field in enumerate(narrow_features):
    # print "before: ", narrow_features[i]
    start = int(field[1])
    end = int(field[2])
    middle = int(math.floor(  (float(start) + float(end)) / 2.0  ))
    newstart = middle - 200
    newend = middle + 200
    field[1] = str(newstart)
    field[2] = str(newend)
    name = field[3]
    fn = name[0:name.find("|")]
    newname = fn+"|"+field[0]+"-"+field[1]+"-"+field[2]
    field[3] = newname
    # print "after: ", narrow_features[i]
    # raw_input()

newbed = pbt.BedTool(narrow_features)

base_url = "http://rest.ensembl.org/sequence/region/human/"
fasta_str = ""
package = []
names = []
plength = 50

for i, feature in enumerate(newbed):
    regions = feature[0]+":"+feature[1]+".."+feature[2]+":1"
    package.append(regions)
    names.append(feature[3])
    if len(package) == 50 or i==len(newbed)-1:
        req = requests.post(base_url+"?coord_system_version=NCBI36", data=json.dumps({"regions":package}),headers={"content-type":"application/json"})
        for j, seq in enumerate(req.json()):
            fasta_str+=">"+ names[j]+"\n"
            fasta_str+= seq["seq"]+"\n"
        package = []
        names = []


# print fasta_str

tf = tempfile.NamedTemporaryFile()
tf.write(fasta_str)

newbed_hse = chrom_utils.fimo_search_for_motif(tf.name)
# print newbed_hse
# print len(newbed)

output_hse_bed = []
good_features = list(set(hse_feature_names)) + list(set(newbed_hse["sequence name"].values.tolist()))

for feature in good_features:
    region = feature[feature.find("|")+1:]
    parts = region.split("-")
    chr = parts[0]
    start = parts[1]
    end = parts[2]
    name = feature
    output_hse_bed.append([chr, start, end, name])

output_hse_bed = pbt.BedTool(output_hse_bed)

output_file_name = os.path.join(data_dir, "hsf1_0.17_merged_regions_hse.bed")
output_hse_bed.saveas(output_file_name)

all_regions_output_file = os.path.join(data_dir, "hsf1_0.17_merged_regions.bed")

all_features = good_initial_regions + [feature.fields for feature in newbed]

for feature in all_features:
    start = int(feature[1])
    end = int(feature[2])
    if abs(start - end) < 400:
        middle = int(math.floor(float(start + end) / 2.0))
        newstart = middle - 200
        newend = middle + 200
        feature[1] = str(newstart)
        feature[2] = str(newend)
        feature[3] = feature[3][0:feature[3].find("|")] + "|" + feature[0]+"-"+feature[1]+"-"+feature[2]+"-1"



for feature in all_features:
    assert abs(int(feature[1]) - int(feature[2])) >= 400, feature

# for feature in all_features:
#     assert abs(int(feature[1]) - int(feature[2])) >= 400, feature

print len(all_features), len(hsf1_bed)

all_features_bed = pbt.BedTool(all_features)

assert len(set([feature.name for feature in all_features_bed])) == len(hsf1_bed)

all_features_bed.sort()
all_features_bed.saveas(all_regions_output_file)