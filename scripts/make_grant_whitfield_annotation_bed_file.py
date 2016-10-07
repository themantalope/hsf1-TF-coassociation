import pandas as pd
import os
import pybedtools as pbt

data_dir = os.path.join("..", "data_files")
annotation_file = os.path.join(data_dir, "grant_whitfield_probe_annotations.xlsx")

annotation = pd.read_excel(open(annotation_file))

chr_col = "CHROMOSOMAL_LOCATION"
annotation_col = "ACCESSION_STRING"

bed_fields = []
for row in annotation.iterrows():
    idx = row[0]
    rowdata = row[1]
    pos = rowdata[chr_col]
    name = rowdata[annotation_col]
    print rowdata[chr_col]
    if pd.isnull(name): continue
    elif "ref" not in name: continue
    elif "chr" not in pos: continue
    chr = pos[0:pos.find(":")]
    start, stop = pos[pos.find(":")+1:].split("-")
    bed_fields.append([chr, start, stop, name])

annobed = pbt.BedTool(bed_fields).sort()
outf = os.path.join(data_dir, "grant_whitfield_bed.bed")
annobed.saveas(outf)

