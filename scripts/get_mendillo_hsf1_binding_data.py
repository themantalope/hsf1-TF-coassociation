import os
import pandas as pd
import pybedtools as pbt
import chrom_utils

data_dir = os.path.join("..", "data_files")
mendillo_output_dir = os.path.join(data_dir, "mendillo")
if not os.path.isdir(os.path.join(data_dir, "mendillo")):
    os.mkdir(os.path.join(data_dir, "mendillo"))

data_file = os.path.join(data_dir, "mendillo_hsf1_binding.xlsx")

mendillo_data = pd.read_excel(open(data_file), sheetname=None)

good_sheets = [k for k in mendillo_data.keys() if "HEATSHOCK" not in k and "Regions" in k]
heatshock_sheets = [k for k in mendillo_data.keys() if "HEATSHOCK" in k and "Regions" in k]
# print good_sheets


#narrow peak format: chr start stop name score strand signalValue pValue qValue peak



for s in good_sheets:
    sheet = mendillo_data[s]
    features = []
    for row in sheet.iterrows():
        idx = row[0]
        rowdata = row[1]
        chr_num = rowdata["CHROM"]

        if chrom_utils.is_number(chr_num): chromosome = "chr"+str(int(rowdata["CHROM"]))
        else: chromosome = "chr"+str(rowdata["CHROM"])


        features.append((chromosome, int(rowdata["START"]), int(rowdata["END"]), ".", "0", ".", int(rowdata["TOTAL TARGET COUNTS"]), "-1", "-1", int(rowdata["PEAK POS"] - rowdata["START"])))

    newbed = pbt.BedTool(features)
    newbed.saveas(os.path.join(mendillo_output_dir, s+"_HSF1.bed"))

mendillo_output_hs_dir = os.path.join(data_dir, "mendillo_HS")

if not os.path.isdir(mendillo_output_hs_dir):
    os.mkdir(mendillo_output_hs_dir)

for s in heatshock_sheets:
    sheet = mendillo_data[s]
    features = []
    for row in sheet.iterrows():
        idx = row[0]
        rowdata = row[1]
        chr_num = rowdata["CHROM"]

        if chrom_utils.is_number(chr_num): chromosome = "chr"+str(int(rowdata["CHROM"]))
        else: chromosome = "chr"+str(rowdata["CHROM"])


        features.append((chromosome, int(rowdata["START"]), int(rowdata["END"]), ".", "0", ".", int(rowdata["TOTAL TARGET COUNTS"]), "-1", "-1", int(rowdata["PEAK POS"] - rowdata["START"])))

    newbed = pbt.BedTool(features)
    newbed.saveas(os.path.join(mendillo_output_hs_dir, s+"_HS_HSF1.bed"))