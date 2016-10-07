import os
import pandas as pd

data_dir = os.path.join("..", "data_files")
grant_whitfield_file = os.path.join(data_dir, "grant_whitfield_cell_cycle_expression.xlsx")
agilent_file = os.path.join(data_dir, "agilent_probe_mapping.xlsx")


agilent = pd.read_excel(open(agilent_file))
grant_whitfield = pd.read_excel(open(grant_whitfield_file),sheetname=None)

cc1 = grant_whitfield["CC1"]
# print cc1.columns.tolist()
# print cc1.iloc[:, 0]

agilent_probe_ids = cc1.iloc[:, 0].values.tolist()
outids = []
for aid in agilent_probe_ids:
    outids.append((aid.replace("AGI_HUM1_OLIGO_", ""), aid))

print outids[0:10]

relevant_ids = agilent[agilent["ProbeID"].isin([x[0] for x in outids])]
print len(outids)
print len(relevant_ids)

relevant_ids["GW_ProbeID"] = relevant_ids.apply(lambda x: next((z[1] for z in outids if z[0] == x["ProbeID"])), axis=1)
print relevant_ids

grant_whitfield["ID_Mapping"] = relevant_ids

writer = pd.ExcelWriter(grant_whitfield_file)
for k in grant_whitfield:
    grant_whitfield[k].to_excel(writer, sheet_name=k)

writer.save()
