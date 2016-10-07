import os
import pandas as pd

data_dir = os.path.join("..", "data_files")
dominguez_file = os.path.join(data_dir, "dominguez_wang_periodic_expression.xlsx")
data = pd.read_excel(open(dominguez_file))

data["Cell Cycle Label"] = data.apply(lambda x: "G1-S" if (x["Seed match"] == 1 or x["Seed match"] == 2) else "G2-M", axis=1)

writer = pd.ExcelWriter(os.path.join(data_dir, "dominguez_wang_periodic_expression_modified.xlsx"))
data.to_excel(writer)
writer.save()



