import os
import pandas as pd
import pickle

data_dir = os.path.join("..", "data_files")
output_file = os.path.join(data_dir, "initial_output.pickle")

data = pickle.load(open(output_file))

rows = []
for x in data:
    if x[0] in rows:
        continue
    else:
        rows.append(x[0])

cols = []
for x in data:
    if x[1] in cols:
        continue
    else:
        cols.append(x[1])


df = pd.DataFrame(index=rows, columns=cols)

for x in data:
    df.set_value(x[0], x[1], x[2])

writer = pd.ExcelWriter(os.path.join(data_dir, "initial_HSF1_output_organized.xlsx"))
df.to_excel(writer)
writer.save()
