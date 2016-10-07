import os

data_dir = os.path.join("..", "data_files")
data_file = os.path.join(data_dir, "human_encode_cell_types.txt")

f = open(data_file)
ctypes = f.readlines()
f.close()

epi = []
for line in ctypes:
    if len(line.split()) > 1 and "+" in line:
        epi.append(line.split()[0])


outfile = os.path.join(data_dir, "human_encode_epithelial_cell_types.txt")
f = open(outfile, "w")
f.write("\n".join(epi))
f.close()