import re
import os

data_dir = os.path.join("..", "data_files", "encode", "narrow_peaks", "human")

cell_pattern = re.compile("(?<=Tfbs)([A-Z-0-9-a-z].*?)(?=[A-Z])")

all_data_files = [f for f in os.listdir(data_dir) if "Histone" not in f and "Tfbs" in f]

cell_types = set()

for f in all_data_files:
    matches = re.findall(cell_pattern, f)
    cell_types.add(matches[0])

# print cell_types

output_file = os.path.join("..", "data_files", "human_encode_cell_types.txt")
with open(output_file, "w") as f:
    f.write("\n".join(cell_types))