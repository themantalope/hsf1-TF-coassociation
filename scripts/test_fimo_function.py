import chrom_utils
import os

data_dir = os.path.join("..", "data_files")
fasta_dir = os.path.join(data_dir, "hsf1_merged_regions_fasta")
lefile = os.path.join(fasta_dir, "hsf1_0.05_merged_regions.fasta")

df = chrom_utils.fimo_search_for_motif(lefile)
print df