import pybedtools as pbt
import re

# cell_type_pattern = re.compile("(?<=Tfbs)([A-Z-0-9-a-z].*?)(?=[A-Z])")
cell_type_tf_pattern = re.compile("(?<=Tfbs)(?P<cell>[A-Z-0-9-a-z].*?)(?=[A-Z])(?P<tf>[A-Z-0-9-a-z].*?)(?=[A-Z])")

def get_cell_type_from_encode_filename(fn):
    return re.findall(cell_type_tf_pattern, fn)[0][0]

def compute_tss_jaccard(tss_bed, a_bed, b_bed):
    a_tss = a_bed.intersect(tss_bed)
    b_tss = b_bed.intersect(b_bed)
    x = pbt.BedTool()
    jd = x.jaccard(a=a_tss, b=b_tss)
    return jd["jaccard"]

def get_tf_type_from_encode_filename(fn):
    return re.findall(cell_type_tf_pattern, fn)[0][1]

