import pybedtools as pbt
import re
import os
import sys
import pandas as pd
import multiprocessing as mp
import subprocess
import StringIO
import scipy.stats as st
import numpy as np
import inspect
import numbers
import math

# cell_type_pattern = re.compile("(?<=Tfbs)([A-Z-0-9-a-z].*?)(?=[A-Z])")
cell_type_tf_exp_pattern = re.compile("(?<=Tfbs)(?P<cell>[A-Z-0-9-a-z].*?)(?=[A-Z])(?P<tf>[A-Z-0-9-a-z].*?)(?=[A-Z])(?P<exp>[A-Z-0-9-a-z].*?)(?=\.bam)")

mendillo_tf_pattern = re.compile("(?<=Regions\_)(.*?)(?=\.)")
mendillo_cell_pattern = re.compile("(?<=^)(.*?)(?=\_Regions)")

chr_start_pattern = re.compile("(?<=^)(chr\w+)(?=$)")

dir_path = os.path.dirname(os.path.realpath(__file__))
chrom_size_file = os.path.join(dir_path, "hg19.chromsizes")

if not os.path.isfile(chrom_size_file):
    print "got dir: ", dir_path
    sys.exit(0)


g = pbt.genome_registry.hg19


hg19_unique_tss_length = 27090


hg19_500_bp_fragments_near_tss = 228109

def compute_tss_jaccard(tss_bed, a_bed, b_bed):
    a_tss = a_bed.intersect(tss_bed)
    b_tss = b_bed.intersect(b_bed)
    x = pbt.BedTool()
    jd = x.jaccard(a=a_tss.sort(), b=b_tss.sort())
    return jd["jaccard"]

def compute_tss_jaccard_intersection(tss_bed, a_bed, b_bed):
    a_tss = a_bed.intersect(tss_bed)
    b_tss = b_bed.intersect(b_bed)
    a_b_tss = a_tss.intersect(b_tss)
    return float(len(a_b_tss)) / float(len(a_tss) + len(b_tss) - len(a_b_tss))


def get_exp_type_from_filename(fn):
    if "Encode" in fn:
        return re.findall(cell_type_tf_exp_pattern, fn)[0][2]
    else:
        return ""


def get_cell_type_from_filename(fn):
    if "Encode" in fn:
        return get_cell_type_from_encode_filename(fn)
    else:
        return get_cell_type_from_mendillo_filename(fn)


def get_tf_type_from_filename(fn):
    if "Encode" in fn:
        return get_tf_type_from_encode_filename(fn)
    else:
        return get_tf_type_from_mendillo_filename(fn)


def get_cell_type_from_encode_filename(fn):
    return re.findall(cell_type_tf_exp_pattern, fn)[0][0]

def get_tf_type_from_encode_filename(fn):
    return re.findall(cell_type_tf_exp_pattern, fn)[0][1]

def get_cell_type_from_mendillo_filename(fn):
    return re.findall(mendillo_cell_pattern, fn)[0]

def get_tf_type_from_mendillo_filename(fn):
    return re.findall(mendillo_tf_pattern, fn)[0]

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def calculate_TSS_co_binding_rate(tss_bed, a_bed, b_bed):
    """
    Calculates the degree of co-binding to a ranged tss between two binding profiles, contained in a_bed and b_bed.
    Make sure that the tss_bed file has a unique "name" attribute for each feature.
    :param tss_bed:
    :param a_bed:
    :param b_bed:
    :return:
    """
    tss_a_intersect = tss_bed.intersect(a_bed)
    tss_b_intersect = tss_bed.intersect(b_bed)
    a_tss_names = set([tss_a_intersect[n].name for n in range(len(tss_a_intersect))])
    b_tss_names = set([tss_b_intersect[n].name for n in range(len(tss_b_intersect))])
    return compute_jaccard_index(a_tss_names, b_tss_names)

# def calculate_TSS_co_binding_minimum_overlap_ratio(tss_bed, a_bed, b_bed):
#     t1 = timeit.default_timer()
#     tss_a_intersect = tss_bed.intersect(a_bed,u=True).sort()
#     tss_b_intersect = tss_bed.intersect(b_bed,u=True).sort()
#     t2 = timeit.default_timer()
#     print "intersect time: ", (t2-t1)
#     print len(tss_a_intersect)
#     print len(tss_b_intersect)
#
#     t1 = timeit.default_timer()
#     a_tss_names = set([tss_a_intersect[n].name for n in range(len(tss_a_intersect))])
#     b_tss_names = set([tss_b_intersect[n].name for n in range(len(tss_b_intersect))])
#     t2 = timeit.default_timer()
#     print "name set time: ", (t2-t1)
#     print len(a_tss_names)
#     print len(b_tss_names)
#     return compute_minimum_overlap_ratio(a_tss_names, b_tss_names)


# def calculate_TSS_co_binding_minimum_overlap_ratio(tss_bed, a_bed, b_bed,slop_width=5000):
#     a_slop_b_intersection_tss_intersection = tss_bed.intersect(a_bed.slop(b=slop_width, g=g).remove_invalid().intersect(b_bed))
#     print a_slop_b_intersection_tss_intersection.head()
#     print float(len(a_slop_b_intersection_tss_intersection))
#     print float(min(len(a_bed), len(b_bed)))
#     return float(len(a_slop_b_intersection_tss_intersection)) / float(min(len(a_bed), len(b_bed)))


def compute_minimum_overlap_ratio(set_a, set_b):
    intersect = float(len(set_a.intersection(set_b)))
    smaller = float(min(len(set_a), len(set_b)))
    return intersect / smaller

def compute_jaccard_index(set_a, set_b):
    n = float(len(set_a.intersection(set_b)))
    return n / float(len(set_a) + len(set_b) - n)


def set_region_start_stop_based_on_peak(bedfile, width=2500, ext=".peak.bed"):
    f= open(bedfile, "r")
    filelines = f.readlines()
    f.close()
    outlines = []
    for line in filelines:
        parts = line.split()
        peak = int(parts[9])
        newstart = str(int(parts[1]) + peak - width)
        newstop = str(int(parts[2]) + peak + width)
        parts[1] = newstart
        parts[2] = newstop
        outlines.append("\t".join(parts))

    dir = os.path.dirname(bedfile)
    base = os.path.basename(bedfile)
    _, bext = os.path.splitext(bedfile)

    outname = os.path.join(dir, base.replace(bext,ext))
    # print outname
    f = open(outname, "w")
    f.write("\n".join(outlines))
    f.close()


def calculate_TSS_co_binding_minimum_overlap_ratio(tss_bed, a_peak_bed, b_peak_bed, useflank=None):

    if useflank > 0 :
        a_peak_bed = a_peak_bed.flank(b=useflank, g=g).remove_invalid().sort()
        b_peak_bed = b_peak_bed.flank(b=useflank, g=g).remove_invalid().sort()

    if len(a_peak_bed) <= b_peak_bed:
        smaller = a_peak_bed
        larger = b_peak_bed
    else:
        smaller = b_peak_bed
        larger = a_peak_bed

    s_tss_i = tss_bed.intersect(smaller).sort().merge()
    s_l_tss_i =  s_tss_i.intersect(larger).sort().merge()

    output = float(len(s_l_tss_i)) / float(len(s_tss_i))

    if output > 1.0:
        # print "a b intersect: ", len(ab_i)
        # print "a: ", len(a_peak_bed)
        # print "b: ", len(b_peak_bed)
        # print "a b t intersect: ", len(tss_bed.intersect(a_peak_bed).intersect(b_peak_bed))
        # print "a t intersect: ", len(a_peak_bed.intersect(tss_bed))
        # print "b t intersect: ", len(b_peak_bed.intersect(tss_bed))
        # print "a.fn: ", a_peak_bed.fn
        # print "b.fn: ", b_peak_bed.fn
        print "smaller: ", len(smaller)
        print "larger: ", len(larger)
        print "smaller tss: ", len(s_tss_i)
        print "smaller tss larger: ", len(s_l_tss_i)
        raw_input()

    return output


def calculate_TSS_co_binding_minimum_overlap_ratio_multiprocessing(tss_bed_filename, bed_file_tuple_list, useflank=None):
    #this will first create a list of tuples for the proper call to the calculate function
    arg_tuple_list = []
    tss_bed = pbt.BedTool(tss_bed_filename)
    for t in bed_file_tuple_list:
        arg_tuple_list.append((tss_bed, t[0], t[1], useflank))

    pool = mp.Pool(mp.cpu_count())
    output = pool.map(_untuple_and_run_co_binding_minimum_overlap, arg_tuple_list)

    return output

def _untuple_and_run_co_binding_minimum_overlap(tup):
    try:
        # print tup
        bt1 = pbt.BedTool(tup[1])
        bt2 = pbt.BedTool(tup[2])
        uf = tup[3]
        # print bt1.fn
        # print bt2.fn
        return calculate_TSS_co_binding_minimum_overlap_ratio(tup[0], bt1, bt2, useflank=uf)
    except TypeError, e:
        print "This was bad."
        print tup


def convert_tuple_list_to_dataframe(tuplelist, row_t_idx = 0, col_t_idx = 1, val_t_idx =2):
    rows = []
    cols = []
    for _,t in enumerate(tuplelist):
        if t[row_t_idx] not in rows:
            rows.append(t[row_t_idx])

        if t[col_t_idx] not in cols:
            cols.append(t[col_t_idx])

    df = pd.DataFrame(index=rows, columns=cols)
    for _, t in enumerate(tuplelist):
        df.set_value(t[row_t_idx], t[col_t_idx], t[val_t_idx])

    return df

motif_database_dir = os.path.join(dir_path, "motif_databases")

def fimo_search_for_motif(fasta_file, motif="MA0486.2", motif_database=os.path.join(motif_database_dir, "JASPAR", "JASPAR_CORE_2016_vertebrates.meme"), thresh=1e-4, verbosity=0):
    if not os.path.isfile(fasta_file):
        raise IOError("{f} is not a file.".format(f=fasta_file))
    if not os.path.isfile(motif_database):
        raise IOError("{f} is not a file.".format(f=motif_database))


    command_text = """fimo -motif {m} -verbosity {v} -thresh {th} -text {db} {f}""".format(m=motif, db=motif_database, f=fasta_file, th=thresh, v=verbosity)
    outtext = subprocess.Popen(command_text, shell=True, stdout=subprocess.PIPE)
    text = outtext.stdout.read()
    if not text:
        return None
    else:
        return pd.read_table(StringIO.StringIO(text))

def calculate_total_basepair_length(bed):
    length_count = 0
    for feature in bed:
        length_count += (int(feature.fields[2]) - int(feature[1]))
    return  float(length_count)

def compute_fisher_exact_TF_cobinding_from_bedfiles(tf1_bedfile, tf2_bedfile, tss_bedfile, tss_length=None, tss_distance_thresh=2000, tf_distance_thresh=2000):
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    argd = {k:values[k] for k in ["tf1_bedfile", "tf2_bedfile", "tss_bedfile"]}

    for k in argd:
        if isinstance(argd[k], pbt.BedTool):
            continue
        elif os.path.isfile(argd[k]):
            argd[k] = pbt.BedTool(argd[k])
        else:
            raise TypeError("{k} is not a bedfile or a BedTool".format(k=k))


    tf1 = argd["tf1_bedfile"]
    tf2 = argd["tf2_bedfile"]
    tss = argd["tss_bedfile"]

    if not isinstance(tss_length, numbers.Number):
        uniques = []
        for feature in tss:
            posstr = ":".join(feature.fields[0:3])
            if posstr not in uniques:
                uniques.append(posstr)

        tss_length = len(uniques)

    #assume that the peak will be the last column of the bedfile
    #otherwise we still need to figure out which column of the bedfile
    #contains the start of the TF1 region

    tss_tf1_close = tf1.sort().closest(tss.sort(), D="a")
    # print len(tss_tf1_close)
    # print len(tf1)
    # raw_input()
    ex_fields = tss_tf1_close[0].fields
    chr_idxs = [i for i in range(len(ex_fields)) if re.match(chr_start_pattern, ex_fields[i])]
    tss_start_field_idx = chr_idxs[1] + 1
    tf1_peak_field_idx = tss_start_field_idx - 2
    tf1_fields = []
    uniques = []
    for feature in tss_tf1_close:
        fields = feature.fields
        if int(fields[-1]) == -1: continue
        tf1_peak_pos = int(fields[1]) + int(fields[tf1_peak_field_idx])
        dist = min( abs( tf1_peak_pos - int(fields[tss_start_field_idx]) ), abs( tf1_peak_pos - int(fields[tss_start_field_idx+1])) )
        btdist = abs(int(fields[-1]))
        if dist <= tss_distance_thresh:
            p = ":".join(fields[tss_start_field_idx:tss_start_field_idx+3])
            if p not in uniques:
                # print btdist, dist
                tf1_fields.append(fields)
                uniques.append(p)


    #repeat for tf2
    tss_tf2_close = tf2.sort().closest(tss.sort(), D="a")
    ex_fields = tss_tf2_close[0].fields
    chr_idxs = [i for i in range(len(ex_fields)) if re.match(chr_start_pattern, ex_fields[i])]
    tss_start_field_idx = chr_idxs[1] + 1
    tf2_peak_field_idx = tss_start_field_idx - 2
    tf2_fields = []
    uniques = []
    for feature in tss_tf2_close:
        fields = feature.fields
        if int(fields[-1]) == -1: continue
        tf2_peak_pos = int(fields[tf2_peak_field_idx]) + int(fields[1])
        dist = min( abs( tf2_peak_pos - int(fields[tss_start_field_idx]) ), abs( tf2_peak_pos - int(fields[tss_start_field_idx+1]) ) )
        if dist <= tss_distance_thresh:
            p = ":".join(fields[tss_start_field_idx:tss_start_field_idx+3])
            if p not in uniques:
                tf2_fields.append(fields)
                uniques.append(p)

    #see how many times the peaks for tf1 and tf2 are within the distance thresh

    tf1_tf2_close = tf1.sort().closest(tf2.sort(), D="a")
    tf1_tf2_close_fields = []
    ex_fields = tf1_tf2_close[0].fields
    chr_idxs = [i for i in range(len(ex_fields)) if re.match(chr_start_pattern, ex_fields[i])]
    tf2_start_idx = chr_idxs[1] + 1
    tf1_peak_idx = chr_idxs[1] - 1
    uniques = []
    for feature in tf1_tf2_close:
        fields = feature.fields
        tf1_peak_pos = int(fields[1]) + int(fields[tf1_peak_idx])
        tf2_peak_pos = int(fields[tf2_start_idx]) + int(fields[-2])
        dist = abs( tf1_peak_pos - tf2_peak_pos )
        if dist <= tf_distance_thresh:
            p = ":".join(fields[0:3])
            if p not in uniques:
                tf1_tf2_close_fields.append(fields)
                uniques.append(p)

    overlap = len(tf1_tf2_close_fields)


    regions_tf1 = len(tf1_fields)
    regions_tf2 = len(tf2_fields)

    # print regions_tf1
    # print regions_tf2
    # print tss_length
    # print overlap

    return compute_fisher_exact_TF_cobinding(overlap, regions_tf1, regions_tf2, tss_length)



def compute_fisher_exact_TF_cobinding(n_regions_overlap,
                                      n_regions_tf1,
                                      n_regions_tf2,
                                      n_total_possible_regions=None,
                                      null_hypothesis="two-sided"):

    n_regions_with_tf1_tf2 = n_regions_overlap
    n_regions_with_tf1_wo_tf2 = n_regions_tf1 - n_regions_overlap
    n_regions_wo_tf1_with_tf2 = n_regions_tf2 - n_regions_overlap
    n_regions_wo_tf1_wo_tf2 = n_total_possible_regions - n_regions_tf1 - n_regions_tf2

    contingency = [[n_regions_with_tf1_tf2, n_regions_wo_tf1_with_tf2],
                   [n_regions_with_tf1_wo_tf2, n_regions_wo_tf1_wo_tf2]]

    # print contingency

    oddsr, pv = st.fisher_exact(np.array(contingency), null_hypothesis)
    return oddsr, pv


def compute_chisquare_with_yates_cobinding(n_regions_overlap,
                                           n_regions_tf1,
                                           n_regions_tf2,
                                           n_total_possible_regions=None):

    n_regions_with_tf1_tf2 = n_regions_overlap
    n_regions_with_tf1_wo_tf2 = n_regions_tf1 - n_regions_overlap
    n_regions_wo_tf1_with_tf2 = n_regions_tf2 - n_regions_overlap
    n_regions_wo_tf1_wo_tf2 = n_total_possible_regions - n_regions_tf1 - n_regions_tf2

    contingency = [[n_regions_with_tf1_tf2, n_regions_wo_tf1_with_tf2],
                   [n_regions_with_tf1_wo_tf2, n_regions_wo_tf1_wo_tf2]]

    chi2, p, dof, ex = st.chi2_contingency(contingency)

    return chi2, p, dof, ex


def rolling_sum(array, n=10):
    if isinstance(array, list):
        array = np.array(array)
    elif not isinstance(array, np.ndarray):
        raise  TypeError("'array' must be a list or numpy array.")

    out = np.zeros((array.shape[0],))
    for i in range(0, len(out)):
        s = np.nansum(array[i:i+n])
        out[i] = s

    return out
