import os
import pybedtools as pbt

data_dir = os.path.join("..", "data_files")
hsf1_basal_bed_file = os.path.join(data_dir, "hsf1_0.17_merged_regions_hg19.refiltered.bed")
hsf1_HS_bed_file = os.path.join(data_dir, "hsf1_HS_merged_regions_hg19.bed")

basal = pbt.BedTool(hsf1_basal_bed_file)
HS = pbt.BedTool(hsf1_HS_bed_file)


basal_exclusive_outfile = os.path.join(data_dir, "hsf1_nonHS_exclusive.bed")
HS_exclusive_outfile = os.path.join(data_dir, "hsf1_HS_exclusive.bed")

#do the same for the HS file

iHS = HS.sort().intersect(basal.sort())
iBasal = basal.sort().intersect(HS.sort())
# print iHS[0].fields
# print len(iHS)
# print len(intersection)
# raw_input()
good_ihs_fields = []
uniques = []
for feature in iHS:
    if feature.name not in uniques:
        good_ihs_fields.append(feature.fields)
        uniques.append(feature.name)

good_ibasal_fields = []
uniques = []
for feature in iBasal:
    if feature.name not in uniques:
        good_ibasal_fields.append(feature.fields)
        uniques.append(feature.name)

print len(good_ihs_fields)
print len(good_ibasal_fields)

ibasal_bed = pbt.BedTool(good_ibasal_fields)
ihs_bed = pbt.BedTool(good_ihs_fields)

basal_exclusive_fields = []
basal_intersection_names = [feature.name for feature in ibasal_bed]

for feature in basal:
    if feature.name not in basal_intersection_names:
        basal_exclusive_fields.append(feature.fields)

basal_exclusive_bed =pbt.BedTool(basal_exclusive_fields).sort()
basal_exclusive_bed.saveas(basal_exclusive_outfile)

hs_exclusive_fields = []
hs_intersection_names = [feature.name for feature in ihs_bed]

for feature in HS:
    if feature.name not in hs_intersection_names:
        hs_exclusive_fields.append(feature.fields)

hs_exclusive_bed = pbt.BedTool(hs_exclusive_fields).sort()
hs_exclusive_bed.saveas(HS_exclusive_outfile)

#merge the intersection files together

i_merged = ihs_bed.cat(*[ihs_bed, ibasal_bed], postmerge=True)
i_merged.sort()

merged_fields = []
for i, feature in enumerate(i_merged):
    name = "intersection-feature{n}|{c}-{s}-{e}-1".format(n=i, c=feature.fields[0], s=feature.fields[1], e=feature.fields[2])
    merged_fields.append(feature.fields + [name])

i_merged = pbt.BedTool(merged_fields).sort()

merged_outfile = os.path.join(data_dir, "hsf1_HS_nonHS_intersection_hg19.bed")
i_merged.saveas(merged_outfile)
