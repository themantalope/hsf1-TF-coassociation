import pybedtools


a = pybedtools.BedTool("data_files/macs2.idrOptimal.bf.dm3.10T-DFD-GFP.3229_dfd_Embryos-16-24-hr_UofC_stn_Rep0_VS_3229-Input-R1_q30_Rep0.narrowPeak")
b = pybedtools.BedTool("data_files/macs2.idrOptimal.bf.dm3.Iso1.CNC-AdultFemale_CNC_AdultFemale_UofC_stn_Rep0_VS_2011-1017_q30_Rep0.narrowPeak")
print a.intersect(b)