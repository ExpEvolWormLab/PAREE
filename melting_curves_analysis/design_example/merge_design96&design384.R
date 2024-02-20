# We need a design data.frame containing details about samples and their position in the 384 qPCR plate

#well_pcr   ratio   calib   sample   pop  generation  plate_qpcr_id   id
#O12         NA     <NA>   G14_70    70         14        plate1      plate1#O12
#N7.         0.10   calib1  calib1   NA         NA        plate1      plate1#N7

# => Must have an "id" column with a unique identifier of each well (well x plate)
# => If several qPCR plate, must have a "plate_qpcr_id" column to identifiy the different qPCR plate
# => Should have "sample" column with sample name
# => Better if there is a "ratio" column and a "calib" indicating calibration samples info (rec-1 mutant freq and calib info, respectively)

# If you transfer the samples in a 96 well plate and then into the 384 qPCR well plate with a multichanel
# you might want to do a design for the 96 well plate containing the samples info
# and a design for the 384 qPCR well plate where you only indicate the corresponding well of the 96-well plate.
# If so, you can use the code below to fuse these two separate design in teh format described above:


### IMPORT FUNCTION & PACKAGES
### (first, be sure than package in "function_melting_analysis.R" are installed)
path = '/Users/tomparee/Desktop/melting_curves_analysis/'
source(paste0(path, "function_melting_analysis.R"))


## import files

design384 <- read_excel(paste0(path,"design_example/example_design_plate384.xlsx"))
# A 384 matrix 16 * (24+1) containing the id of the well in plate96 (A1; H12; etc.) or a calibration details (calib1_0.10)
# First column id the row.names (on a 384 well plate)
# colnames = column names (on a 384 well plate)


plate96 = read_excel(paste0(path,"design_example/example_design_plate96.xlsx"))
# data.frame with sample details
# Must have column "well96" (A1; H12), "plate96" (96 well plate ID); "sample" (sample name)
# And any other column you wish.


# format.design fuse the design384 and plate96 details in a single data.frame:
# Goal associate a 
design = format.design(design384, plate96, row.names = T, calibId384 = "calib", keep.info.well96 = F)
#head:
#well_pcr ratio  calib    sample   pop generation
#O12       NA   <NA>      G14_70    70         14
#M7        0.10  calib1   calib1.   NA        NA
#=> well_pcr = well id in the 384 qPCR plate
#=> ratio = rec-1 mutant proportion; NA if a not a calibration samples
#=> calib = calib ID (usually: calib1; calib2; calib3)
#=> other columns are optional


# If several plate you may want to do rbind(design1, design2) with "plate_qpcr_id" column.
design$plate_qpcr_id = "plate1"  #Not necessary if only one 384 qPCR plate 

#Important: id is a unique identifier of each qPCR well
# well_pcr id x plate_qpcr id (just well_pcr is fine if only one qPCR plate)
# an "id" column must be present in the design for subsequent functions
design$id = paste0(design$plate_qpcr_id,"#", design$well_pcr)


write_xlsx(design, path = paste0(path,"design_example/DESIGN_example.xlsx"))



grid.arrange(pcurve, pp, heights =c(1,2))



