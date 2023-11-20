#DESIGN_example.xlsx 
samples info with their corresponding well on the qPCR plate. 

well_pcr   ratio   calib   sample   pop  generation  plate_qpcr_id   id
O12         NA     <NA>   G14_70    70         14        plate1      plate1#O12
N7.         0.10   calib1  calib1   NA         NA        plate1      plate1#N7
 => Must have an "id" column with a unique identifier of each well (well x plate)
 => If several qPCR plates, must have a "plate_qpcr_id" column to identifiy the different qPCR plate
 => Should have "sample" column with sample name
 => Better if there is a "ratio" column and a "calib" indicating calibration samples info (rec-1 mutant freq and calib info, respectively)

With many samples, I usually transfer them in a 96-well plate before transferring them to the 384 qPCR plate. So you may want to do a 96 well plate design and a 384 well plate design and merge the two. This can be done with the following script and files.


# example_design_plate96.xlsx
The sample information and their corresponding well on the 96-well plate
=> Must have column "well96" (A1; H12), "plate96" (96 well plate ID); "sample" (sample name)
=> And any other column you wish (here, generation and population)

# example_design_plate384.xlsx
Indicate the 96-well within the 384 qPCR plate

#merge_design96&design384.R
Fuse ‘example_design_plate96.xlsx’ and ‘example_design_plate384.xlsx’ into DESIGN_example.xlsx 

