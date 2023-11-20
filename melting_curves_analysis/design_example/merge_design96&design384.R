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

# path = './Data/Direct_Fitness_cost/Competition_nacl/march2023/'
# 
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






rawd1 = readHRM( paste0(path, "raw_data/melting_Rpoly2_RP2#LSC1.txt"))
rawd2 = readHRM( paste0(path,  "raw_data/melting_Rpoly2_RP2#LSC2.txt"))
ampli1 = readAmpli( paste0(path, "raw_data/ampli_Rpoly2_RP2#LSC1.txt"))
ampli2 = readAmpli( paste0(path,  "raw_data/ampli_Rpoly2_RP2#LSC2.txt"))

rawd1$plate_qpcr = 1
rawd2$plate_qpcr = 2
ampli1$plate_qpcr = 1
ampli2$plate_qpcr = 2

rawd1$plate_qpcr_id = "RP2#LSC1"
rawd2$plate_qpcr_id = "RP2#LSC2"
ampli1$plate_qpcr_id = "RP2#LSC1"
ampli2$plate_qpcr_id = "RP2#LSC2"

rawd=rbind(rawd1, rawd2)
ampli = rbind(ampli1, ampli2)

rawd$id = paste0(rawd$plate_qpcr_id,rawd$well_pcr)
ampli$id = paste0(ampli$plate_qpcr_id,ampli$well_pcr)

normd = HRM(rawd, ampli)

meltinfo = normd[[2]]
meltinfo= merge(meltinfo, design)

meltinfo=badscore(meltinfo)

badcurve = meltinfo$id[meltinfo$badscore>0 & is.na(meltinfo$calib)]
plot.bad(meltinfo)

ggplot(subset(meltinfo, badscore == 0),
       aes(generation, normEBS.decrease, group=replicate))+
  stat_summary(geom="line")+geom_point()+
  facet_wrap(~replicate)


unique(meltinfo$replicate)
ggplot(subset(meltinfo, badscore == 0 & replicate %in% c(49, 51, 56, 61, 69, 73, 76, 79, 85, 93)),
       aes(generation, normEBS.decrease, group=replicate))+
  stat_summary(geom="line")+stat_summary(size=0.1)+ylim(0.125, 0.5)


ggplot(subset(meltinfo, badscore == 0 & !(replicate %in% c(49, 51, 56, 61, 69, 73, 76, 79, 85, 93))),
       aes(generation, normEBS.decrease, group=replicate))+
  stat_summary(geom="line")+stat_summary(size=0.1)+ylim(0.125, 0.5)

ggplot(subset(meltinfo, badscore == 0 & (replicate %in% c(49, 51, 56, 61, 69, 73, 76, 79, 85, 93)) &
                replicate %in% 53:72),
       aes(generation, normEBS.decrease, group=replicate))+
  stat_summary(geom="line")+stat_summary(size=0.1)+ylim(0.125, 0.5)+
  geom_smooth(method = "lm")

ggplot(subset(meltinfo, badscore == 0 & !(replicate %in% c(49, 51, 56, 61, 69, 73, 76, 79, 85, 93)) &
                replicate %in% 73:76),
       aes(generation, normEBS.decrease, group=replicate))+
  stat_summary(geom="line")+stat_summary(size=0.1)+ylim(0.125, 0.5)


ggplot(subset(meltinfo, badscore == 0 & replicate %in% c(73, 76, 79, 85, 93)),
       aes(generation, normEBS.decrease, group=replicate))+
  stat_summary(geom="line")+stat_summary(size=0.1)+ylim(0.125, 0.5)

ggplot(subset(meltinfo, badscore == 0 & replicate %in% c( 56, 61, 69)),
       aes(generation, normEBS.decrease, group=replicate))+
  stat_summary(geom="line")+stat_summary(size=0.1)+ylim(0.125, 0.5)


unique(meltinfo$replicate)
ggplot(subset(meltinfo, badscore == 0 & replicate %in% c(49, 51, 56, 61, 69, 73, 76, 79, 85, 93)),
       aes(generation, normEBS.decrease, group=replicate))+
  stat_summary(geom="line")+stat_summary()+ylim(0.125, 0.5)

ggplot(meltinfo, aes(ratio, normEBS.decrease))+geom_point()


#normd=subset.ID(normd, c(badcurve, "DS#TP1E19", "DS#TP2P15", "DS#TP2P12"),out=T)

meltinfo = normd[[2]]
meltinfo= merge(meltinfo, design)


save(meltinfo, file = paste0(path, "melting_info_Rpoly2_LSC.Rdata"))





library(gam)


gamRP2 <-gam(ratio~s(normEBS.decrease, df=4)  ,data = subset(meltinfo, !is.na(ratio) & plate_qpcr_id=='RP2#LSC2'))
summary(gamRP2)

pgam =data.frame(xx=seq(0.0,1,0.01),
                 yy=predict(gamRP2, newdata = data.frame(normEBS.decrease=seq(0,1,0.01))) )



ggplot()+
  geom_point(data=subset(meltinfo, !is.na(ratio) & plate_qpcr_id=='RP2#LSC2'),
             aes(normEBS.decrease, ratio, color=calib), shape=1, size=3,stroke=1)+
  theme_Publication2()+
  xlab(expression(d[mut]))+ylab("rec-1 mutant frequency")+
  theme(axis.title = element_text(size=17,colour = "black", face = "bold"),    
        axis.text = element_text(size=15, colour = "black"),
        strip.text = element_text(size = 15, color = "black", face = "bold"),
        plot.title = element_text(size = 17, hjust = -0.1, face = "bold"),
        legend.position = "none")+
  geom_line(data=pgam, aes(x=xx, y=yy),size=1)



pp = plot.melting.phase2(normd, ids=c("RP2#LSC2E9","RP2#LSC2L8", "RP2#LSC2A7"), titles = paste0(c(0,50,100), "%"))

fluo = merge(normd[[1]], design)

source("./basics_TP.R")

pr = ggplot(subset(fluo, calib == 'calib2' & temp >75 & temp<85), aes(temp, fluo, color=ratio, group=id))+
  geom_line()+scale_color_continuous(low='blue', high='red')+theme_Publication2()+
  ylab("fluorescence")+xlab("Temperature")+
  theme(legend.position = "none")


pn = ggplot(subset(fluo, calib == 'calib2' & temp >75 & temp<85), aes(temp, normfluoEBS, color=ratio, group=id))+
  geom_line()+scale_color_continuous(low='blue', high='red')+theme_Publication2()+
  ylab("Normalized fluo. (EBS)")+xlab("Temperature")+
  theme(legend.position = "none")


pd = ggplot(subset(fluo, calib == 'calib2' & temp >75 & temp<85), aes(temp, deriv, color=ratio, group=id))+
  geom_line()+scale_color_continuous(low='blue', high='red', name="mutant \n proportion \n")+theme_Publication2()+
  ylab("Derivartive")+xlab("Temperature")

pcurve = grid.arrange(pr, pn, pd, nrow=1, widths = c(1,1,1.1))


grid.arrange(pcurve, pp, heights =c(1,2))



