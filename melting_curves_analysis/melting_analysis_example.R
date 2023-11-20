
path = '/Users/tomparee/Desktop/melting_curves_analysis/'

#################################################################################
####################### ANALYSE THE MELTING CURVES AND EXTRACT METRICS ##########
#################################################################################

# Import function; all package on the function_melting_analysis.R must be installed
source(paste0(path, "function_melting_analysis.R"))

# Import design (samples details and their corresponding qPCR well)
design = read_xlsx(paste0(path,"design_example/DESIGN_example.xlsx"))


# Import melting curve data (melting_example.txt; fluorescence ~ temperature, for each well) and format them.
# To obtain the melting_example.txt file in the Roche480 software, you go into melting curve analysis and right-click on the melting curve figures => export data => .txt
rawd = readHRM( paste0(path, "rawdata_example/melting_example.txt"))

# Optional:
# Import and format the aplification phase data (ampli_example.txt)
# This file is only used to filter out samples with falied/weird amplification
# To obtain "ampli_example.txt" in the Roche480 software, you right-click on the amplification phase figure => export data => .txt
ampli = readAmpli( paste0(path, "rawdata_example/ampli_example.txt"))

# Not necessary if only one qPCR plate
rawd$plate_qpcr_id = "plate1"
ampli$plate_qpcr_id = "plate1"

# to bind several plates (If so, must have a plate_qpcr_id):
#rawd=rbind(rawd1, rawd2)
#ampli = rbind(ampli1, ampli2)

# Must have a unique identifier for each well (well x plate)
rawd$id = paste0(rawd$plate_qpcr_id,"#",rawd$well_pcr)
ampli$id = paste0(ampli$plate_qpcr_id,"#",ampli$well_pcr)


#### HRM is the main function
# It normalize the raw melting curve
# Search for the different melting phase
# Extract few metrics from the curve that can be then used to infer the rec-1 mutant frequency and to filter abnormal curves
# (!) A lot of warnings are expected (empry well or failed amplification usually). Don't worry.
normd = HRM(rawd, ampli)

#### normd is a list of three elements:

#### normd[[1]] contain the melting curves (raw and normalized).
#temp   well_pcr     fluo      plate_qpcr_id        id       normfluo    normfluoEBS      deriv
#59.51    A1       10.75183        plate1        plate1#A1      100           100      0.08154002
#fluo = raw fluorescence;
#normfluo = normalized fluorescence (basic); normfluoEBS = normalized fluo by EBS (better in theory) 
#=> see "High resolution melting curve analysis with MATLAB-based program", Li, 2016 for details of these methods.
# deriv in the derivative of the raw fluo (i.e., melting peaks)

#### normd[[2]] contains the extracted metrics from the curves
#id         Tbf   Taf   raw.decrease  norm.decrease   normEBS.decrease   het.decrease   plateau_slope      maxfluo
#plate1#A1 73.40 83.52    0.2921659     0.1983239        0.2587710       7.389785      0.0000000000    10.629590
#Taf & Tbf: are the temperatures before and after melting detected by the algorithm
# xxx.decrease: The proportion of a fluorescence decrease during melting attributed to the rec-1 mutant by the algorithm.
# xxx.decrease can be obtained from the different curves (raw or normalized). Usually I used "normEBS.decrease" to infer rec-1 mutant frequency
# het.decrease is the decrease of fluorescence during the heteroduplex melting. This is only used to detect outlier curves.
# plateau_slope is the slope at the end of amplification (did amplification reach a plateau). This is only used to detect outlier curves.
# maxfluo is the fluo after amplification / before melting.  This is only used to detect outlier curves.

#### normd[[3]] contains temperatures delimiting each melting phase (heteroduplux, wt, mut) for each curve.
# start  inflex.bf  peak    inflex.af   end   phase.genotype phase.genotype.score raw.decrease norm.decrease   normEBS.decrease         id
# 73.40     75.59   76.30     76.79     77.08      het                    1        1.3088392      5.818178        10.75188        plate1#A1
# five first column are different landmark of each melting peaks (not very important for you)
# phase.genotype: the identified melting phase (het, wt, mut)


# What you care about is normd[[2]]:
meltinfo = normd[[2]]
meltinfo= merge(meltinfo, design) # Add design information

#badscore function attribute a score indicating outlier curves (0 to 3; more is worst).
# +1 is being outside the 99% interval for the plateau_slope, max.fluo or het_decrease
meltinfo=badscore(meltinfo)
# To plot the metric distribution ~ of normEBS.decrease, use plot.bad
# Outlier in red; other samples in blue; calibration samples in black.
plot.bad(meltinfo) 

# You should get rid of the curve with a high badscore
badcurve = meltinfo$id[meltinfo$badscore>0 & is.na(meltinfo$calib)]
normd=subset.ID(normd, badcurve,out=T)
meltinfo = normd[[2]]
meltinfo= merge(meltinfo, design)



#You may want to look at the melting phase decomposition, especially for outliers
# It is possible that the algorithm did not decompose properly the melting curve
# If it is the case, you should also filter out these samples.
# Example: to look at the curve with id = "plate1#A19", you do:

#plot.melting.phase(normd = normd, id = "plate1#A19")
plot.melting.phase2(normd = normd, id = "plate1#A19") #Purple = heteroduplex; red = mutant; blue = wild-type.


save(meltinfo, file = paste0(path, "melting_info_example.Rdata"))


#################################################################################
################### INFER REC-1 MUTANT FREQUENCY FROM DECREASE.NORMEBS ##########
#################################################################################

library(gam) # Generalized additive model 
load(paste0(path, "melting_info_example.Rdata"))


data = subset(meltinfo, is.na(calib)) #  samples
calib = subset(meltinfo, !is.na(calib)) # calibration 
# Train the model on the calibration samples
gam_model <-gam(ratio~s(normEBS.decrease, df=4)  ,data = calib)
summary(gam_model)

# If you want to plot the model
pgam =data.frame(xx=seq(0.0,1,0.01),
                 yy=predict(gam_model, newdata = data.frame(normEBS.decrease=seq(0,1,0.01))) )


# Black line = model
# Dot are calibration samples; colors are independent calibration curves
ggplot()+
  geom_point(data=subset(meltinfo, !is.na(calib)),
             aes(normEBS.decrease, ratio, color=calib), shape=1, size=3,stroke=1)+
  theme_Publication2()+
  xlab(expression(d[mut]))+ylab("rec-1 mutant frequency")+
  theme(axis.title = element_text(size=17,colour = "black", face = "bold"),    
        axis.text = element_text(size=15, colour = "black"),
        strip.text = element_text(size = 15, color = "black", face = "bold"),
        plot.title = element_text(size = 17, hjust = -0.1, face = "bold"),
        legend.position = "none")+
  geom_line(data=pgam, aes(x=xx, y=yy),size=1)



data$predicted_freq = predict(gam_model, newdata=data)





