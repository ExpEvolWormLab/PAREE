library(segmented)
library(readxl)
require(plyr)
require(data.table)
library(reshape)
library(ggplot2)
library(ggResidpanel)
library("birk")
library(fitdistrplus)
library(lme4)
library(lmerTest)
library("optimx")
library(emmeans)

format.design = function(design384, plate96, row.names = T, calibId384 = "calib", keep.info.well96 = F){
  #design384 is the design of the 384 plate (see github for example)
  # => row = row of the qPCR plate: column = column of the qPCR plare
  # => If first column is the row names => option: row.names = T
  # value of each cell for a sample = "well96_plate96" (example: "H3_plateA"; "A1_plate2"; ...)
  # value of each cell for a calib = "calibId384_freq" (example: "calib1_0.10"; "calib2_0"; "calibTest_0.5"; ...)
  # (except first row if row.names = T)
  # Note: if calibration info are contained in plate96 file, calib cells can have the same format as a sample (ex: "A1_plateCalib")
  # => if so, calibId384 = NULL
  # Note= calibId384 should not contain "_"
  
  
  #plate96 is the design of the 96 well
  # column contain at least "well96", "sample" and "plate96" (!) with these name
  # other columns may contain as many info as wanted ("generation", "strain", "treatment", etc.)
  # "well96" should be a well id corresponding to the cell in the design384 (example: "H3", "A1", ...)
  # "plate96" should be a plate id corresponding to the cell in the design384 (example: "plateA", "plate2", ...)
  
  
  
  if(!is.null(calibId384)){ if(grepl("_", calibId384)){ stop("calibId384 should not contain a underscore") }}
  if(sum(c("well96", "sample", "plate96") %in% colnames(plate96))!=3){stop("plate96 column names should include 'well96', 'sample','plate96'")}
  
  design384 = as.matrix(design384)
  if(row.names){
    # if the first column is the row names
    rownames(design384)= design384[,1]
    design384 = design384[,-1]
  }
  
  # melt 
  design384 = reshape2::melt(design384)
  
  rowisnum = sum(design384$Var1 %in% 1:24) == length(design384$Var1)
  colisnum =  sum(design384$Var2 %in% 1:24) == length(design384$Var2)
  
  if(!(rowisnum | colisnum) | (rowisnum & colisnum)) stop("in .xlsx file, row should be letters from A to P and column be numbers from 1 to 24,  or inversely")
  
  # Create a column containing 384 qPCR well identifier ("A1", "P12", etc.)
  if(rowisnum) design384$well_pcr = paste0(design384$Var1, design384$Var2)
  if(colisnum) design384$well_pcr = paste0(design384$Var1, design384$Var2)
  
  
  # Create a column "ratio" containing the rec-1 mutant frequency if this info is available in the design384 for the calib
  # For the sample, the value is NA
  if(!is.null(calibId384)){
    
    design384$ratio = data.table::tstrsplit(design384$value, '_')[[2]]
    design384$ratio[!grepl(calibId384,design384$value)]=NA
    
    # in case number are in format "0,9" instead of "0.9"
    if(sum(grepl(",",design384$ratio[grepl(calibId384,design384$value)]))>0){
      
      design384$ratio[grepl(calibId384,design384$value)] = sapply(design384$ratio[grepl(calibId384,design384$value)], function(x){
        if(grepl(",",x)){
          x = data.table::tstrsplit(x,",", fixed=T)
          x = paste0(x[[1]], '.', x[[2]])
        }
        x
      })
      
    }
    
    
    design384$ratio = as.numeric(design384$ratio)
    
    design384$plate96[is.na(design384$ratio)]=data.table::tstrsplit(design384$value[is.na(design384$ratio)], '_')[[2]]
    design384$calib[!is.na(design384$ratio)]=data.table::tstrsplit(design384$value[!is.na(design384$ratio)], '_')[[1]]
    design384$well96=NA
    design384$well96[is.na(design384$ratio)]=data.table::tstrsplit(design384$value[is.na(design384$ratio)], '_')[[1]]
    
    
    design384 = design384[,c('well96','well_pcr', 'ratio', 'plate96', 'calib')]
    
    calib = subset(design384, is.na(well96))
    
  }else{
    
    design384$well96=data.table::tstrsplit(design384$value, '_')[[1]]
    design384$plate96=data.table::tstrsplit(design384$value, '_')[[2]]
    
    design384 = design384[,c('well96','well_pcr', 'plate96')]
    
    
  }
  
  
  
  
  #merge design384 with plate96 
  design384 = merge(subset(design384, !is.na(well96)), plate96,all=T)
  
  
  if(!is.null(calibId384)){
    #bind df conatining sample and df containing calib info
    misscol = colnames(design384)[is.na(match(colnames(design384), colnames(calib)))]
    calib = cbind(calib, setNames( as.data.frame(matrix(NA, ncol=length(misscol), nrow=nrow(calib))), misscol))
    calib$sample = calib$calib
    design384=rbind(design384,calib)
  }
  
  if(keep.info.well96==F){
    design384=design384[,-which(colnames(design384) %in% c("well96", "plate96"))]
  }
  
  
  design384
  
}



# format.design = function(design, plate96){
#   
#   design = as.matrix(design)
#   rownames(design)= design[,1]
#   design = design[,-1]
#   design = reshape2::melt(design)
#   
#   design$well_pcr = paste0(design$Var2, design$Var1)
#   
#   design$plate96=as.numeric(substr(design$value,1,1))
#   
#   design$well96=NA
#   design$well96=sapply(design$value, function(x){ as.character(substr(x, 2, length(data.table::tstrsplit(x, ''))))})
#   design = design[,c('well96','well_pcr', 'plate96')]
#   
#   if(!("well96" %in% colnames(plate96) & "plate96"  %in% colnames(plate96))){
#     stop("96-well plate design must contain column 'well96' and 'plate96'")
#   }
#   
#   design = merge(subset(design, !is.na(well96)), plate96,all=T)
#   
#   design
# }

#### FUNCTION ####
readHRM <- function(datf){
  rawd = read.delim2(datf)
  rawd = rawd[,-seq(3, ncol(rawd), by=2)]
  rawd = rawd[!duplicated(rawd[,1]),]
  rawd= melt(as.data.table(rawd),1)
  names(rawd) = c('temp', 'well_pcr', 'fluo')
  rawd$temp = as.numeric(rawd$temp)
  rawd$well_pcr = as.character(rawd$well_pcr)
  rawd$fluo = as.numeric(rawd$fluo)
  rawd$well_pcr <- tstrsplit(as.character(rawd$well_pcr), "..Sample")[[1]]
  #rawd = subset(rawd, between(temp, 76.5, 84))
  return(as.data.frame(rawd))
}


readAmpli = function(datf, ncycle=NULL){
  require(reshape)
  require(data.table)
  ampli <- read.delim2(datf)
  ampli = ampli[,-1]
  dt = diff(ampli[,1])
  if(is.null(ncycle)){
    cycles = which(dt > dt[1] - 10)
    cspan = unique(diff(cycles))
    if(!(length(cspan) == 1 & cspan[1]==1)){stop("Provide the number of PCR cycle in the ncycle argument")}
    ncycle = cycles[length(cycles)] + 1
    ampli = ampli[1:ncycle,]
    
    print(paste0("Found ", ncycle, " PCR cycles. Is it correct? If not, provide it in the ncycle argument"))
  }else{
    ampli = ampli[1:(ncycle), ]
  }
  
  ampli = ampli[,-seq(3,ncol(ampli),2)]
  colnames(ampli)[2:ncol(ampli)] = data.table::tstrsplit(colnames(ampli)[2:ncol(ampli)], "..", fixed=T)[[1]]
  ampli[,1] = 1:ncycle
  ampli= melt(as.data.table(ampli),1)
  names(ampli) = c('cycle', 'well_pcr', 'fluo')
  return(as.data.frame(ampli))
  
}


slopePlateau = function(ampli){
  plateau = do.call(rbind, lapply(split(ampli, ampli$id), function(amplix){
    #amplix = subset(ampli, id=='D14')
    amplix$relfluo = (amplix$fluo - min(amplix$fluo))/(max(amplix$fluo)-min(amplix$fluo))
    #plot(amplix$cycle, amplix$fluo)
    slope = log(amplix$relfluo[nrow(amplix)-0]) - log(amplix$relfluo[nrow(amplix)-2])
    if(slope<0) slope=0; if(max(amplix$fluo)<3) slope=NA
    data.frame(id = amplix$id[1], plateau_slope = slope, maxfluo = max(amplix$fluo))
  }))
  
  plateau
  
}

############################################

#data = data.frame(temp, fluo) for one well

#data = subset(rawd, id=='TP1A6')
#peak.delimitation(data)

peak.delimitation = function(data, trange.mut = c(77.5,79.75), trange.wt = c(79.75,82)){
  
  #####################################################################
  #### Data is the raw data (temp,fluo) for one single curve (ex: data=subset(rawd, id=="G12"))
  #### The algorithm use the first curve derivative (melting peaks), the second and the thirs
  #### To detect sign of peak or bump which correspond to the melting of different alleles
  #### 1) It first detect the main peak and delimit it using derivative first and second
  #### Note: the peak delimitation is done with the peak.limits function
  #### 2) If other peaks (outside the range of the first one), delimit them too
  #### 3) Because peaks are against each other and have different size, some peak can be just "bumps" (no local maxima), 
  #### =>find them thks to deriv second
  #### 4) Do a first filtration based one the sd(deriv2), area(deriv1) to remove false peak (noise)
  #### when % mutant is really low, its melting is really hard to distinguish
  #### => but possible to find on deriv 3 using notably a segmentation algorithm to confirm
  #### 5) Attribute the genotype of peaks based on provided melting range (trange.wt, trange.mut) and order
  #### If confident , score = 1, if maybe false score < 1
  
  
  ### Prepare the the melting data 
  # No duplicated temperature
  data = data[!duplicated(data$temp),]
  #Smooth the curve
  data$fluo = predict(loess(fluo~temp, data, span=0.1), newdata=data)
  # Calculate the inverse derivative
  datad = derivative(data, xname="temp", yname="fluo", smooth=0.1, return.df =T)
  
  # Find temp where the (inverse) derivative is at its max (the peak) and keep temp range +10 to -10
  tmax = datad$temp[which.max(datad$deriv)]
  temprange = tmax + c(-10,+10) 
  datad = datad[datad$temp>temprange[1] & datad$temp<temprange[2],]
  
  # Do the second and third derivative 
  datad$deriv2 = -derivative(datad, xname="temp", yname="deriv",smooth=0.1)
  datad$deriv3 = -derivative(datad, xname="temp", yname="deriv2",smooth=0.1)
  #plot(data$temp, data$fluo)
  
  ### Refine the melting range based on the background noise
  ## deriv second for background before and after melting
  bf.background = datad$deriv2[datad$temp<(tmax-9)]
  af.background = datad$deriv2[datad$temp>(tmax+9)]
  
  if(length(bf.background)==0|length(af.background)==0) stop("No melting range found. Likely due to not enough signal. Empty well?")
  
  #Background should be close to 0
  # If too high, it can be that there is a noise ""peak"" in the background
  # So try to trim the temperature range 0.5 by 0.5 degree until decent background noise
  if(abs(mean(bf.background)) > 0.01){
    
    high.bf.noise = T
    tx = c(tmax-10, tmax-9)+0.5
    bf2= datad$deriv2[datad$temp>tx[1] & datad$temp<tx[2]]
    bf2 = mean(abs(bf2))
    bf1 = mean(abs(bf.background))
    
    while(mean(abs(bf2)) > 0.01 | (bf2-bf1)>0.01){
      tx = tx+0.5
      bf1 = bf2
      bf2= datad$deriv2[datad$temp>tx[1] & datad$temp<tx[2]]
      bf2 = mean(abs(bf2))
      if(tx[2]>=tmax) stop("Background noise too high")
    }
    
    bf.background = datad$deriv2[datad$temp>tx[1] & datad$temp<tx[2]]
    
  }else{high.bf.noise = F}
  
  
  # which are over the noise threshold (10 x background noise)
  wx = which((abs(datad$deriv2) > 10*mean(abs(bf.background)) & datad$temp<tmax & datad$deriv2>mean(bf.background))|(abs(datad$deriv2) > 10*mean(abs(af.background)) & datad$temp>tmax))
  
  # Before and after melting temperatures
  Tbf = min(datad$temp[wx])- 2# Let a 2 degree marge before 
  Taf = max(datad$temp[wx])+0.5 # Let a smaller marge after (less noise after so sufficent)
  
  # If Tbf and Taf not found, stop here
  if(is.infinite(Tbf)|is.infinite(Taf)){stop("Not enough signal to find melting range")}
  
  # cut the melting data in the temperature range
  datad2 = datad[datad$temp>Tbf & datad$temp<Taf,]
  
  
  #plot(datad2$temp, datad2$deriv)
  
  ### Find strateging point on peaks (maxima and inflexion point)
  # Find the melting peak (local maxima of the second derivative)
  peaks = find.extremes(datad2$deriv)[[2]]
  
  # Find the inflexion points == local max and min in second derivative
  inflex.up = find.extremes(datad2$deriv2)[[2]]
  inflex.down = find.extremes(datad2$deriv2)[[1]]
  inflex.down=inflex.down[inflex.down>min(inflex.up)]
  inflex.down=inflex.down[1:length(inflex.up)]
  inflex=data.frame(up=inflex.up,down=inflex.down)
  inflex$temprange = datad2$temp[inflex$down]-datad2$temp[inflex$up]
  
  ### First control for noise
  # Deriv 3 min for noise sorting and later
  d3.min=find.extremes(datad2$deriv3)[[1]]
  
  noise = which(inflex$temprange<0.5) # If the range between a min and max is too low 
  if(length(noise)>0){
    
    # If to "noise" are against each other , merge their range
    tomerge = noise[diff(noise)==1]
    tomerge = tomerge[!is.na(inflex$down[tomerge+1])]
    
    while(length(tomerge)>0){
      inflex$down[tomerge[1]] = inflex$down[tomerge[1]+1]
      inflex = inflex[-(tomerge[1]+1),]
      inflex$temprange[tomerge] = datad2$temp[inflex$down[tomerge]]-datad2$temp[inflex$up[tomerge]]
      noise = which(inflex$temprange<0.5)
      tomerge = noise[diff(noise)==1]
      tomerge = tomerge[!is.na(inflex$down[tomerge+1])]
    }
    
    # Calculate the sd deriv 3 and the diff for each noise suspicion
    stat.noise = do.call(rbind, lapply(noise, function(i){
      x= inflex[i,]
      
      data.frame(diff.deriv2 = max(datad2$deriv2[x$down:x$up]) - min(datad2$deriv2[x$down:x$up]),
                 sd.deriv3 = sd(datad2$deriv3[x$down:x$up]))
     
    }))
    
    
    
    
    # If really a peak, the sd should be somewhat high 
    noise = noise[stat.noise$sd.deriv3<0.01 | stat.noise$diff.deriv2<0.01]
    if(length(noise)>0){
      inflex = inflex[-noise,]
    }
  }
  
  if(sum(inflex$temprange<0.1,na.rm=T)>0){inflex = inflex[-which(inflex$temprange<0.1), ]}
  
  
  
  #ggplot(datad2, aes(temp, deriv))+
  #  geom_point()+
  #  geom_vline(xintercept = datad2$temp[inflex$up], color='red')+
  #  geom_vline(xintercept = datad2$temp[inflex$down], color='blue',linetype='dashed')+
  #  geom_vline(xintercept = datad2$temp[peaks], color='green',linetype='dashed')
  
  
  # Only take the peaks that are after the first inflexion point (if before, small noise, not reel peaks)
  peaks = peaks[peaks>min(inflex.up) & peaks<max(inflex.down,na.rm =T) ]
  peaks = data.frame(peak = peaks, deriv = datad2$deriv[peaks])
  
  ### Begin by deltimiting the main peak (the highest)
  main.peak = peaks$peak[which.max(peaks$deriv)] # max derivative
  main.inflex = which(inflex$up<main.peak & inflex$down>main.peak) # The peak is between two inflexion point
  main.peak = peak.limits(inflex, datad2, which.row.inflex=main.inflex,which.is.peak=main.peak) # Find peak delimitation
  
  # Take the main peak info out of the "remaining" peaks
  peaks = peaks[peaks$peak<main.peak$start | peaks$peak>main.peak$end,]
  
  ##########
  #### If there are more peaks, delimit them
  other.peaks = NULL
  if(nrow(peaks)>0){
    
    # Find for each which inflexion point correspond
    peaks$whichinflex = unlist(lapply(1:nrow(peaks), function(i){
      x=peaks[i,]
      winflex = which(inflex$up<x$peak & inflex$down>x$peak)
      if(length(winflex)==0) winflex =NA
      winflex 
    }))
    
    peaks = peaks[!is.na(peaks$whichinflex),]
    
    if(nrow(peaks)>0){
      # If there are more than one peak between two inflexion point, the true one is the higest, trash the other
      dup = sapply(unique(peaks$whichinflex), function(i) sum(peaks$whichinflex==i)) > 1
      trash=NULL
      for(i in unique(peaks$whichinflex)[dup]){
        wdup = which(peaks$whichinflex == i)
        trash = c(trash, wdup[wdup!=wdup[which.max(peaks$deriv[wdup])]])
      }
      
      if(!is.null(trash)){peaks=peaks[-trash,]}
      
      # Now delimit for all the other peaks
      
      other.peaks = do.call(rbind, lapply(1:nrow(peaks), function(i){
        
        x=peaks[i,]
        peak.limits(inflex, datad2, which.row.inflex=x$whichinflex,which.is.peak=x$peak)
        
      }))
    }
    
  }else{
    
    other.peaks=NULL
  }
  
  
  
  peaks = rbind(main.peak, other.peaks)
  
  
  #### Then the remaining inflexion point correspont to "Bumps" 
  #### => peak but without maxima because it is against another one
  #### Find the bumps:
  inflex=inflex[!is.na(inflex$up) & !is.na(inflex$down),]
  inflex$p = 'bump'
  for(i in 1:nrow(peaks)){
    #i=1
    inflex[((inflex$up < peaks[i,]$end & inflex$up > peaks[i,]$start) 
            | (inflex$down < peaks[i,]$end & inflex$down > peaks[i,]$start)) ,]$p = "peak"
    
  }
  
  
  bumps = do.call(rbind, lapply(which(inflex$p == 'bump'), function(i){
    peak.limits(inflex, datad2, which.row.inflex=i,which.is.peak=NA)
    
  }))
  
  peaks = rbind(peaks, bumps)
  
  #ggplot(datad2, aes(temp, deriv))+
  #  geom_point()+
  #  geom_vline(xintercept = datad2$temp[unlist(c(peaks[1,]))], color='red')+
  #  geom_vline(xintercept = datad2$temp[unlist(c(peaks[2,]))], color='blue',linetype='dashed')+
  #  geom_vline(xintercept = datad2$temp[unlist(c(peaks[3,]))], color='green',linetype='dashed')
  
  
  
  #### covert "which" to temperature values 
  xx = as.matrix(expand.grid(1:nrow(peaks),1:ncol(peaks)))
  peaks[xx] = datad2$temp[peaks[xx] ]
  
  
  #### Second noise control
  ## Calcularte area
  peaks$raw.area = unlist(lapply(1:nrow(peaks), function(i){
    
    d = peaks[i,]
    d = datad2[datad2$temp>d$start & datad2$temp<d$end,]
    d$deriv[d$deriv<0]=0
    d$deriv = d$deriv - min(d$deriv)
    
    measureIntergral(d, xname='temp', yname='deriv')
    
  }))
  
  # Not more than three peaks, => the ones with the biggest area
  peaks = peaks[order(as.numeric(peaks$raw.area), decreasing = T),]
  
  #ggplot(datad2, aes(temp, deriv))+
  #  geom_point()+
  #  geom_vline(xintercept = unlist(c(peaks[1,1:5])), color='red')+
  #  geom_vline(xintercept = unlist(c(peaks[2,1:5])), color='blue',linetype='dashed')+
  #  geom_vline(xintercept = unlist(c(peaks[3,1:5])), color='green',linetype='dashed')
  
  n2 = ifelse(nrow(peaks)>=3, 3, nrow(peaks))
  peaks = peaks[1:n2,]
  
  peaks = cbind(peaks,do.call(rbind,lapply(1:nrow(peaks), function(i){
    x= peaks[i,]
    wx=datad2$temp>x$start & datad2$temp<x$end
    data.frame(sd.deriv2 =  sd(datad2$deriv2[wx]), 
               diff.deriv = max(datad2$deriv[wx])-min(datad2$deriv[wx]),
               mean.deriv = mean(datad2$deriv[wx]),
               mean.deriv2 = mean(datad2$deriv2[wx]))
  })))
  
  #bad = peaks$raw.area/max(peaks$raw.area) < 0.005 | peaks$sd.deriv2 < 0.01 | peaks$diff.deriv <0.01 | peaks$mean.deriv < 0.01
  bad = peaks$raw.area/max(peaks$raw.area) < 0.01 | peaks$sd.deriv2 < 0.01 | peaks$diff.deriv <0.01 | peaks$mean.deriv < 0.01
  
  # for low wt frequency the wt peak can be under bad threshold
  # Small verification
  if(bad[which.max(peaks$start)]){
    
    if(peaks$raw.area[which.max(peaks$start)]>0.01 &
       peaks$mean.deriv2[which.max(peaks$start)]< -0.01){
      bad[which.max(peaks$start)]=F
    }
    
  }
  
  
  peaks = peaks[!bad,]
  
  #ggplot(datad2, aes(temp, deriv))+
  #  geom_point()+
  #  geom_vline(xintercept = unlist(c(peaks[1,1:5])), color='red')+
  #  geom_vline(xintercept = unlist(c(peaks[2,1:5])), color='blue',linetype='dashed')+
  #  geom_vline(xintercept = unlist(c(peaks[3,1:5])), color='green',linetype='dashed')
  
  
  #ggplot(datad2, aes(temp, deriv3))+
  #  geom_point()+
  #  geom_vline(xintercept = unlist(c(peaks[1,1:5])), color='red')+
  #  geom_vline(xintercept = unlist(c(peaks[2,1:5])), color='blue',linetype='dashed')+
  #  geom_vline(xintercept = datad2$temp[d3.min], color='green',linetype='dashed')
  
  if(nrow(peaks)==0) stop("No peak found: Maybe an empty well?")
  
  ##############################################################################################
  #### Look for missed bumps (notably, very very slight mutant phase when less than 10% mutant)
  #### Use of deriv 3 and segmentation model
  if(nrow(peaks)<3 & nrow(peaks)>0){
    candidate = datad2$temp[d3.min]
    for(i in 1:nrow(peaks)){
      
      # If missed peak, we expect it to be between the begining/end and inflexion point of anothe rone
      c1 = candidate>peaks[i,]$start & candidate<peaks[i,]$inflex.bf
      c2 = candidate>peaks[i,]$inflex.af & candidate<peaks[i,]$end
      
      ## If possibly two peaks in the start => inflex.bf interval
      if(sum(c1)>0){
        # See if there are several slope in deriv 3
        wx = datad2$temp>peaks[i,]$start & datad2$temp< peaks[i,]$inflex.bf
        x=datad2$temp[wx]
        y=datad2$deriv3[wx]
        wx2 = 1:which.max(y)
        x=x[wx2]
        y=y[wx2]
        #plot(x,y)
        d=NULL
        d=data.frame(x=x,y=y)
        fit=lm(y~x,data=d)
        # Find breakpoints
        breakpoints = try(segmented(fit, seg.Z = ~x, npsi = 2)$psi[1:2,2])
        
        if(class(breakpoints)[1]=='try-error'){
          pass=FALSE
        }else{
          slp1=slp2=slp3=NULL
          slp1 = lm(y~x, d[d$x<breakpoints[1],])$coef[2]
          slp2 = lm(y~x, d[d$x>breakpoints[1] & d$x<breakpoints[2],])$coef[2]
          slp3 = lm(y~x, d[d$x>breakpoints[2],])$coef[2]
          
          # Noise control
          pass = slp3>10*slp2 & slp1>10*slp2 & sum(candidate[c1] > breakpoints[1] & candidate[c1] < breakpoints[2])>0
        }
        
        
        if(pass){
          # If pass, add the new 'peak' and adjust
          newtemp = c(peaks[i,]$start,mean(candidate[c1]))
          peaks[i,]$start = newtemp[2]
          toadd=NULL
          toadd = data.frame(start = newtemp[1],
                             inflex.bf = NA, 
                             peak=NA, 
                             inflex.af = NA,
                             end = newtemp[2])
          
          toadd = cbind(toadd, matrix(rep(NA, ncol(peaks)-ncol(toadd)), nrow=1))
          colnames(toadd) = colnames(peaks)
          
          peaks = rbind(peaks, toadd)
          
        }
        
      }
      
      
      
      ## If if possibly two peaks in the inflex.af => end interval
      if(sum(c2)>0){
        wx = datad2$temp < peaks[i,]$end & datad2$temp > peaks[i,]$inflex.af
        x=datad2$temp[wx]
        y=datad2$deriv3[wx]
        wx2 = which.max(y):length(y)
        x=x[wx2]
        y=y[wx2]
        #plot(x,y)
        d=NULL
        d=data.frame(x=x,y=y)
        fit=lm(y~x,data=d)
        breakpoints =try(segmented(fit, seg.Z = ~x, npsi = 2)$psi[1:2,2])
        
        if(class(breakpoints)[1]=='try-error'){
          pass=FALSE
        }else{
          slp1=slp2=slp3=NULL
          slp1 = lm(y~x, d[d$x<breakpoints[1],])$coef[2]
          slp2 = lm(y~x, d[d$x>breakpoints[1] & d$x<breakpoints[2],])$coef[2]
          slp3 = lm(y~x, d[d$x>breakpoints[2],])$coef[2]
          
          pass = slp3<10*slp2 & slp1<10*slp2 & sum(candidate[c2] > breakpoints[1] & candidate[c2] < breakpoints[2])>0
          
        }
        
        if(pass){
          newtemp = c(mean(candidate[c1]), peaks[i,]$end)
          peaks[i,]$start = newtemp[1]
          toadd=NULL
          toadd = data.frame(start = newtemp[1],
                             inflex.bf = NA, 
                             peak=NA, 
                             inflex.af = NA,
                             end = newtemp[2])
          
          toadd = cbind(toadd, matrix(rep(NA, ncol(peaks)-ncol(toadd)), nrow=1))
          colnames(toadd) = colnames(peaks)
          
          peaks = rbind(peaks, toadd)
          
        }
        
      }
      
    }
    
  }
  
  ##################################### 
  ### Finally, attribute genotypes to peaks
  peaks = peaks[order(peaks$start),]
  peaks = peaks[,1:5]
  
  peaks$phase.genotype=NA
  peaks$phase.genotype.score=1
  
  # First attribute it to the main peak based on teh provided melting range
  mainpeak = which.max(approx(x=datad2$temp, y=datad2$deriv, peaks$peak)$y)
  if(peaks$peak[mainpeak]>trange.mut[1] & peaks$peak[mainpeak]<trange.mut[2]){peaks$phase.genotype[mainpeak]="mut"}
  if(peaks$peak[mainpeak]>trange.wt[1] & peaks$peak[mainpeak]<trange.wt[2]){peaks$phase.genotype[mainpeak]="wt"}
  
  # If more than one peak
  if(!is.na(peaks$phase.genotype[mainpeak]) & nrow(peaks)>1){
    
    if(nrow(peaks)==2){
      # if after mut main peak, wt, if before het
      if( peaks$phase.genotype[mainpeak]=="mut"){
        peaks$phase.genotype[which(c(1,2)!=mainpeak)] = ifelse(which(c(1,2)!=mainpeak)==1, 'het', 'wt')}
      
      if(peaks$phase.genotype[mainpeak]=="wt"){
        # if before wt main peak, het, if before: should not happen, probably noise or peaks are shifter
        peaks$phase.genotype[which(c(1,2)!=mainpeak)] = ifelse(which(c(1,2)!=mainpeak)==1, 'het', NA) # should look at melting range to be sure differenciate het and mut
        peaks$phase.genotype.score[which(c(1,2)!=mainpeak)] = peaks$phase.genotype.score[which(c(1,2)!=mainpeak)] - ifelse(is.na(peaks$phase.genotype[which(c(1,2)!=mainpeak)]), 0.5, 0)
      }
    }
    
    if(nrow(peaks)==3){
      # If three peaks give the genotype with order c("het", "mut", "wt")
      #if genotype does not fit the melting range of the main peak, lower than he score
      score.decrease = ifelse(peaks$phase.genotype[mainpeak] == c("het", "mut", "wt")[mainpeak], 0, 0.5)
      peaks$phase.genotype.score = peaks$phase.genotype.score - score.decrease
      peaks$phase.genotype = c("het", "mut", "wt")

      
    }
    
  }
  
  #No space between phase
  if(nrow(peaks)>1){
    for(i in 1:(nrow(peaks)-1)){
      tjunction = mean(peaks$end[i], peaks$start[i+1])
      peaks$end[i] = peaks$start[i+1] = tjunction
    }
  }
  
  Taf = max(peaks$end)
  Tbf = min(peaks$start)
  return(list(peaks=peaks, melt.range = c(Tbf=Tbf,Taf=Taf)))
  
}









# Function that goes into peak.delimitation
peak.limits = function(inflex,data.melt,which.row.inflex,which.is.peak){
  
  i= which.row.inflex
  
  #Infor for the peak (peak, inflex point)
  thispeak=data.frame(start = ifelse(length(inflex$down[i-1])==0,NA,inflex$down[i-1]),
                      inflex.bf = inflex$up[i], 
                      peak=which.is.peak, 
                      inflex.af = inflex$down[i],
                      end = ifelse(length(inflex$up[i+1])==0,NA,inflex$up[i+1]))
  
  
  # Define start and end
  
  if(is.na(thispeak$start)){
    #x = data.melt$temp[1:thispeak$inflex.bf]
    #y = data.melt$deriv[1:thispeak$inflex.bf]
    
    #fit = lm(y~x)
    #breakpoint =segmented(fit, seg.Z = ~x, npsi = 1)$psi[2]
    
    tmax = data.melt$temp[thispeak$inflex.bf]
    tstart =  max(data.melt[data.melt$temp<tmax,]$temp[which.closest(data.melt$deriv2[data.melt$temp<tmax],0)],
                  data.melt[data.melt$temp<tmax,]$temp[which.min(data.melt$deriv[data.melt$temp<tmax])])
    
    thispeak$start= which.closest(data.melt$temp, tstart)
    
    
  }else{
    newstart=NULL
    start.peak = thispeak$start:thispeak$inflex.bf
    
    newstart = max(start.peak[which.closest(data.melt$deriv2[start.peak],0)],
                   start.peak[which.min(data.melt$deriv[start.peak])])
    
    
    if(thispeak$start < newstart) thispeak$start = newstart
    
  }
  
  
  
  if(is.na(thispeak$end)){
    #x = data.melt$temp[thispeak$inflex.af:length(data.melt$temp)]
    #y = data.melt$deriv[thispeak$inflex.af:length(data.melt$deriv)]
    #fit = lm(y~x)
    #breakpoint =segmented(fit, seg.Z = ~x, npsi = 1)$psi[2]
    
    tmin = data.melt$temp[thispeak$inflex.af]
    tend =  min(data.melt[data.melt$temp>tmin,]$temp[which.closest(data.melt$deriv2[data.melt$temp>tmin],0)],
                data.melt[data.melt$temp>tmin,]$temp[which.min(data.melt$deriv[data.melt$temp>tmin])])
    
    
    thispeak$end= which.closest(data.melt$temp, tend)
    
  }else{
    
    end.peak = thispeak$inflex.af:thispeak$end
    
    newend = NULL
    
    newend = min(end.peak[which.closest(data.melt$deriv2[end.peak],0)],
                 end.peak[which.min(data.melt$deriv[end.peak])])
    
    if(thispeak$end > newend)  thispeak$end = newend
    
  }
  
  
  thispeak
  
  
}












melting.range = function(data){
  # No duplicated temps
  data = data[!duplicated(data$temp),]
  #Smooth the curve
  data$fluo = predict(loess(fluo~temp, data, span=0.1), newdata=data)
  #plot(data$temp, data$fluo)
  datad = derivative(data, xname="temp", yname="fluo", smooth=0.1, return.df =T)
  #plot(datad$temp, datad$deriv)
  
  # Find temp where the derivative is at its max (the peak) and keep temp range +10 to -10
  tmax = datad$temp[which.max(datad$deriv)]
  temprange = tmax + c(-10,+10)
  datad = datad[datad$temp>temprange[1] & datad$temp<temprange[2],]
  
  #plot(datad$temp, abs(datad$deriv))
  
  # Do the absolute of the second derivative 
  datad$deriv2 = -derivative(datad, xname="temp", yname="deriv", smooth=0.1)
  ##plot(datad$temp, datad$deriv2)
  #plot(datad$temp, abs(datad$deriv2))
  #datad$deriv2 = abs(datad$deriv2)
  
  ## deriv second for background before and after melting
  bf.background = datad$deriv2[datad$temp<(tmax-9)]
  #sd(bf.background)
  af.background = datad$deriv2[datad$temp>(tmax+9)]
  #sd(af.background)
  
  #plot(datad$temp, abs(datad$deriv))
  
  wx = which((abs(datad$deriv2) > 10*mean(abs(bf.background)) & datad$temp<tmax & datad$deriv2>mean(bf.background))|(abs(datad$deriv2) > 10*mean(abs(af.background)) & datad$temp>tmax))
  #seg = wx[which(diff(wx)>1)]
  #pass = 0
  #i=1
  #while(pass==0 & i<=length(seg)& length(seg)>0){
  #  pass = sum(datad$deriv[1:seg[i]] >= mean(bf.background)*50, na.rm=T)
  #  i = i +1
  #}
  #i=i-2
  
  #datad = datad[wx[wx>ifelse(i>0, seg[i], 0)],]
  
  datad=datad[wx,]
  
  #unlist(lapply(split(datad, rep(1:ceiling(nrow(datad)/10), each=10)[1:nrow(datad)]), function(x) sd(x$deriv2)))
  
  #Temperature before and after melting
  return(data.frame(Tbf = min(datad$temp)-0.5, Taf = max(datad$temp)))
}





derivative=function(data, xname="temp", yname="normfluo", smooth=NULL, return.df =F){
  data <-data[order(data[,colnames(data)==xname]),]
  x = data[,colnames(data)==xname]
  y = data[,colnames(data)==yname]
  #if(!is.null(smooth)) y=predict(loess(y~x, span=smooth), newdata=x)
  #plot(x,y)
  #data = subset(data, between(temp, 78.5, 825))
  data$deriv = c(0,-diff(y)) / c(1, diff(x))
  #plot(x,data$deriv)
  if(!is.null(smooth)) data$deriv=predict(loess(deriv~temp, data, span=smooth), newdata=data)
  if(return.df){
    xy = data.frame(data$temp, data$deriv)
    names(xy)=c(xname,'deriv')
    return(xy)
  }else{return(data$deriv)}
}

############################################

find.extremes=function(y){
  sidi = sign(diff(y))
  sidi[sidi==0]= sidi[which(sidi==0) + 1]
  local.max=which(diff(sidi)==-2)+1
  local.min=which(diff(sidi)==2)+1
  return(list(local.min=local.min,local.max=local.max))
}

############################################

measureIntergral = function(data, xname="tempshift", yname="normfluo"){
  #data$id=paste0(data$well_pcr, data$plate_qpcr)
  x = data[,which(colnames(data)==xname)]
  y = data[,which(colnames(data)==yname)]
  data=data.frame(x,y)
  data=data[order(x),]
  data=subset(data, !is.na(x) & !is.na(y))
  diffx= diff(data$x)
  if(length(which(diffx==0))>0){
    data$y[which(diffx==0)]= (data$y[which(diffx==0)+1] + data$y[which(diffx==0)])/2
    data=data[-which(diffx==0),]
    diffx= diff(data$x)
  }
  data=data.frame(y1=data$y[-1], y2=data$y[-length(data$y)], h=diffx) 
  data$it= 0.5*(data$y1+data$y2)*data$h
  integral= sum(data$it)
  #len=aggregate(flu~id,data=data,length)
  #if(length(unique(len$flu))>1) print("WARNING:MissingTemp")
  #int = aggregate(flu~id,data=data,sum)
  #colnames(int)[which(colnames(int)=='flu')]='integral'
  return(integral)
}

############################################

tempshift=function(normd){
  #Find melting peaks and assign wild-type or mutant
  
  info = normd[[2]]
  peaks = peaks2= normd[[3]]
  normd = normd[[1]]
  
  
  peaks = subset(peaks,phase.genotype.score ==1 & phase.genotype  != 'het')
  meanpeaks=setNames(aggregate(peak~phase.genotype, peaks2, mean, na.rm=T), c('p','temp'))
  
  pwt=subset(peaks, phase.genotype=="wt")
  pwt=data.frame(id=pwt$id, peakwt=pwt$peak)
  pmut=subset(peaks, phase.genotype=="mut")
  pmut=data.frame(id=pmut$id, peakmut=pmut$peak)
  
  meandiff = mean(pwt[,2],na.rm=T) - mean(pmut[,2],na.rm=T)
  
  peaks=merge(pwt,pmut, all=T)
  peaks$peaksdiff=peaks$peakwt - peaks$peakmut
  peaks$shift=(peaks$peakwt + peaks$peakmut - sum(meanpeaks$temp))/2
  peaks$shift[is.na(peaks$shift)]  = ifelse(is.na(peaks$peakmut[is.na(peaks$shift)]),
                                            peaks$peakwt[is.na(peaks$shift)]-meanpeaks$temp[meanpeaks$p=="wt"],
                                            peaks$peakmut[is.na(peaks$shift)]-meanpeaks$temp[meanpeaks$p=="mut"])
  
  #peaks$shift2 = peaks$peaksdiff - meandiff
  #peaks$shift2[is.na(peaks$shift2)]=0
  normd=merge(normd,peaks)
  normd$tempshift = normd$temp - normd$shift
  #normd$tempshift2 = normd$tempshift - (normd$shift2/2)
  
  if(sum(colnames(normd)=='bk1')>0){
    
    startmelt = aggregate(bk1~id, normd, mean, na.rm=T)
    meanstartmelt= mean(startmelt[,2],na.rm=T)
    peaks=merge(peaks, startmelt)
    
    pref=as.vector(c(meanstartmelt,meanpeaks$temp[meanpeaks$p=="wt"]))
    normd=bind_rows(lapply(split(normd, normd$id), function(x){
      i = x$id[1]
      print(i)
      ptar = unlist(peaks[peaks$id==i,c(6,2)])
      if(length(ptar)==2){
        m=lm(pref~ptar)
        x$tempshift2= x$temp * m$coefficients[2] + m$coefficients[1]
      }else{x$tempshift2=x$tempshift}
      if(sum(is.na(x$tempshift2))==length(x$tempshift2)){x$tempshift2=x$tempshift}
      return(x)
    }))
    
  }
  
  
  return(list(norm=normd, melt.info=info, peaks = peaks2))
}





######################################


normalization.base = function(data, Tbf, Taf, correction=0){
  lmbf = lm(fluo~temp,subset(data, temp<Tbf & temp>(Tbf-1) ))
  lmaf = lm(fluo~temp,subset(data, temp<(Taf+1) & temp>(Taf)))
  
  #print(data$id[1])
  #l=nrow(data)
  
  pred=function(model,newdata, correction=0){
    b = model$coef[1]
    a=model$coef[2] + correction
    pr=a*newdata + b
    return(pr)
  }
  
  #L1=predict(lmbf, newdata=data.frame(temp=data$temp))
  #L0=predict(lmaf, newdata=data.frame(temp=data$temp))
  L1=pred(lmbf, data$temp, correction=correction)
  L0=pred(lmaf, data$temp,correction=0)
  
  data$normfluo = (data$fluo - L0)/(L1-L0)
  #data$normfluo[data$temp<Tbf]=1
  #data$normfluo[data$temp>Taf]=0
  #data$normfluo=data$normfluo*100
  minf = min(data$normfluo[data$temp>Tbf & data$temp<Taf], na.rm=T)
  maxf= max(data$normfluo[data$temp>Tbf & data$temp<Taf], na.rm=T)
  data$normfluo = ((data$normfluo -minf)*100)/(maxf-minf)
  data$normfluo[data$temp<Tbf]=100
  data$normfluo[data$temp>Taf]=0
  
  #ggplot(data, aes(temp, normfluo, color=id))+geom_point()+xlim(74,85)
  
  #print(length(data$normfluo)-l)
  
  return(list(data$normfluo, 
              breakpoints= setNames(c(Tbf,NA,NA,NA, Taf), c("Tbf", 'bk1', 'bk2', 'bk3', 'Taf')),
              slpbf=setNames(c(lmbf$coefficients[2:1]), c("slope", 'intercept')), 
              slpaf=setNames(c(lmaf$coefficients[2:1]), c("slope", 'intercept')), 
              slphet=setNames(c(NA,NA), c("slope", 'intercept'))))
}


######################################

#data = subset(rawd, id == "RP1#TP1A1")
#plot(data$temp, data$fluo)
#Tbf = 73
#Taf = 84

# From: High resolution melting curve analysis with MATLAB-based program
# Li, Huaizhong, Lan, Ruiting, Peng, Niancai, Sun, Jing, Zhu, Yong
# 2016


#C = 1
#a = -0.5
#t = seq(0,10,0.1)
#Tbf = 1; Taf = 9
#ft = C*exp(1)^(a*t)
#plot(t, ft)

#data = data.frame(temp = t, fluo = ft)

#lm(fluo~temp,subset(data, temp<Tbf & temp>(Tbf-1)))
#lm(fluo~temp,subset(data, temp<(Taf+1) & temp>(Taf)))

normalizationEBS=function(data, Tbf, Taf, correction=0){
  loglmbf = lm(log(fluo)~temp,subset(data, temp<Tbf & temp>(Tbf-1)))
  loglmaf = lm(log(fluo)~temp,subset(data, temp<(Taf+1) & temp>(Taf)))
  
  expcurve = function(loglm, x){
    e=exp(1)
    B=loglm$coefficients[1]
    beta = loglm$coefficients[2]
    k= e^B
    # y = k*(e^(beta*x)) 
    # (derivative) dy = beta*k*(e^(beta*x)) 
    dy=beta*k*(e^(beta*x)) 
    return(dy)
  }
  
  
  dbf=expcurve(loglmbf, Tbf) + correction
  daf=expcurve(loglmaf, Taf)
  
  a = (log(daf/dbf))/(Taf-Tbf)
  C= dbf/a
  
  data$normfluo=data$fluo - ( C*exp(1)^(a*(data$temp - Tbf)) )
  
  ggplot()+
    geom_point(data = subset(data, temp>Tbf+1 & temp<Taf-1), aes(temp, normfluo))
  
  #plot(data$temp, data$fluo)
  #plot(data$temp, data$normfluo)
  
  minf = min(data$normfluo[data$temp>Tbf & data$temp<Taf], na.rm=T)
  maxf= max(data$normfluo[data$temp>Tbf & data$temp<Taf], na.rm=T)
  data$normfluo = ((data$normfluo -minf)*100)/(maxf-minf)
  data$normfluo[data$temp<Tbf]=100
  data$normfluo[data$temp>Taf]=0
  
  return(data$normfluo)
  
}




############################


HRM = function(rawd, ampli=NULL){
  
  if(sum(colnames(rawd) == "id",na.rm=T)==0){
    if(sum(grepl("plate", colnames(rawd)), na.rm = T)==0){
      warning("Variable 'id' is missing (unique identifier per curve). Create id <- well_pcr. OK if only one plate")
      rawd$id = rawd$well_pcr
    }else{
      stop("Variable 'id' is missing (= unique identifier per curve). Add it before using the function")
    }
  }
  
  
  normd = lapply(split(rawd, rawd$id), function(x){
    
    print( x$id[1])
    
    #x=subset(rawd, id=="1#A1")
    #plot(x$temp, x$fluo)
    ID = x$id[1]
    x = x[!duplicated(x$temp),]
    x = x[order(x$temp),]
    p = try(peak.delimitation(x))
    
    if(class(p)[1]=='try-error'){print(paste0("Analyze failed for ", ID, " => Cause: ", p[[1]])); return(NULL)}else{
      
      Tbf = p[[2]][1]
      Taf = p[[2]][2]
      peaks = p[[1]]
      
      nn = try(normalization.base(x, Tbf=Tbf, Taf=Taf))
      if(class(nn)[1]=='try-error'){ x$normfluo=NA}else{x$normfluo = nn[[1]]}
      
      nEBS = try(normalizationEBS(x, Tbf=Tbf, Taf=Taf))
      if(class(nEBS)[1]=='try-error'){ x$normfluoEBS=NA}else{x$normfluoEBS = nEBS}

      
      peaks = cbind(peaks, do.call(rbind, lapply(1:nrow(peaks), function(i){
        melt.range = subset(x, temp > peaks[i,]$start & temp < peaks[i,]$end)
        melt.range = melt.range[c(1,nrow(melt.range)),]
        raw.decrease = melt.range$fluo[1] - melt.range$fluo[2]
        norm.decrease = melt.range$normfluo[1] - melt.range$normfluo[2]
        normEBS.decrease = melt.range$normfluoEBS[1] - melt.range$normfluoEBS[2]
        
        data.frame(raw.decrease=raw.decrease, norm.decrease =norm.decrease,normEBS.decrease=normEBS.decrease )
      })))
      
      weight = sapply(peaks$phase.genotype, function(x){
        if(!is.na(x)){
          c(0.5,1,0)[which(x==c('het','mut','wt'))]
        }else{ NA }})
      
      mutant.contribution = apply(peaks[,(ncol(peaks)-2):ncol(peaks)], 2, function(x){
        sum(x*weight) / sum(x)
      })
      
      
      
      smooth.fluo = predict(loess(fluo~temp, x, span=0.1), newdata=x)
      
      x$deriv =  derivative(data.frame(temp=x$temp, fluo=smooth.fluo), xname="temp", yname="fluo", smooth=0.1, return.df =F)
      
      #plot(x$temp, x$deriv)
      melt.info = cbind(data.frame(id=ID, Tbf=Tbf, Taf=Taf), t(as.matrix(mutant.contribution)))
      

     if(sum(peaks$phase.genotype=='het', na.rm=T)>0){
       phet = peaks[peaks$phase.genotype=='het',]
        if(!is.na(phet$peak)){
           melt.range.het = subset(x, temp > phet$peak & temp <phet$end)
        }else{
          melt.range.het = subset(x, temp > phet$inflex.bf & temp <phet$end)
        }
       
        melt.range.het = melt.range.het[c(1,nrow(melt.range.het)),]
        melt.info$het.decrease = (melt.range.het$normfluoEBS[1] - melt.range.het$normfluoEBS[2])/(melt.range.het$temp[2] - melt.range.het$temp[1])
     }else{
       melt.info$het.decrease = 0
     }
     
      
      
      peaks$id = ID
      
      return(list(norm=x, melt.info=melt.info, peaks.decomposition = peaks))
      
    }
    
  })
  
  
  normd = list(norm=do.call(rbind, lapply(normd, function(x) x[[1]])), 
               melt.info=do.call(rbind, lapply(normd, function(x) x[[2]])),
               peaks.decomposition = do.call(rbind, lapply(normd, function(x) x[[3]])))
  
  
  ## Some curve characteristic to compute badscore
  if(!is.null(ampli)){
    
    pslope = slopePlateau(ampli)
    normd[[2]] = merge(normd[[2]], pslope)
    
    
  }else{
    
    maxfluo = do.call(rbind, lapply(split(rawd, rawd$id), function(x){
      data.frame(id=x$id[1], maxfluo=max(x$fluo))}))
    
    normd[[2]] = merge(normd[[2]], maxfluo)
    
  }
  
  
  return(normd)
  
}





###########################################
############################################

subset.ID = function(normd, id, out=F){
  ID=id
  
  n = normd[[1]]
  x = normd[[2]]
  p = normd[[3]]
  
  if(out==T){
    n=n[!(n$id %in% ID),]
    x=x[!(x$id %in% ID),]
    p=p[!(p$id %in% ID),]
  }else{
    n=n[n$id %in% ID,]
    x=x[x$id %in% ID,]
    p=p[p$id %in% ID,]
  }
  
  normd[[1]] = n
  normd[[2]] = x
  normd[[3]] = p
  
  return(normd)
  
}






#################################
#################################

plot.melting.phase = function(normd,id,title="none"){
  
  require(gridExtra)
  
  ID=id
  
  n = normd[[1]]
  x = normd[[2]]
  p = normd[[3]]
  
  n=n[n$id==ID,]
  x=x[x$id==ID,]
  p=p[p$id==ID,]
  
  Tbf = x$Tbf
  Taf = x$Taf
  
  if(title=="none"){tID = id}else{tID=title}
  
  deriv = ggplot()+
    geom_rect(data=p, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, color=phase.genotype, fill=phase.genotype), alpha=0.5)+
    geom_point(data=subset(n, temp > Tbf -3 & temp < Taf + 3) , aes(temp, deriv), color='grey')+
    geom_point(data=subset(n, temp > Tbf & temp < Taf) , aes(temp, deriv))+theme_minimal()+
    ggtitle(tID)+
    xlab("Temperature")+
    ylab("Derivative")
  
  
  raw = ggplot()+
    geom_rect(data=p, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, color=phase.genotype, fill=phase.genotype), alpha=0.5)+
    geom_point(data=subset(n, temp > Tbf -3 & temp < Taf + 3) , aes(temp, normfluoEBS), color='grey')+
    geom_point(data=subset(n, temp > Tbf & temp < Taf) , aes(temp, normfluoEBS))+theme_minimal()+
    xlab("Temperature")+
    ylab("Normalized fluo. (EBS)")
  
  
 grid.arrange(deriv, raw)
 
 
  

  
}



plot.melting.phase2 = function(normd,ids,titles="none"){
  
  require(gridExtra)
  
  source("./basics_TP.R")
  
  ps = lapply(ids, function(id){
    
    i = which(ids==id)
    
    ID=id
    
    n = normd[[1]]
    x = normd[[2]]
    p = normd[[3]]
    
    n=n[n$id==ID,]
    x=x[x$id==ID,]
    p=p[p$id==ID,]
    
    Tbf = x$Tbf
    Taf = x$Taf
    
    if(titles[1]=="none"){tID = id}else{tID=titles[i]}
    
    deriv = ggplot()+
      geom_rect(data=p, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill=phase.genotype),color=NA, alpha=0.3)+
      geom_point(data=subset(n, temp > Tbf -3 & temp < Taf + 3) , aes(temp, deriv), color='grey')+
      geom_point(data=subset(n, temp > Tbf & temp < Taf) , aes(temp, deriv))+
      theme_Publication2()+
      scale_fill_manual(breaks = c("het", "mut", "wt"), values = c("darkorchid2","red3","steelblue2"))+
      ggtitle(tID)+
      xlab("Temperature")+
      ylab("Derivative")+
      theme(legend.position = "none")
    
    
    raw = ggplot()+
      geom_rect(data=p, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, color=NA, fill=phase.genotype),color=NA, alpha=0.3)+
      geom_point(data=subset(n, temp > Tbf -3 & temp < Taf + 3) , aes(temp, normfluoEBS), color='grey')+
      geom_point(data=subset(n, temp > Tbf & temp < Taf) , aes(temp, normfluoEBS))+
      theme_Publication2()+
      scale_fill_manual(breaks = c("het", "mut", "wt"), values = c("darkorchid2","red3","steelblue2"))+
      xlab("Temperature")+
      ylab("Normalized fluo. (EBS)")+
      theme(legend.position = "none")
    
    
    px=grid.arrange(deriv, raw)
    
    px
  })
  
  
  pp = grid.arrange(grobs=ps, nrow = 1)
  
  pp
  
}



######################################################
### Curve Quality control
### Attribute a score given three curve characteristic that I notice produce outliers
### 1) The maximum fluoresence (maxfluo) = raw fluo after amplification/ before melting
### 2) the slope at the heteroduplex melting peak: it increase when closer to 50% mutant (more heterozygote)
### => But when higher than expected => result in faster melting = more mutant
### 3) The slope at the end of PCR cycle: does it reach a plateau or not?
### Each time a curve is outside the 0.99 distribution for one of these charateristic: +1 to the score
### The higher score = the more absnormal curve 
badscore = function(melting.info){
  
  if("plateau_slope" %in% colnames(melting.info)){pslope = T}else{pslope=F}
  if("plate_qpcr_id" %in% colnames(melting.info)){multiplate = T}else{multiplate=F}
  
  if(!multiplate){
    warning("No plate_qpcr_id variable, assuming there only one qpcr plate. If not, it's better to add a plate_qpcr_id variable")
    melting.info$plate_qpcr_id = 'A'
  }
  
  melting.info = do.call(rbind, lapply(split(melting.info, melting.info$plate_qpcr_id), function(x){
    maxfluorange = mean(x$maxfluo, na.rm=T) + c(-1,1)*sd(x$maxfluo, na.rm=T)*2.33 # z score 0.99
    hetrange = mean(x$het.decrease, na.rm=T) + c(-1,1)*sd(x$het.decrease, na.rm=T)*2.33
    
    if(pslope){prange = mean(x$plateau_slope, na.rm=T) + c(-1,1)*sd(x$plateau_slope, na.rm=T)*2.33}else{
      prange = rep(0, length(hetrange))
    }
    
    
    badscore = cbind((x$maxfluo < maxfluorange[1] | x$maxfluo > maxfluorange[2]),
                     (x$het.decrease < hetrange[1] | x$het.decrease > hetrange[2]),
                     (x$plateau_slope < prange[1] | x$plateau_slope > prange[2]))
    
    badscore = apply(badscore, 1, sum, na.rm=T)
    
    x$badscore = badscore
    x
    
  }))
  
  if(!multiplate){melting.info = melting.info[, -which(colnames(melting.info)=='plate_qpcr_id')]}
  
  melting.info
  
}




#############################################################
##### Plot bad



plot.bad = function(meltinfo, design = NULL, calibname = "calib", highlight = NULL){
  
  require(gridExtra)
  
  if(!is.null(design)){
    meltinfo = try(merge(meltinfo, design))
  }
  
  
  if(!('plate_qpcr_id' %in% colnames(meltinfo))){
    meltinfo$plate_qpcr_id='plate1'
  }
  
  meltinfo = subset( meltinfo, !is.na(plate_qpcr_id))
  mhl = subset(meltinfo, id %in% highlight)
  
  
  if(calibname %in% colnames(meltinfo)){
    mcalib = subset(meltinfo, !is.na(calib))
    meltinfo = subset(meltinfo, is.na(calib))}else{
      mcalib = meltinfo[0,]
    }
  
  
  if(!('maxfluo' %in% colnames(meltinfo))){pmaxfluo = NULL}else{
    pmaxfluo =  ggplot()+
      facet_grid(plate_qpcr_id~., scale='free')+
      theme_minimal()+
      geom_point(data = mcalib, aes(normEBS.decrease, maxfluo), color='black')+
      geom_point(data = meltinfo, aes(normEBS.decrease, maxfluo, color=badscore))+
      geom_point(data = mhl, aes(normEBS.decrease, maxfluo), color='green')+
      scale_color_gradient(low='blue', high='red')+theme(legend.position = "none")+
      ggtitle("Fluorescence before melting")
  }
  
  
  if(!('plateau_slope' %in% colnames(meltinfo))){pplat = NULL}else{
    pplat = ggplot()+
      facet_grid(plate_qpcr_id~., scale='free')+
      theme_minimal()+
      geom_point(data = mcalib, aes(normEBS.decrease, plateau_slope), color='black')+
      geom_point(data = meltinfo, aes(normEBS.decrease,plateau_slope, color=badscore))+
      geom_point(data = mhl, aes(normEBS.decrease, plateau_slope), color='green')+
      scale_color_gradient(low='blue', high='red')+theme(legend.position = "none")+
      ggtitle("Slope at amplification plateau")
  }
  
  if(!('het.decrease' %in% colnames(meltinfo))){phet = NULL}else{
    phet = ggplot()+
      facet_grid(plate_qpcr_id~., scale='free')+
      theme_minimal()+
      geom_point(data = mcalib, aes(normEBS.decrease, het.decrease), color='black')+
      geom_point(data = meltinfo, aes(normEBS.decrease,het.decrease, color=badscore))+
      geom_point(data = mhl, aes(normEBS.decrease, het.decrease), color='green')+
      scale_color_gradient(low='blue', high='red')+theme(legend.position = "none")+
      ggtitle("Slope at heterozygous melting phase")
  }
  
  
  
  plist = list(pmaxfluo, phet, pplat)
  plist = plist[!unlist(lapply(plist, is.null))]
  do.call("grid.arrange", c(plist, ncol=length(plist)))
  
}
