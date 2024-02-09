setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/")

generalpath = "./Data/rec1EE/SLiM/" # Change the path (sufficient for this script if you keep the same directory structure)
# + You also need to change in the SliM scripts the Setwd line to your your path to transitory

library(data.table)
library(dplyr)
library(ggplot2)
library(parallel)


armbias =  F # If T: No fitness (fitness) loci sampled within the chromosomal centers.

transitpath = paste0(generalpath, "transitory/") # Some directory where the transitory file can be saved
slimpath = paste0(generalpath, "scripts/") # Where is the .slim file containing simulation script

source(paste0("./Data/Experimental_Evolution_R/EE_Rhom/Pool_seq_R_clean/","utils_poolseq_analysis.R"))

######################################################
######### IMPORT CEMEE GENOTYPE DATA #################

cemeesnps <- fread(paste0(generalpath, "input_data/CeMEEv2_RIL_snps_ws245.csv.gz"))
cemeegt <- fread(paste0(generalpath, '/input_data/CeMEEv2_RIL_geno.csv.gz'))
cemeegt = as.matrix(cemeegt)
cemeegt = cemeegt[,grepl('A6140',colnames(cemeegt))]

###############################################
######### IMPORT LINKAGE MAPS #################

#=> Code for to obtain smoothmaps at the end of the scipt
load(paste0(generalpath, "input_data/smoothmaps.Rdata"))

### Subset the simulated chromsomes
#tarchrom = target chromosome
tarchrom = c("I", "II")
if(armbias==F){
  # Keep the target chromosomes
  tarchrom = cemeesnps$chrom %in% tarchrom
}else{
  # Keep only the arms of the target chromosomes
  library(readxl)
  boundaries <- read_excel(paste0(generalpath, "domain_boundaries.xlsx"))
  tarchrom=unlist(lapply(tarchrom, function(chr){
    center = subset(boundaries, domain=='center' & chrom==romtonum(chr))
    which(cemeesnps$chrom == chr & (cemeesnps$pos<center$start_pos | cemeesnps$pos>center$end_pos))
  }))
  
}

cemeegt = cemeegt[tarchrom,]
cemeesnps = cemeesnps[tarchrom,]

#Define a line as the EEV1403 line used to introgress the rec-1 mutant
# It will represent half of the chromosome in the ancestral population
recline = cemeegt[,colnames(cemeegt)=='A6140L234']
cemeegt = cemeegt[, -which(colnames(cemeegt)=='A6140L234')]



#################################################################
############ FUNCTIONS ##########################################
#################################################################

####################################################################
#### I) Functions to create the initial populations (sampled CeMEE genotype with/without burn-in) 
#### and annex files (files containing selection coeficient and recombinationr ates at each locus)


### The main function is "architechture" and call most of other function
# It creates the initial population and auxiliary files (rrmap; selection) a call SLiM for the burn-in if needed. 
# First,  the genotype matrix (gmatrix; contain genomes of the populationloci x haploid genotype, ~  )
# and auxiliary files are created on R (file containing selection coeficient and recombinationr ates at each locus)
# (= > call sample.loci, make.rrmap; make.scoef; make.genomatrix). 
# It also transforms the gmatrix in an initial population in a format readable by SLiM 
# (= > call genoSLiM, gvar.correction). 
# and run a burn-in if necessary (call burnslim,parse_slim,. gvar.correction)
# Note: The "gid", a random identifier is attributed to the gmatrix/initial pop

#### architechture function is itself called by INITIAL_POPS function (at the end of this section)
# It calls architechture in multiple iteration to create multiple gmatrix in accordance to the different parameters and replication






# Samples fitness loci among SNV position (in cemeesnps) and stores info in a data.frame
# n loci per chromosome = lociperchr
#return:
#    chrom   pos     row.in.cemeesnps    type
#1      I   717890             <NA>          rec
#2      I  1649959             4640 fitness.loci
#3      I  1807883             5407 fitness.loci
sample.loci = function(lociperchr, cemeesnps, recpos = 717890, recchrom='I'){
  # Toy parameters:
  #lociperchr=20 # The number of loci per chromosome (fixed)
  #cemeesnps # the snps of the cemee (df of the same name loaded above)
  # recpos = 717890 # rec-1 position (bp); NULL if monomorphic simulations
  # recchrom='I' # chromosome containing rec-1 
  
  randomloci = unlist(lapply(split(1:length(cemeesnps$chrom), cemeesnps$chrom), function(x){
    sample(x, lociperchr)
  }))
  
  loci = randomloci[order(randomloci)]
  loci = data.frame(chrom = cemeesnps$chrom[loci], 
                    pos = as.numeric(cemeesnps$pos[loci]), 
                    row.in.cemeesnps = loci,
                    type = "fitness.loci")
  
  loci$pos=as.numeric(loci$pos)
  
  # Add rec-1 locus if necessary
  if(!is.null(recpos)){
    
    loci = rbind(loci, c(recchrom, recpos, NA, "rec" ))
    loci = loci[order(romtonum(loci$chrom), as.numeric(loci$pos)),]
  }
  
  row.names(loci) = 1:nrow(loci)
  
  loci
  
}



#loci = sample.loci(lociperchr=20, cemeesnps, recpos = 717890, recchrom='I')


#  Create a rrmap.txt file with the recombination rate info for slim
# row i corresponds locus i: 3 columns: "wt", "mut", 'chrom' (rr wild-type, rr mutant, rr chromosome)
# Estimate rr from smoothmaps
make.rrmap = function(loci, smoothmaps, path=NULL, cemeesnps, id=NULL){
  # Toy parameters:
  # loci = loci # The data.frame returned by sample.loci
  # smoothmaps = smoothmaps # the df of the same name containing rec-1 wt and mut genetic positions
  # path = where to save (.txt) the rrmap # If NULL, simply returned by the function
  #cemeesnps # the snps of the cemee (df of the same name loaded above)
  # id = "gujk5k6a6bw" a unique identifier # During simulations, it corresponds to the "gid"
  
  rweight = do.call(rbind, lapply(unique(loci$chrom), function(chr){
    
    rr = do.call(cbind, lapply(c('wt', 'mut'), function(recgeno){
      
      l= subset(loci, chrom==chr)
      g=subset(smoothmaps, chrom==romtonum(chr) & rec==recgeno)
      genetic = approx(g$pos, g$smoothgenetic, xout = l$pos, rule=2)$y
      r = c(0,diff(genetic)) #The first loci value won't be taken into account in SLiM
      r
    }))
    
    colnames(rr)=c('wt', 'mut')
    
    rr=cbind(rr,chrom=rep(as.numeric(romtonum(chr)), nrow(rr)))
    
    #Warning: the order: "wt", "mut", "chrom is important in the SLiM code
    rr
  }))
  
  loci$gpos.wt = as.numeric(rweight[,1])
  loci$gpos.mut = as.numeric(rweight[,2])
  
  #The recombination rates between sampled loci according the smoothed map
  if(!is.null(path)){
    # The rr matrix wt,mut,chrom and the genomatrix locixlines need to be saved in a txt file
    if(is.null(id)){xs =""}else{xs="_"}
    write.table(t(rweight), file = paste0(path, "rrmap",xs,id,".txt"), sep = ",", col.names = F, row.names = F)
    return(loci)
  }else{
    return(loci)
  }
}






#  Create a selection.txt file with the selection and dominance coef info for slim
# row i corresponds locus i: 3 columns: "wt", "mut", 'chrom' (rr wild-type, rr mutant, rr chromosome)
# s are sampled from a normal distribution (sd does not really matter because s are adjusted after with gvar.correction)
# If dominance = T: h are sampled from uniform distrib; if F: codominance (h=0.5)
  
make.scoef=function(loci, s.sd=1, path=NULL, id=NULL, dominance =T, neutralProportion = 0){
  
  # Toy parameters:
  # loci = loci # The data.frame returned by sample.loci
  # s.sd = 1 # sd of the normal distribution used to sampled s (sd does not really matter because s are adjusted after with gvar.correction)
  # path = where to save (.txt) the rrmap # If NULL, simply returned by the function
  # id = "ghjk666hjk" a unique identifier # During simulations, it corresponds to the "gid"
  # dominance = F 
  # neutralProportion = 0, the proportion of neutral loci (0 to 1)
  
  nloci = nrow(loci)
  s = rnorm(nloci, sd = s.sd)
  s[sample(1:nloci, neutralProportion*nloci)] = 0
  
  if(dominance == T){h = runif(nloci)}else{h=rep(0.5, nloci)}
  
  
  s=cbind(s=s, h=h)
  
  # Selection = 0 at the rec-1 locus if present
  recl = loci$type == "rec"
  s[recl,1] = 0
  s[recl,2] = 0.5
  
  
  if(!is.null(path)){
    # The rr matrix wt,mut,chrom and the genomatrix locixlines need to be saved in a txt file
    if(is.null(id)){xs =""}else{xs="_"}
    write.table(t(s), file = paste0(path, "selection",xs,id,".txt"), sep = ",", col.names = F, row.names = F)
  }else{
    loci = cbind(loci, s)
    return(loci)
  }
}






# This function create a Genotype Matrix "genomatrix" : loci x haploid genomes
#  with 0 = ref , 1 = alt
# WARNING: The EEV1401 (here the L234) lines must be the last line
make.genomatrix = function(loci, N=1000, path=NULL, cemeegt, recline, id=NULL, burn=F){
  #loci = result of sample.loci
  #N = population size
  #Path = where to save the genomeatrix (.txt). If NULL, just return the genomatrix as a result of this function
  # cemeegt = cemeegt, the genotype matrix corrsponding to cemee lines
  # recline = recline, one column of the cemeegt chosen to represent half of haplotype (mimick the construction of the R populations)
  # id = an identifier for the genotype matrix that will be put in the name of the .txt file
  # before_rec= position before the rec-1 locus in loci; result of the function whereisrec
  # burn = F, = sample alleles at rec-1 lous
  # burn = T, rec-1 = always 0 for the burn-in. The allele will be sampled after the burn-in
  
  # Create a genomatrix from cemeegt, with [loci, n column from cemeegt]
  wloci = as.numeric(loci$row.in.cemeesnps)
  genomatrix = cbind(cemeegt[wloci, sample(1:ncol(cemeegt), N, replace = T)], recline[wloci])
  genomatrix = round(genomatrix)
  
  
  #Then we want to change 0/1=ref/alt to 0/1=neutral/selected
  # The beneficial allele is assigned randomly to the ref of alt
  # Note: probably wrong assumption because N2 allele (ref) are often beneficial
  switch = sample(c(0,1), nrow(genomatrix), replace = T)
  genomatrix = abs(genomatrix - switch)
  
  # Add the rec-1 locus if needed
  
  if(sum(loci$type=="rec")>0){
    
    recl = which(loci$type=="rec")
    
    if(burn){
      genomatrix[recl,] = 0
    }else{
      genomatrix[recl,] = sample(c(0,2), N+1, replace=T)
    }
    
    
  }
  
  
  
  if(!is.null(path)){
    # The rr matrix wt,mut,chrom and the genomatrix locixlines need to be saved in a txt file
    if(is.null(id)){xs =""}else{xs="_"}
    write.table(genomatrix, file = paste0(path, "genomatrix",xs,id,".txt"), sep = ",", col.names = F, row.names = F)
    
    
  }else{
    return(genomatrix)
  }
}


# This function is called after the creation of the gmatrix and annex files 
# The goal of this function is to scale the selection coeficient of the selection_ID.txt file to match the wanted fitness variance
# This is done by measuring individuals' genetic value (sum of selection coeficient, = log fitness when no epistasis)
# Then, a scaling factor is approximated given the desired beta and fitness variance:
# This is done by calculating population fitness variance given beta with a series of scaling factor
# and then intrepolating the one corresponding to the desired fitness varaince
# If there is a burn-in , the fitness variance will be corrected after the burn-in.

gvar.correction = function(fit.vars, gid, path, betas, N, loci, gvalue = NULL, addinfilename = NULL, selectionfile=NULL){
  
  # fit.vars: the different values genetic variance for fitness (if several, will create several corrected selection_ID.txt files for a given gmatrix )
  # gid: unique identifier
  # path: path containing the selection files
  # betas:  a parameter that define the curvature of the fitness function during evolution (0 = no epistasis; 1= stabilizing selection)
  # N popsize (corresponds to Nbase in other funtions)
  # loci: output of sample.loci
  # gvalue: the vector of gvalue (sum of indiviual coeficient). If not NULL, skip their calculation from the gmatrix. This is notably used when the function is called after the burn-in.
  # addinfilename: something to add in the file name or NULL
  # selectionfile: only use after when the function is called after burn-in
  
  
  fit.vars = unique(fit.vars)
  betas = unique(betas)
  ##### Estimate the add genetic variance and adjust to get the one wanted
  
  if(is.null(gvalue)){
    
    ### Reupload selection coeficients
    ss = as.matrix(fread(paste0(path, "selection_",gid,".txt"), header = F))
    
    ### Reupload the gmatrix to estimate population fitness
    gmatrix = as.matrix(fread(paste0(path, "genomatrix_",gid,".txt"), header = F))
    
    if(N == ncol(gmatrix)-1){
      ### last column is the rec1 line
      gmatrix = cbind(gmatrix[,1:N], gmatrix[,rep(N+1, N)])
    } ### else ncol should be 2N
    
    
    ### Allow for free recombination between chromosome
    gmatrix = do.call(rbind, lapply(split(as.data.frame(gmatrix), loci$chrom), function(x){
      x[,sample(1:(N*2),N*2)]
      as.matrix(x)
    }))
    
    ### Sum the genotype from two haploid genomes = 0,1,2 = ref, het, hom
    gmatrix = gmatrix[,1:N] +  gmatrix[,(N+1):(2*N)]
    
    ### Replace genotype by their selection coefficient and h
    
    hh = as.matrix(ss[2,], ncol=1)[,rep(1, N)]
    gmatrix[gmatrix==1] = hh[gmatrix==1]
    gmatrix[gmatrix==2] = 1
    gmatrix =  gmatrix * as.matrix(ss[1,], ncol=1)[,rep(1, N)]
    # genetic value for fitness (sum of selection coeficient = log fitness)
    gvalue = apply(gmatrix, 2, sum)
    #hist(exp(gvalue))
  }else{
    ss = as.matrix(fread(paste0(path, selectionfile), header = F))
  }
  
  
  # Approx the wanted g variance 
  
  xx=lapply(betas, function(beta){
    
    scalefactors = seq(0.0001, 1, 0.005)
    fitvars = unlist(lapply(scalefactors, function(s){
      g = gvalue*s
      g = g - min(g)
      grel = g/max(g)
      grel = g - beta*g^2 
      g = grel * max(g)
      var(g)
    }))
    
    fitvars = cbind(scalefactors, fitvars)
    fitvars = fitvars[!is.infinite(fitvars[,2]),]
    #plot(fitvars[,1], fitvars[,2])
    
    scalefactor = approx(fitvars[,2],fitvars[,1], xout=fit.vars)[[2]]
    
    
    for(i in 1:length(fit.vars)){
      newss = ss
      newss[1,] = newss[1,]*scalefactor[i]
      write.table(newss, file = paste0(path, "selection_",addinfilename, "fitvar",fit.vars[i],'_beta',beta,'_', gid,".txt"), sep = ",", col.names = F, row.names = F)
    }
    
  })
  
}



# Transform the genomatrix.txt file in a slim.txt file (Use SLiM for that: call genomatrixToSlim.slim)
# the pop in slim.txt file can be reuploaded into slim after, 
# and it's much faster than recreating the genomatrix from the genomatrix.txt file
genoSLiM = function(slimpath, rrpath, genopath, outpath, nloci, N){
  #path: path to genomatrixToSlim.slim
  #rrpath: path to rrmap_ID.txt
  # outpath: output path
  # nloci the total number of loci (nchrom * lociperchrom +1 if rec-1 locus)
  # N: initial pop size (corresponds to Nbase in come other functions)
  
  param_rr = paste0("-d rrpath=\"'", rrpath, "'\"")
  param_geno = paste0("-d genomatrixpath=\"'", genopath, "'\"")
  param_out = paste0("-d outpath=\"'", outpath, "'\"")
  param_nloci =  paste0("-d nloci=",nloci)
  param_N =  paste0("-d N=",N)
  
  # Command to write in the terminal to call slim
  slimcmd =paste("slim",
                 param_rr,
                 param_geno,
                 param_out,
                 param_N,
                 param_nloci,
                 slimpath,
                 " ")
  
  # Call slim in the terminal
  slim_out = system(slimcmd, intern=T)
  
  #return(paste0("saved ", outpath))
}





#Call slim to do the burn-in. 
#“stab” when beta = 1 (negative epistasis) and drift when no selection (as a control). 
#Both burn-in types can be done or only one. If both are done, the simulation will be run for each one.  
#Return info about fitness variance
burnslim = function(nrun, sfile, selfingrate, beta,nloci, N,rrfile,initialpopfile,outpoppath, slimpathfile, afc_outid=NULL){
  
  ## Toy parameter
  #nrun=1
  #scoef="\"0.1 0.2 0\""
  #selfingrate=0.2
  #nloci=15
  #N=1000
  #rec='wt'
  #rrfile="rrmap.txt"
  #genofile='genomatrix.txt'
  #slimfile="Data/Experimental_Evolution_R/EE_Rhom/SLiM/uniform/uniform_selection.slim"
  
  param_selfrate = paste0("-d self_rate=",selfingrate)
  param_nloci =  paste0("-d nloci=",nloci)
  param_N =  paste0("-d N=",N)
  param_beta =  paste0("-d beta=",beta)
  param_rr = paste0("-d rrpath=\"'", rrfile, "'\"")
  param_s = paste0("-d spath=\"'", sfile, "'\"")
  param_geno = paste0("-d initialpoppath=\"'", initialpopfile, "'\"")
  param_out = paste0("-d outid=\"'", afc_outid, "'\"")
  param_outpop = paste0("-d outpath=\"'", outpoppath, "'\"")
  
  slimcmd =paste("slim",
                 param_selfrate,
                 param_beta,
                 param_s,
                 param_nloci,
                 param_N,
                 param_rr,
                 param_geno,
                 param_out,
                 param_outpop,
                 slimpathfile,
                 " ")
  
  
  slim_out = system( slimcmd,intern=TRUE)
  
  slim_out
}


#This function parse the slim print output (from burnslim or runslimpolym) into a data.frame
parse_slim = function(slim_out, full.numeric=T){
  
  #slim_out=outputs[[1]]
  seed = as.numeric(slim_out[2])
  slim_out = slim_out[(which(grepl("#OUT:", slim_out))+1):length(slim_out)]
  slim_out = strsplit(slim_out, " ")
  outnames = slim_out[[1]]
  outnames=tstrsplit(outnames, '\"')[[2]]
  
  # Remains of another version of simulation
  # In the current state, full.numeric should always be T
  if(full.numeric){
    values = matrix(as.numeric(unlist(slim_out[2:length(slim_out)])), ncol=length(outnames),byrow=T)
  }else{
    
    values = do.call(rbind, lapply(2:length(slim_out), function(i){
      
      if(grepl('\"', slim_out[[i]][1])){
        tstrsplit(slim_out[[i]], '\"')[[2]]
      }else{
        as.numeric(unlist(slim_out[[i]]))
      }
      
    }))
    
  }
  values = as.data.frame(values)
  colnames(values)=outnames
  
  #values$generation = (1:(nrow(values)*2))[((1:(nrow(values)*2))%% 2)==1]
  values$runseed = seed
  values
  
}







# This function creates the initial population and auxiliary files (rrmap; selection) a call SLiM for the burn-in if needed. 
# First,  the genotype matrix (gmatrix; contain genomes of the populationloci x haploid genotype, ~  )
# and auxiliary files are created on R (file containing selection coeficient and recombinationr ates at each locus)
# (= > call sample.loci, make.rrmap; make.scoef; make.genomatrix). 
# It also transforms the gmatrix in an initial population in a format readable by SLiM 
# (= > call genoSLiM, gvar.correction). 
# and run a burn-in if necessary (call burnslim,parse_slim,. gvar.correction)
# Note: The "gid", a random identifier is attributed to the gmatrix/initial pop
architechture = function(slimpath, transitpath, lociperchr, N=1000, cemeesnps, smoothmaps, fit.vars,betas,
                         recpos = 717890, recchrom="I", burn=F, burn.fitvar=0.15, burntype = c("stab", 'drift'),dominance =F){
  
  Nbase = N
  if(recpos == "not_explicit") recpos = NULL
  
  nchrom = length(unique(cemeesnps$chrom))
  nloci = nchrom*lociperchr
  
  ### Make the genotype matrix (nloci x (cemme lines + L234 line))
  
  #1) random identifier to tag the genomatrix (needed if parallelisation)
  gid = paste(sample(c(letters, 1:9), 10), collapse='')
  while(sum(grepl(gid,list.files(transitpath)))>0){
    gid = paste(sample(c(letters, 1:9), 10), collapse='')
  }
  
  system(paste0("mkdir ",transitpath,gid))
  transitpath2 = paste0(transitpath,gid,"/")
  
  #2) make the files for slim in the path 
  loci = sample.loci(lociperchr=lociperchr, cemeesnps, recpos, recchrom) #sample loci in cemeesnps and store them in a df
  loci = make.rrmap(loci=loci, smoothmaps=smoothmaps, path=transitpath2, cemeesnps=cemeesnps, id=gid) # Add recombination rates (rr) info in df + save rr in .txt file for SLiM
  make.scoef(loci=loci, s.sd=1, path=transitpath2, id=gid, dominance =dominance) # draw selection coeficient and save them as a .txt file
  make.genomatrix(loci, N=Nbase, cemeegt=cemeegt, recline=recline,path=transitpath2, id=gid, burn=burn) # Create a genotype matrix for SLiM to create an initial population
  
  # 3) Calculate the selection needed to get some genetic variance for fitness and correct the selection file
  
  if(burn){fv = burn.fitvar;bb = 1}else{fv = fit.vars; bb = betas}
  gvar.correction(fit.vars=fv , 
                  gid=gid, path=transitpath2, 
                  betas=bb, N=Nbase, loci=loci)
  
  # 4) Create the initial population in slim format
  
  nloci = nrow(loci)
  
  # The genomatrix files created above is going to read by SliM to generate an initial population in the good format for slim
  
  genoSLiM(slimpath=paste0(slimpath,"genomatrixToSlim.slim"),
           rrpath=paste0(transitpath2, "rrmap_",gid,".txt"),
           genopath=paste0(transitpath2, "genomatrix_",gid,".txt"),
           outpath=paste0(transitpath2, "initialPop_",gid,".txt"),
           N=Nbase, nloci=nloci)
  
  
  
  system(paste0("rm ",transitpath2,"genomatrix_",gid,".txt")) # delete genomatrix file, won't be needed anymore

  
  
  # 5) Now, if there is a burn in phase:
  # Note: rec-1 locus not explicitely modeled during burn-in, just a fixed reference neutral locus
  
  if(burn == T){
    
    OUT=do.call(rbind, lapply(burntype, function(type){
      
      outburn = burnslim(sfile=paste0(gid, "/selection_fitvar0.15_beta1_",gid,".txt"), # s not taken into account if burn==drift
                         selfingrate=0.1,
                         beta=ifelse(type=='stab',1,0), #if 0 the s will also be set to 0 in the .slim
                         nloci=nloci,
                         N=N,
                         rrfile= paste0(gid,"/rrmap_",gid,".txt"),
                         initialpopfile = paste0(gid,"/initialPop_",gid,".txt"),
                         outpoppath = paste0(gid,"/initialPop_burnin_",type, '_',gid,".txt"),
                         slimpathfile=paste0(slimpath,"negative_epistasis_burn.slim"), 
                         afc_outid=NULL)
      
      outburn = parse_slim(outburn)
      outburn$run=1
      outburn$burn_in = "burn"
      outburn$burntype = type
      
      print('burn_in_done')
      
      # Recalculate the selection coefficient for the wanted fit.var taking account the beta 
      #ss = as.matrix(fread(paste0(transitpath2, "selection_fitvar0.15_beta1_",gid,".txt"), header = F))
      
      # Estimate the fitness variance given the additive genetic variance var.gabs
      
      
      
      #  Draw a distribution of "genetic value" given the genetic variance after the burn-in
      # This is treated as if it was the real gvalue distrubution of the populaion
      # calculate fitness:
      gdis = rnorm(1000, 0, sqrt(outburn$var.gabs[outburn$generation==99]))
      
      gvar.correction(fit.vars=fit.vars , 
                      gid=gid, path=transitpath2, 
                      betas=betas, N=Nbase, loci=loci,
                      gvalue = gdis, addinfilename = paste0("burn",type, "_"),
                      selectionfile = paste0("selection_fitvar0.15_beta1_",gid,".txt"))
    }))
    
    system(paste0("rm ",transitpath2,"selection_fitvar", burn.fitvar, "_beta1_",gid,".txt"))
    system(paste0("rm ",transitpath2,"initialPop_",gid,".txt"))
    
  }
  
  system(paste0("rm ",transitpath2,"selection_",gid,".txt"))
  
  
  loci$gid = gid
  return(loci)
}


# Calls architechture in multiple iteration to create multiple gmatrix in accordance to the different parameters and replication
INITIAL_POPS = function(ngmatrix, parameters, slimpath, transitpath, cemeesnps=cemeesnps,
                        smoothmaps=smoothmaps, burn=F, burntype = c("stab", 'drift'), burn.fitvar=0.15,
                        recchrom="I", Nbase = 1000,dominance =F){
  #Toy parameters:
  # parameters = paramdf.polym(nloci=c(300),N=c(1000), fit.vars = c(0.02), recpos =c(717890), selfing=c(0, 0.5))
  #ngmatrix = 2
  #burn = T
  #burntype = "drift"
  #burn.fitvar=0.15
  #recchrom = "I"
  #Nbase = 1000
  #dominance = 1000
  
  
  # gparam = df with parameters that require to create a new genotype matrix
  # i.e. the position of rec-1 and the number of loci
  # => ex: different fit.var can be done for the same genotype matrix
  gparam = unique(parameters[, c("recpos", "lociperchr")])
  
  do.call(rbind, lapply(1:nrow(gparam), function(i){
    #for each set of unique paramters:
    recpos = gparam[i,]$recpos  #position of rec-1
    lociperchr = gparam[i,]$lociperchr # number of loci per chromosome
    
    #px = subset of paramters corresponding for this position of rec-1 and this number of loci
    px = parameters[parameters$recpos==recpos & parameters$lociperchr==lociperchr, ] 
    
    print(paste0("Prepare ",ngmatrix," gmatrix with lociperchr=",lociperchr," and recpos=", recpos))
    
    
    
    # Create gmatrix correspnding to gparam[i,] (times ngmatrix)
    gidloci = lapply(1:ngmatrix, function(ng){
      
      # architechture function create all the gmatrix, the file containing recombination maps and selection coeficients 
      #info = list[gid=identifier of gamtrix, loci, reclocus,regions]
      loci = architechture(slimpath=slimpath, transitpath=transitpath, 
                           lociperchr=lociperchr, N=Nbase, cemeesnps=cemeesnps,
                           smoothmaps=smoothmaps,
                           fit.vars=unique(px$fitvar), betas=unique(px$beta), 
                           recpos = recpos, recchrom=recchrom, burn=burn, burntype = burntype, burn.fitvar=burn.fitvar)
      
      loci
      
      
    })
    
    
    gidinfo = do.call(rbind, lapply(gidloci, function(x){
      data.frame(gid=x$gid [1],nloci = nrow(x), reclocus = ifelse(gparam[i,]$recpos == "not_explicit","not_explicit", which(x$type=='rec')))
      }))
    
    
    tx = expand.grid(1:nrow(px), 1:nrow(gidinfo))
    px = cbind(px[tx[,1],],  gidinfo[tx[,2],])
    
    return(list(info = px, loci = do.call(rbind, gidloci)))
    
  }))
}




############################################################################################
#### II) Functions to simulate the evolution of the initial populations

### Run slim in the terminal nrun times and return output
### (!) slim file is the whole path (or need to setwd() to the file containing the .slim first)
### rrfile and genofile must be located in the same directory as slim
### nloci and N must match with genomatrix and rrmap files
runslimpolym = function(nrun, sfile, selfingrate,beta, nloci, N,rrfile,initialpopfile,slimpathfile, HighAndLowMod=0, reclocus, ngeneration = 100, outpath,outputMSgenerations=NULL){
  
  ## Toy parameter
  #nrun=1
  #scoef="\"0.1 0.2 0\""
  #selfingrate=0.2
  #nloci=15
  #N=1000
  #rec='wt'
  #rrfile="rrmap.txt"
  #genofile='genomatrix.txt'
  #slimfile="Data/Experimental_Evolution_R/EE_Rhom/SLiM/uniform/uniform_selection.slim"
  
  param_selfrate = paste0("-d self_rate=",selfingrate)
  param_nloci =  paste0("-d nloci=",nloci)
  param_reclocus =  paste0("-d reclocus=",reclocus-1) 
  param_N =  paste0("-d N=",N)
  param_rr = paste0("-d rrpath=\"'", rrfile, "'\"")
  param_s = paste0("-d spath=\"'", sfile, "'\"")
  param_beta =  paste0("-d beta=",beta)
  param_geno = paste0("-d initialpoppath=\"'", initialpopfile, "'\"")
  param_ngeneration = paste0("-d ngeneration=",ngeneration)
  param_mod =  paste0("-d highAndLowMod=",HighAndLowMod)
  param_outpath = paste0("-d outpath=\"'", outpath, "'\"")
  
  #Null = no output
  if(is.null(outputMSgenerations)){outputMSgenerations = -1} # This will ask SLiM to output at G = -1, which never happen.
  param_outGenerations = paste0("-d outGenerations=\"'", paste(outputMSgenerations, collapse = ";"), "'\"")
  
  slimcmd =paste("slim",
                 param_selfrate,
                 param_s,
                 param_beta,
                 param_nloci,
                 param_N,
                 param_rr,
                 param_geno,
                 param_mod,
                 param_reclocus,
                 param_ngeneration,
                 param_outpath,
                 param_outGenerations,
                 slimpathfile,
                 " ")
  
  
  slim_out = lapply(1:nrun, function(i){system( slimcmd,intern=TRUE)})
  
  slim_out
}





# Call runslimpolym for each parameters for a given initial pop 
# Return the parsed basic output (fitness, rec1 freq)
EVOLUTION_POLYM = function(gid, slimpath, slimfile,transitpath, nrun, selfingrate, nloci, reclocus, N,  
                fit.var,beta, HighAndLowMod,ngeneration=100, outputMSgenerations=NULL, burn=F){
  
  #gid: the gid of the target genomatrix/initialpop
  #slimpath: Where is the .slim file containing simulation script
  #slimfile: the slimfile for evolution
  #transitpath: the transitory file directory
  #nrun: number of simulation run per initial pop
  #selfingrate: selfing rate
  #nloci: number of loci
  #reclocus: position of the rec-1 locus (which locus; fisrt, second,.. ? )
  #N: population size during evolution
  #fit.var: fitness variance
  # beta:  a parameter that define the curvature of the fitness function during evolution (0 = no epistasis; 1= stabilizing selection)
  #HighAndLowMod:If 0: rec-1 WT vs Mut; If 1: WT vs Low rr modifier (25cM; codominant; WT rr landscape)
  #ngeneration: number of generation
  #outputMSgenerations:Generation at which we want to output genotype of the populations to extract metrics (fitness, rec-1 freq already given in a basic output)
  # burn: T or F
  
  
  ### Run n times the slim code for this specific genomatrix and parameters for each of rec-1 rrmap
  
  sfiles = list.files(paste0(transitpath, gid))
  sfiles = sfiles[grepl(paste0("fitvar",fit.var, '_'), sfiles) & grepl(paste0("beta",beta,'_'), sfiles)]
  
  if(burn){ sfiles = sfiles[grepl('burn', sfiles)]}else{sfiles = sfiles[!grepl('burn', sfiles)]}
  
  if(length(sfiles)>1 & sum(grepl('burn', sfiles))==0) stop("Duplicated selection files")
  if(length(sfiles)==0) stop("No selection files found")
  
  
  OUT = do.call(rbind, lapply(sfiles, function(sfile){
    
    if(burn){
      
      burntype = strsplit(sfile, "burn")[[1]][2]
      burntype = strsplit(burntype, "_")[[1]][1]
      
    }else{burntype = "noburn"}
    
    
    outputs = runslimpolym(nrun=nrun, sfile=paste0(gid,'/',sfile), beta=beta,
                           selfingrate=selfingrate, nloci=nloci, N=N,
                           rrfile= paste0(gid,"/rrmap_",gid,".txt"),
                           initialpopfile = paste0(gid,"/initialPop_", ifelse(burn, paste0("burnin_", burntype, '_') , NULL) ,gid,".txt"),
                           slimpathfile=paste0(slimpath,slimfile),
                           HighAndLowMod=HighAndLowMod, ngeneration=ngeneration, outpath=paste0(gid, "/"),
                           reclocus=reclocus,
                           outputMSgenerations=outputMSgenerations)
    
    
    
    outputs = do.call(rbind, lapply(1:length(outputs), function(i){ 
      x=outputs[[i]]
      x=parse_slim(x, full.numeric = T)
      x$run = i
      x
    }))
    
    outputs$burn_in = 'after_burn'
    outputs$burntype = burntype
    outputs$nloci = nloci -1 
    outputs$N = N
    outputs$fit.var = fit.var
    outputs$beta = beta
    outputs$selfing = selfingrate
    outputs$modifier = ifelse(HighAndLowMod==1, 'lowrecombiningallele', 'rec1')
    
    outputs
    
    
  }))
  
  
  
  OUT$gid = gid
  
  OUT
}





### Run slim ... in the terminal nrun times and return output
### (!) slim file is the whole path (or need to setwd() to the file containing the .slim first)
### rrfile and genofile must be located in the same directory as slim
### nloci and N must match with genomatrix and rrmap files

runslimmono = function(nrun, sfile, selfingrate, beta,nloci, N, rec,rrfile,ncrossover=1,initialpopfile,slimpathfile, ngeneration = 100, outpath,outputMSgenerations=NULL){
  
  ## Toy parameter
  #nrun=1
  #scoef="\"0.1 0.2 0\""
  #selfingrate=0.2
  #nloci=15
  #N=1000
  #rec='wt'
  #rrfile="rrmap.txt"
  #genofile='genomatrix.txt'
  #slimfile="Data/Experimental_Evolution_R/EE_Rhom/SLiM/uniform/uniform_selection.slim"
  
  param_selfrate = paste0("-d self_rate=",selfingrate)
  param_nloci =  paste0("-d nloci=",nloci)
  param_N =  paste0("-d N=",N)
  rec = c(0,1,2)[which(c("wt", "mut", "low")==rec)]; param_rec = paste0("-d rec=",rec) # WARNING: 0 = wt, 1 = mut
  param_nCO = paste0("-d nCO=",ncrossover)
  param_rr = paste0("-d rrpath=\"'", rrfile, "'\"")
  param_s = paste0("-d spath=\"'", sfile, "'\"")
  param_beta =  paste0("-d beta=",beta)
  param_geno = paste0("-d initialpoppath=\"'", initialpopfile, "'\"")
  param_ngeneration = paste0("-d ngeneration=",ngeneration)
  param_outpath = paste0("-d outpath=\"'", outpath, "'\"")
  
  #Null = no output
  if(is.null(outputMSgenerations)){outputMSgenerations = -1} # This will ask SLiM to output at G = -1, which never happen.
  param_outGenerations = paste0("-d outGenerations=\"'", paste(outputMSgenerations, collapse = ";"), "'\"")
 
  
  
  slimcmd =paste("slim",
                 param_selfrate,
                 param_s,
                 param_beta,
                 param_nloci,
                 param_N,
                 param_rec,
                 param_rr,
                 param_nCO,
                 param_geno,
                 param_ngeneration,
                 param_outpath,
                 param_outGenerations,
                 slimpathfile,
                 " ")
  
  
  slim_out = lapply(1:nrun, function(i){system( slimcmd,intern=TRUE)})
  
  slim_out
}




# Call runslimpolym for each parameters for a given initial pop 
# Return the parsed basic output (fitness, rec1 freq)
EVOLUTION_MONO = function(gid, slimpath, slimfile,transitpath, nrun, selfingrate, nloci, N,  
                           fit.var,beta, recombination,ngeneration=100, outputMSgenerations=NULL, burn=F){
  
  
  #gid: the gid of the target genomatrix/initialpop
  #slimpath: Where is the .slim file containing simulation script
  #slimfile: the slimfile for evolution
  #transitpath: the transitory file directory
  #nrun: number of simulation run per initial pop
  #selfingrate: selfing rate
  #nloci: number of loci
  #reclocus: position of the rec-1 locus (which locus; fisrt, second,.. ? )
  #N: population size during evolution
  #fit.var: fitness variance
  # beta:  a parameter that define the curvature of the fitness function during evolution (0 = no epistasis; 1= stabilizing selection)
  #recombination: c("wt", "mut", "high", "low"), see below
  #ngeneration: number of generation
  #outputMSgenerations:Generation at which we want to output genotype of the populations to extract metrics (fitness, rec-1 freq already given in a basic output)
  # burn: T or F
  
  
  # The "type" of recombination c("wt", "mut", "high", "low"), corresponds to different rr landscapes and number of CO
  wrecombination = which(c('wt', 'mut', 'low', 'high') %in% recombination)
  rlandscapes = c('wt', 'mut', 'low', 'wt')[wrecombination]
  numbercrossover = c(1,1,0,3)[wrecombination]
  
  
  ### Run n times the slim code for this specific genomatrix and parameters for each of rec-1 rrmap
  
  sfiles = list.files(paste0(transitpath, gid))
  sfiles = sfiles[grepl(paste0("fitvar",fit.var, '_'), sfiles) & grepl(paste0("beta",beta,'_'), sfiles)]
  
  if(burn){ sfiles = sfiles[grepl('burn', sfiles)]}else{sfiles = sfiles[!grepl('burn', sfiles)]}
  
  if(length(sfiles)>1 & sum(grepl('burn', sfiles))==0) stop("Duplicated selection files")
  if(length(sfiles)==0) stop("No selection files found")
  
  
  OUT = do.call(rbind, lapply(sfiles, function(sfile){
    
    if(burn){
      
      burntype = strsplit(sfile, "burn")[[1]][2]
      burntype = strsplit(burntype, "_")[[1]][1]
      
    }else{burntype = "noburn"}
    
    
    b_outputs = do.call(rbind, lapply(1:length(recombination), function(n){
    
    outputs = runslimmono(nrun=nrun,
                          sfile=paste0(gid,'/',sfile),
                          beta = beta,
                          selfingrate=selfingrate,
                          nloci=nloci,
                          N=N,
                          rec=rlandscapes[n],
                          rrfile=paste0(gid,"/rrmap_",gid,".txt"),
                          ncrossover=numbercrossover[n],
                          initialpopfile=paste0(gid,"/initialPop_", ifelse(burn, paste0("burnin_", burntype, '_') , NULL) ,gid,".txt"),
                          slimpathfile=paste0(slimpath,slimfile),
                          ngeneration = ngeneration,
                          outpath=paste0(gid, "/"),
                          outputMSgenerations=outputMSgenerations)
    
    
    
    outputs = do.call(rbind, lapply(1:length(outputs), function(i){ 
      x=outputs[[i]]
      x=parse_slim(x, full.numeric = T)
      x$run = i
      x
    })) # each run
    
    outputs$burn_in = 'after_burn'
    outputs$burntype = burntype
    outputs$nloci = nloci 
    outputs$N = N
    outputs$fit.var = fit.var
    outputs$beta = beta
    outputs$selfing = selfingrate
    outputs$rec = recombination[n]
    
    outputs
    
    })) # for each recombinaion
    
    b_outputs
    
  }))
  
  
  
  OUT$gid = gid
  
  OUT
}

##########################################################################################################
#### III) Functions calculate metrics from the genotype data outputed during simulation by slim (MS format)


extractDiffLogWFromOutputMS=function(loci, transitpath, fitnessout, burn = F){
  gid = unique(fitnessout$gid)
  fit.var = unique(fitnessout$fit.var)
  beta = unique(fitnessout$beta)
  
  files = list.files(paste0(transitpath, gid))
  files = files[grepl("outputMS", files)] 
  runseedfiles = tstrsplit(files, "_")[[3]]
  runseedfiles = tstrsplit(runseedfiles, ".txt")[[1]]
  
  locix = subset(loci, loci$gid== unique(fitnessout$gid))
  reclocus = which(locix$type=="rec")
  
  if(burn){burntype = "burntype"}else{burntype = NULL}
  
  uparam = unique(fitnessout[,c("fit.var", "beta",burntype)])
  
  
  diffsel = do.call(cbind, lapply(1:nrow(uparam), function(i){
    
    if(burn){
      fitvar = uparam$fit.var[i]
      type = uparam$burntype[i]
      beta = uparam$beta[i]
      ss = as.matrix(fread(paste0(transitpath, gid, "/selection_burn", type,"_fitvar", fitvar, "_beta", beta,"_",gid,".txt")))
      ss = ss[,-reclocus]
      runseed = as.character(unique(subset(fitnessout, fit.var==fit.var & burntype == type & beta == beta)$runseed))
    }else{
      fitvar = uparam$fit.var[i]
      beta = uparam$beta[i]
      ss = as.matrix(fread(paste0(transitpath, gid, "/selection_fitvar", fitvar, "_beta", beta,"_",gid,".txt")))
      ss = ss[,-reclocus]
      runseed = as.character(unique(subset(fitnessout, fit.var==fit.var & beta == beta)$runseed))
    }
    
    
    
    do.call(cbind, lapply(files[runseedfiles %in% runseed], function(ff){
      
      nloci = nrow(locix)-length(reclocus)
      outputMS <- fread(paste0(transitpath, gid, "/", ff))
      
      relpos = as.character(outputMS[2,1])
      relpos = as.numeric(strsplit(relpos," ")[[1]][-1])
      relpos = (relpos * nloci)+1
      if(sum(duplicated(relpos))){warning("rec-1 allele and selected allele at the same position "); reclocus = reclocus+1}
      
      geno = outputMS[-1:-2,1]
      geno = sapply(as.vector(geno), function(x){as.numeric(unlist(strsplit(as.character(x),"")))})
      geno = t(matrix(geno, ncol = length(relpos), byrow = T))
     
      
      mx = match(round(1:(nloci+length(reclocus)), digits = 0), round(relpos, digits = 0))    
      geno = geno[mx,]
      geno[is.na(geno)]=0
      
      if(length(reclocus)>0){
        rec1 = geno[reclocus,]
        geno = geno[-reclocus,]
      }
      
      
      
      
      diff.sel.unscaled =  apply(geno[,rec1==1], 1, mean) - apply(geno[,rec1==0], 1, mean)
      diff.sel.unscaled = as.vector(diff.sel.unscaled)
      
      freq = apply(geno, 1, sum)/ncol(geno)
      diff.sel.unscaled[freq == 0 | freq == 1]=NA # NA for fitness loci fixed; all Nan when rec-1 is fixed
      
      x = data.frame(diff.sel.unscaled*ss[1,])
      names(x) = paste0(gid,'_', tstrsplit(tstrsplit(ff, ".txt")[[1]], "outputMS_")[[2]])
      x
    }))
    
    
  }))
  
  return(cbind(locix[-reclocus,c("chrom", "pos")], diffsel))
  
  
}













path = paste0(transitpath, gid)
filename = "outputMS_G10_1927470421569.txt"

extractOutputMS=function(filename, nloci, transitpath, gid){
  
  outputMS <- fread(paste0(transitpath, gid, "/", filename))
  
  relpos = as.character(outputMS[2,1])
  relpos = as.numeric(strsplit(relpos," ")[[1]][-1])
  relpos = (relpos * nloci)+1
  if(sum(duplicated(relpos))){warning("rec-1 allele and selected allele at the same position ")}
  
  geno = outputMS[-1:-2,1]
  geno = sapply(as.vector(geno), function(x){as.numeric(unlist(strsplit(as.character(x),"")))})
  geno = t(matrix(geno, ncol = length(relpos), byrow = T))
  
  
  mx = match(round(1:nloci, digits = 0), round(relpos, digits = 0))    
  geno = geno[mx,]
  geno[is.na(geno)]=0
  
  return(geno)
}


GetDiffLogW = function(geno, recallele, ss){
  
  diff.sel.unscaled =  apply(geno[,recallele==1], 1, mean) - apply(geno[,recallele==0], 1, mean)
  diff.sel.unscaled = as.vector(diff.sel.unscaled)
  
  freq = apply(geno, 1, sum)/ncol(geno)
  diff.sel.unscaled[freq == 0 | freq == 1]=NA # NA for fitness loci fixed; all Nan when rec-1 is fixed
  
  data.frame(diff.sel.unscaled*ss[1,])
}

# For monomophic experiment, measure the allele frequencies and LogDiff Fitness
stat.mono = function(loci, transitpath, fitnessout, burn = F){
  
  loci = subset(loci, loci$gid== unique(fitnessout$gid))
  gid = unique(fitnessout$gid)
  
  files = list.files(paste0(transitpath, gid))
  files = files[grepl("outputMS", files)] 
  
  fit.var = unique(fitnessout$fit.var)
  beta = unique(fitnessout$beta)
  
  if(burn){burntype = "burntype"}else{burntype = NULL}
  
  uparam = unique(fitnessout[,c("fit.var", "beta",burntype)])
  
  INFO = do.call(c, lapply(1:nrow(uparam), function(i){
    
    if(burn){
      fitvar = uparam$fit.var[i]
      type = uparam$burntype[i]
      beta = uparam$beta[i]
      ss = as.matrix(fread(paste0(transitpath, gid, "/selection_burn", type,"_fitvar", fitvar, "_beta", beta,"_",gid,".txt")))
      #if(!is.null(reclocus)){ss = ss[,-reclocus]}
      runseeds = unique(subset(fitnessout, fit.var==fit.var & burntype == type & beta == beta)[,c("rec", "runseed", "run")])
    }else{
      fitvar = uparam$fit.var[i]
      beta = uparam$beta[i]
      ss = as.matrix(fread(paste0(transitpath, gid, "/selection_fitvar", fitvar, "_beta", beta,"_",gid,".txt")))
      #if(!is.null(reclocus)){ss = ss[,-reclocus]}
      runseeds = unique(subset(fitnessout, fit.var==fit.var & beta == beta)[,c("rec", "runseed", "run")])
    }
    
    
    out = lapply(unique(runseeds$run), function(runx){
      
      runseed = subset(runseeds, run == runx)[,c("rec", "runseed")]
      runseed$runseed = as.character(runseed$runseed)
      
      
      tarfiles = files[runseedfiles %in% runseed$runseed]
      
      generations = tstrsplit(tarfiles, "_G")[[2]]
      generations = as.numeric(tstrsplit(generations, "_")[[1]])
      
      genos = lapply(tarfiles, function(ff){
        gx = extractOutputMS(filename=ff, nloci=nrow(loci), transitpath=transitpath, gid=gid)
        gx
      })
      
      
      afc = do.call(rbind, lapply(runseed$runseed, function(r){
        
        genostart = genos[[which(grepl(paste0("_G", min(g), "_"), tarfiles) & grepl(r, tarfiles))]]
        genoend = genos[[which(grepl(paste0("_G", max(g), "_"), tarfiles) & grepl(r, tarfiles))]]
        
        freqstart = apply(genostart, 1, sum)/ncol(genostart)
        freqend = apply(genoend, 1, sum)/ncol(genoend)
        
        af = data.frame(initial.freq=freqstart, end.freq = freqend )
        af = cbind(loci, af)
        
      }))
      
      afc$beta = beta
      afc$fit.var = fit.var
      afc$burntype = type
      afc$run = runx
      
      
      diffsel = do.call(rbind, lapply(sort(unique(generations[generations != 0])), function(g){
        
        genowt = genos[[which(grepl(paste0("_G", g, "_"), tarfiles) & grepl(runseed$runseed[runseed$rec=="wt"], tarfiles))]]
        genomut = genos[[which(grepl(paste0("_G", g, "_"), tarfiles) & grepl(runseed$runseed[runseed$rec=="mut"], tarfiles))]]
        
        geno = cbind(genowt, genomut)
        
        recallele = c(rep(0, ncol(genowt)), rep(1, ncol(genomut)))
        
        DiffLogW = setNames(GetDiffLogW(geno, recallele, ss), "diffsel")
        DiffLogW$chrom = loci$chrom
        DiffLogW$pos = loci$pos
        DiffLogW$generation = g
        DiffLogW
      }))
      
      
      diffsel=diffsel[!(is.na(diffsel$diffsel) & !is.nan(diffsel$diffsel)),]
      diffsel = aggregate(diffsel~chrom+pos,data=diffsel, mean)
      diffsel=diffsel[order(diffsel$chrom, diffsel$pos),]
      
      diffsel$beta = beta
      diffsel$fit.var = fit.var
      diffsel$burntype = type
      diffsel$run = runx
      
      #ggplot(diffsel, aes(pos, diffsel))+geom_point()+facet_wrap(~chrom)
      
      
      return(list(afc, diffsel))
      
    })
    
    return(out)
  }))
  
  
  INFO = list(afc = do.call(rbind, lapply(INFO, function(x){x[[1]]})), diffsel = do.call(rbind, lapply(INFO, function(x){x[[2]]})))
  INFO[[1]] = merge(INFO[[1]], unique(fitnessout[,c("runseed","run","burn_in","burntype","nloci","N","fit.var","beta","selfing","rec","gid" )]))
  INFO[[2]] = merge(INFO[[2]], unique(fitnessout[,c("run","burn_in","burntype","nloci","N","fit.var","beta","selfing","gid" )]))
  
  return(INFO)
  

}
  




stat.poly = function(loci, transitpath, fitnessout, burn = F){
  
  loci = subset(loci, loci$gid== unique(fitnessout$gid))
  gid = unique(fitnessout$gid)
  
  files = list.files(paste0(transitpath, gid))
  files = files[grepl("outputMS", files)] 
  
  fit.var = unique(fitnessout$fit.var)
  beta = unique(fitnessout$beta)
  
  if(burn){burntype = "burntype"}else{burntype = NULL}
  
  reclocus = which(loci$type=="rec")
  
  uparam = unique(fitnessout[,c("fit.var", "beta",burntype)])
  
  INFO = do.call(rbind, lapply(1:nrow(uparam), function(i){
    
    if(burn){
      fitvar = uparam$fit.var[i]
      type = uparam$burntype[i]
      beta = uparam$beta[i]
      ss = as.matrix(fread(paste0(transitpath, gid, "/selection_burn", type,"_fitvar", fitvar, "_beta", beta,"_",gid,".txt")))
      if(!is.null(reclocus)){ss = ss[,-reclocus]}
      runseeds = as.character(unique(subset(fitnessout, fit.var==fit.var & burntype == type & beta == beta)$runseed))
    }else{
      fitvar = uparam$fit.var[i]
      beta = uparam$beta[i]
      ss = as.matrix(fread(paste0(transitpath, gid, "/selection_fitvar", fitvar, "_beta", beta,"_",gid,".txt")))
      if(!is.null(reclocus)){ss = ss[,-reclocus]}
      runseeds = as.character(unique(subset(fitnessout, fit.var==fit.var & beta == beta)$runseed))
    }
    
    
    out = do.call(rbind, lapply(runseeds, function(runseed){
      
      tarfiles = files[runseedfiles == runseed]
      
      generations = tstrsplit(tarfiles, "_G")[[2]]
      generations = as.numeric(tstrsplit(generations, "_")[[1]])
      
      genos = lapply(tarfiles, function(ff){
        gx = extractOutputMS(filename=ff, nloci=nrow(loci), transitpath=transitpath, gid=gid)
        gx
      })
      
      
      diffsel = do.call(rbind, lapply(sort(unique(generations[generations != 0])), function(g){
        
        geno = genos[[which(grepl(paste0("_G", g, "_"), tarfiles))]]
        
        recallele = geno[reclocus, ]
        geno = geno[-reclocus, ]
        
        DiffLogW = setNames(GetDiffLogW(geno, recallele, ss), "diffsel")
        DiffLogW$chrom = loci[-reclocus,]$chrom
        DiffLogW$pos = loci[-reclocus,]$pos
        DiffLogW$generation = g
        DiffLogW
      }))
      
      
      diffsel=diffsel[!(is.na(diffsel$diffsel) & !is.nan(diffsel$diffsel)),]
      diffsel = aggregate(diffsel~chrom+pos,data=diffsel, mean)
      diffsel=diffsel[order(diffsel$chrom, diffsel$pos),]
      
      diffsel$beta = beta
      diffsel$fit.var = fit.var
      diffsel$burntype = type
      diffsel$runseed = runseed
      
      #ggplot(diffsel, aes(pos, diffsel))+geom_point()+facet_wrap(~chrom)
      
      
      return(diffsel)
      
    }))
    
    return(out)
  }))
  
  INFO = merge(INFO, unique(fitnessout[,c("runseed","run","burn_in","burntype","nloci","N","fit.var","beta","selfing","modifier","gid","recpos" )]))
  return(INFO)
  
  
}





#### To do add info on afc 
#### return(list(afc, diffsel))
#### bind all run
#### Change a bit for the polymorphic
#### Incorporate in the main function
#### Save on github




# Calculate the allele frequency change



# Bind the geno of wt and mut and calculate the LogDiffFitness


















