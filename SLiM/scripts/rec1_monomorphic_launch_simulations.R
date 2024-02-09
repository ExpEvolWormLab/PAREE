#########################################################################
#### IMPORT DATA AND FUNCTIONS ##########################################
#########################################################################

source("/Users/tomparee/Desktop/SLiM/scripts/rec1_EE_simulation_sources.R")

slimoutpath = paste0(generalpath, "monomorphic_outputs/rawdata/test/") # Where to save 
slimfile = "rec1_monomorphic_simulation.slim" # The main slim file (NOTE: there is a second annex slim file "genomatrixToSlim.slim")

################################################################################
##################### FUNCTIONS ################################################
################################################################################

### Create a data.frame containing the parameters of evolution
paramdf.mono = function(N=1000,fit.vars=0.01,nloci=100, selfing=0, beta=0){
  
  #nloci = Number of loci PER CHROMOSOME (!!!) (without counting the rec-1 locus)
  #N = Pop size during evolution
  #fit.vars = Genetic variance for fitness
  #selfing = Selfing rate
  #beta =  a parameter that define the curvature of the fitness function during evolution (0 = no epistasis; 1= stabilizing selection)
  #recombination 
  
  
  parameters = setNames(expand.grid(list(nloci, fit.vars, N, beta, selfing)),
                        c('lociperchr', 'fitvar','popsize','beta','selfrate'))
  
  
  
  parameters
  
}







### The function to launch simulations with rec-1 polymorphic 
### This function call function stored in rec1_polymorphic_simulation_sources.R
simulation_mono = function(ngmatrix,nrunpergmatrix, parameters, slimpath, slimfile, transitpath, cemeesnps, smoothmaps, ngeneration=100, 
                            Nbase=1000,burn.fitvar=0.15,recombination = c('wt', 'mut', 'low', 'high'), burn=T, burntype=c("stab", "drift"), outputMSgenerations = NULL, ncores=3, dominance = F){
  
  #Toy parameters:
  # ngmatrix=2
  # nrunpergmatrix = 2
  # parameters = paramdf.mono(N=1000, fit.vars=c(0.01,0.02), nloci=100, selfing=0.15, beta=0)
  # ngeneration = 100
  # Nbase = 1000
  # burn.fitvar=0.15
  # burn=T
  # burntype=c("stab")
  # outputMSgenerations = seq(10,100,10)
  #recombination = c("wt", "mut")
  # ncores=3
  
  #####################################################################################################################
  ### 0) Subdivide the parameters depending if if require to create a different initial population (ex: number of loci)
  ### or not (i.e. can differ during evolution of the same initial population; ex: selfing rate)
  #####################################################################################################################
  
  # gparam = df with parameters that require to create a new genotype matrix
  # i.e. the position of rec-1 and the number of loci
  # => ex: different selfing rate can be done for the same genotype matrix
  
  nloci_param = unique(parameters[, "lociperchr"])
  
  lapply(nloci_param, function(lociperchr){
    #for each number of loci
  
    # Note: we could directly feed the entirety parameters in prepare.initial.pop
    # => in this case, all initial pop would be created before the simulations
    # Here, we first subset so when there is a bug in the simu, we don't wait the creation of all populations to notice it
    paramx = parameters[parameters$lociperchr==lociperchr, ] 
    paramx$recpos = "not_explicit"
    
    ####################################################################################################################
    ### I) Prepare the all initial populations (=ngmatrix) for one set of genotype matrix parameter (lociperchr, recpos)
    ####################################################################################################################
    
    # Function below creates initial pop .slim files and annex .txt files + return some info:
    info = INITIAL_POPS(ngmatrix=ngmatrix, #number of initial population (~genotype matrix) for a set of gmatrix param (nlociperchrom, recpos)
                        parameters=paramx, #The dataframe containing the paramters and created by paramdf.polym function
                        slimpath=slimpath, # The path to the directory containing the slim scripts
                        transitpath=transitpath, # The path to a directory to store tampôrary files created during simu
                        cemeesnps=cemeesnps, # The snps of the cemee (chrom, pos, etc.)
                        smoothmaps=smoothmaps, # The rec-1 WT & mut recombination maps (smoothed)
                        burn=burn, # if T, burn in phase of 100 generations following the creation of populatopn
                        burntype = burntype, # c("stab", "drift") or "stab" or "drift"; the type of burn in: stabilizing selection (i.e. negative epistasis) and/or drift.
                        burn.fitvar=burn.fitvar, # The genetic variance for fitness during burn-in (stab)
                        recchrom=NULL, # The chromosome containing rec-1
                        Nbase = Nbase, # The pop size of the initial population (just matter to sampled the cemee genetic var and the burn-in; but for evolution N can be adjusted)
                        dominance = dominance) # If F, codominance (h=0.5). If T, h sampled from uniform distrib
    
    paramx = info[[1]] # combination of parameters and gid (gmatrix identifier)
    loci = info[[2]] # Loci and their positon
    
    
    ############################################################################################################
    ### II) For these each of these initial pop created above run simulations for each set of "evolution" param
    ### (fitvar, popsize, beta, selfrate, modifier) * nrun
    ############################################################################################################
    #=> i.e.  For each paramter combination, call SLIM
    mclapply(1:nrow(paramx),mc.cores = ncores, function(p){
      
      print(paramx$gid[p])
      
      
      ### Call SLiM for simulation
      ### Return basic outputs (fitness, rec1 mutant frequency)
      fitnessout = EVOLUTION_MONO(gid = paramx$gid[p], # The identifier of the initial pop
                                   slimpath=slimpath, # The path to the directory containing the slim scripts
                                   slimfile=slimfile, # The .slim file encoding the simulation
                                   transitpath=transitpath,
                                   nrun=nrunpergmatrix,# The number of simulation run per initial pop
                                   recombination = recombination, # The position of the rec-1 locus (in number of loci)
                                   fit.var = paramx$fitvar[p], # Genetic variance for fitness
                                   beta = paramx$beta[p], # a parameter that define the curvature of the fitness function (0 = no epistasis; 1= stabilizing selection)
                                   selfingrate=paramx$selfrate[p], # The selfing rate
                                   nloci=paramx$nloci[p], # The number of loci (nlociperchr * nchrom +1) 
                                   N=paramx$popsize[p], # Population size during evolution
                                   ngeneration = ngeneration, # Number of generation of evolution
                                   outputMSgenerations = outputMSgenerations, #Vector of generation at which a MS output in created (all genotypes in the pop; useful to subsequently calculate metrics)
                                   burn=burn) # T or F
      
      gid=fitnessout$gid[1]
      
      
      ### Write the .txt file with fitness
      write.table(fitnessout, file = paste0(slimoutpath, "fitness_",ifelse(paramx$modifier[p]==1, 'low', 'rec1'),
                                            "_pos", paramx$recpos[p],
                                            "_nloci", paramx$lociperchr[p],"_fitvar", paramx$fitvar[p],
                                            "_selfing", paramx$selfrate[p],
                                            "_N", paramx$popsize[p], '_',gid,".txt"), 
                  col.names = T, row.names = F)
      
      
      
      ############################################################################################################
      ### III) Calculate metrics from the genotype output (if they exist)
      ### You may want to change this part to calculate different metrics
      ############################################################################################################
      if(!is.null(outputMSgenerations)){
        stat = stat.mono(loci = loci, transitpath = transitpath, fitnessout = fitnessout, burn=burn)
        allele.freq = stat[[1]]
        diffsel = stat[[2]]
        
        #ggplot(diffsel, aes(pos, diffsel))+geom_point()+facet_grid(run~generation)
        #ggplot(aggregate(diffsel~run+generation,data=diffsel, sum), aes(generation,diffsel))+geom_point()+facet_wrap(~run)+geom_hline(yintercept = 0)
        #ggplot(aggregate(diffsel~run+pos,data=diffsel, mean), aes(pos,diffsel))+geom_point()+facet_wrap(~run)+geom_hline(yintercept = 0)+geom_smooth(method = loess, span = 0.3)
        #aggregate(diffsel~run+pos,data=diffsel, mean)
        
        write.table(diffsel, file = paste0(slimoutpath, "deltaLogFitness_",
                                           "_nloci", paramx$lociperchr[p],"_fitvar", paramx$fitvar[p],
                                           "_selfing", paramx$selfrate[p],
                                           "_N", paramx$popsize[p], '_',gid,".txt"), 
                    col.names = T, row.names = F)
        
        write.table(allele.freq, file = paste0(slimoutpath, "AlleleFreq_",
                                           "_nloci", paramx$lociperchr[p],"_fitvar", paramx$fitvar[p],
                                           "_selfing", paramx$selfrate[p],
                                           "_N", paramx$popsize[p], '_',gid,".txt"), 
                    col.names = T, row.names = F)
      }
      
    })
    
    
    # Delete all transitory files for this gparam[i,]
    lapply(unique(paramx$gid), function(gid){
      system(paste0("rm -r ",transitpath,gid))
    })
    
    
  })
  
}



####################################################################################
######## LAUNCH SIMULATIONS ########################################################
####################################################################################


# 1) Set the parameters
parameters = paramdf.polym(nloci=c(250), # Number of loci PER CHROMOSOME (!!!) (without counting the rec-1 locus)
                           N=c(1000), # Pop size during evolution
                           fit.vars = c(0.02), # Genetic variance for fitness
                           recpos =c(717890), # Position of rec-1 in bp
                           selfing=c(0.15), # Selfing rate
                           beta = 0, # a parameter that define the curvature of the fitness function during evolution (0 = no epistasis; 1= stabilizing selection)
                           modifier=0) # # If 0: rec-1 WT vs Mut; If 1: WT vs Low rr modifier (25cM; codominant; WT rr landscape)

# 2) Launch simulations
simulation_polym(ngmatrix=2, # The number of gmatrix per set of parameters
                 nrunpergmatrix=2, # The number of replicate simulation per initial population (~gmatrix)
                 parameters=parameters, # The df containing parameters combinations created by paramdf.polym
                 slimpath=slimpath, # path to .slim scripts
                 slimfile=slimfile, # The .slim file of evolution 
                 transitpath= transitpath, # Some directory where the transitory file can be saved
                 cemeesnps=cemeesnps, #The snps of the cemee (chrom, pos, etc.)
                 smoothmaps= smoothmaps, # The rec-1 WT & mut recombination maps (smoothed)
                 ngeneration=100, # Number of generation of evolution
                 Nbase=1000, # The size of the initial pop (Sampling CeMEE variation + burn-in) !!!: not the population size during evolution (contained in parameters df)
                 burn.fitvar=0.15, # Fitness variance during burn-in with negative epistasis (stabilizing selection)
                 recchrom=c("I"), # chromosome containing rec-1
                 burn=T, # T if require a burn-in
                 burntype=c("stab"), # c("stab", "drift") or "stab" or "drift"; the type of burn in: stabilizing selection (i.e. negative epistasis) and/or drift.
                 outputMSgenerations = seq(10,100,10), # Generation at which we want to output genotype of the populations to extract metrics (fitness, rec-1 freq already given in a basic output)
                 ncores=3, # number of core for parallelisaion (i.e. mclapply(mc.cores = ncores))
                 dominance = dominance) # If F, codominance (h=0.5). If T, h sampled from uniform distrib)





