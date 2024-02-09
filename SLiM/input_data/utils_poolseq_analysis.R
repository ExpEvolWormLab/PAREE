

assign.info.melt = function(melted.freq.matrix, haplo = F){
  require(data.table)
  names(melted.freq.matrix) = c('marker', 'sample', 'freq')
  melted.freq.matrix$chrom = tstrsplit(melted.freq.matrix$marker, '_')[[1]]
  melted.freq.matrix$pos = as.numeric(tstrsplit(melted.freq.matrix$marker, '_')[[2]])
  
  if(grepl("_g", melted.freq.matrix$sample[1])){
    melted.freq.matrix$generation = as.numeric(tstrsplit(melted.freq.matrix$sample, '_g')[[2]])
    melted.freq.matrix$pop = tstrsplit(melted.freq.matrix$sample, '_g')[[1]]
  }else{names(melted.freq.matrix)[2]='pop'}
  
  melted.freq.matrix = assign.infopop(melted.freq.matrix)
  
  if(haplo == T){names(melted.freq.matrix)[1] = 'haplo'
  melted.freq.matrix = melted.freq.matrix[,-which(names(melted.freq.matrix)=="pos")] 
  }
  
  return(melted.freq.matrix)
}

#
arraytomatrix = function(arr,dimtomerge=NULL){
  require(abind)
  if(is.null(dimtomerge)) dimtomerge= 2:length(dim(arr))
  dims=dim(arr)[dimtomerge]
  dims=lapply(dims,function(x){1:x})
  
  dcomb = expand.grid(dims[1:(length(dims)-1)])
  tardim = dimtomerge[1:(length(dimtomerge)-1)]
  dn = dimnames(arr)
  
  mx = NULL
  for(i in 1:nrow(dcomb)){
    tar = dcomb[i,]
    tarl=list()
    for(n in 1:length(tar)){tarl[n]=dn[[tardim[n]]][unlist(tar[n])]}
    x=asub(arr,tarl,tardim)
    colnames(x) = paste0(colnames(x),"_",paste(unlist(tarl), collapse='_'))
    mx=cbind(mx,as.matrix(x))
  }
  return(mx)
}

# Asign populations informations
assign.infopop = function(data, col.pop.name='pop'){
  x = which(colnames(data)==col.pop.name)
  data$rec[grepl('LR',data[,x])]='mut'
  data$rec[grepl('HR',data[,x])]='wt'
  data$env[grepl('S',data[,x])]='salt'
  data$env[!grepl('S',data[,x])]='ngm'
  return(data)
}


interpolateGeneticPos <- function(x, gmap, xChrom='chrom', xPos='pos', gChrom='chrom', gPos = 'pos', gGpos='genetic', F2=FALSE){
  
  rename <- function(df, namein, nameout){
    stopifnot(namein %in% names(df))
    names(df)[names(df)==namein] <- nameout
    df
  }
  
  if(xChrom!='chrom') x = rename(x, xChrom, 'chrom')
  if(xPos!='pos') x = rename(x, xPos, 'pos')
  if(gChrom!='chrom') gmap = rename(gmap, gChrom, 'chrom')
  if(gPos!='pos') gmap = rename(gmap, gPos, 'pos')
  if(gGpos!='genetic') gmap = rename(gmap, gGpos, 'genetic')
  stopifnot(all(sum(c('chrom', 'pos') %in% intersect(names(x), names(gmap))), 'genetic' %in% names(gmap)))
  
  x <- x[order(x$chrom, x$pos),]
  g <- split(gmap, gmap$chrom)
  x$genetic <- unlist(lapply(split(x, x$chrom), function(i) approx(g[[i$chrom[1]]]$pos, g[[i$chrom[1]]]$genetic, xout = i$pos, rule=2)$y))
  if(F2) x <- do.call(rbind, lapply(split(x, x$chrom), function(i) {i$genetic = i$genetic * 50/max(i$genetic); i}))
  
  if(xChrom!='chrom') x = rename(x, 'chrom', xChrom)
  if(xPos!='pos') x = rename(x, 'pos', xPos)
  x
}



mutant.gmap = function(x){
  do.genetic.cM = function(target){
    minpos=min(target$POS[target$domain %in% c('arm', 'center')],na.rm=T)
    maxpos=max(target$POS[target$domain %in% c('arm', 'center')],na.rm=T)
    rr=50/(maxpos-minpos)
    mut.genetic=rep(NA, length(target$POS))
    mut.genetic[target$domain %in% c('arm', 'center')] = rr*(target$POS[target$domain %in% c('arm', 'center')]-minpos)
    mut.genetic[target$POS<minpos] = 0
    mut.genetic[target$POS>maxpos] = 50
    return(mut.genetic)
  }
  
  return(unlist(lapply(split(snps, snps$CHROM), function(i) do.genetic.cM(i))))
}




#Just a function which do all combination of a vector c(1,2,3) => 1-1,1-2,1-3, 2-2,2-3, 3-3

vector.combination <- function(t, with.itself=T, param.name=c('r1','r2')){
  if(with.itself){
    r1=NULL
    r2 = NULL
    for(n in 1:length(t)){
      r1 = c(r1, rep(t[n],((length(t)+1)-(1:length(t)))[n]))
      r2=c(r2, t[n:length(t)])
    }
  }else{
    r1= (length(t))-(1:length(t))
    r1=NULL
    r2 = NULL
    for(n in 1:length(t)){
      r1 = c(r1, rep(t[n],(length(t)-(1:length(t)))[n]))
      r2=c(r2, t[(n+1):length(t)])
    }
    r2=r2[1:length(r1)]
  }
  
  r = data.frame(r1,r2)
  colnames(r)=param.name
  
  return(r)
}




do.bins = function(snps, binlength=100000, bp=T){
  binsnp=function(x,nsnp){
    #nsnp =300
    x$bin=NA
    for(n in 1:(floor(nrow(x)/nsnp))){
      x$bin[(((n-1)*nsnp)+1):(n*nsnp)]=n
    }
    x$bin[which(is.na(x$bin))]=max(x$bin, na.rm=T)
    return(x)
  }
  
  binbp=function(x,nbp){
    require(birk)
    size=max(x$POS, na.rm=T)
    x$bin=NA
    for(n in 1:(floor(size/nbp))){
      start=(n-1)*nbp
      end=n*nbp
      x$bin[x$POS>=start & x$POS<end]=n
    }
    x$bin[which(is.na(x$bin))]=max(x$bin, na.rm=T)
    return(x)
  }
  
  
  if(bp){x <- do.call(rbind, lapply(split(snps, snps$chrom), function(i) {i=binbp(i,binlength); i}))
  }else{
    x <- do.call(rbind, lapply(split(snps, snps$chrom), function(i) {i=binsnp(i,binlength); i}))
  }
  
  x
  
}


# Change chromosomes id Roman <-> num

numtorom = function(vector_chromid){
  vector_chromid = as.character(as.roman(vector_chromid))
  vector_chromid[vector_chromid=="VI"]="X"
  vector_chromid = factor(vector_chromid, levels = c("I","II","III","IV","V","X"))
  return(vector_chromid)
}

romtonum = function(vector_chromid){
  vector_chromid[vector_chromid=="I"]=1
  vector_chromid[vector_chromid=="II"]=2
  vector_chromid[vector_chromid=="III"]=3
  vector_chromid[vector_chromid=="IV"]=4
  vector_chromid[vector_chromid=="V"]=5
  vector_chromid[vector_chromid=="X"]=6
  return(vector_chromid)
}






