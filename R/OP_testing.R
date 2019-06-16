#' Generate Open Pollinated Families
#'
#' This function simulates the creation of an Open Pollinated design and can make selections of which Parents perform the best
#' @param map.info Object returned from create_genetic_map()
#' @param parent.info Object returned from create_parents()
#' @param parent.phenos Object returned from sim_parents_phenos()
#' @param parents.TGV Object returned from calc_parents_TGV()
#' @param cross.prog The number of crosses to be made between each founder parent
#' @param num.select The number of parents to select from the OP cross design
#' @param dom.coeff The dominance cofficient used in estimating genetic values
#' @param A The value assigned to the major allele
#' @param a The value assigned to the minor allele
#' @param h2 Heritability to be used in estimating phenotypes
#' @param run.parallel logical. Set TRUE to run in parallel.
#' @param n.cores The number of cores which will be used if ran in parallel
#' @keywords Generate OP families
#' @export
#' @examples
#' op.families <- OP_testing( map.info = the.map, parent.info = the.parents, parent.phenos = parent.PHENOS,
#' parents.TGV = parent.TGV, cross.prog = 1, num.select = 64, dom.coeff = 1, A = 1, a = -100, h2 = .3, 
#' run.parallel = T,n.cores = 3)


#Extract OP genotypic values####
OP_testing <- function(map.info,parent.info,parent.phenos,parents.TGV,num.select,cross.prog = 1,dom.coeff,A,a,
                       h2, E.var = NULL, run.parallel=F, n.cores= NULL ) {
  library(abind)
  num.pars <- as.numeric(length(colnames(parent.info$genos.3d)))                         
  parent1 <-  rep(1:num.pars,each=(num.pars-1)); parent2 <- vector()
  for(i in 1:num.pars){
    all.parents <- 1:num.pars
    parent2 <- c(parent2,all.parents[-i])}
  progeny <- rep(cross.prog,length(parent2))
  OP.crosses <- data.frame(par1=parent1,par2=parent2,num.prog=progeny)
  
  if(run.parallel){
    library(parallel)
    g.val <- mclapply(1:nrow(OP.crosses),mc.cores = n.cores,FUN = function(x,crosses=OP.crosses){
      
      y= crosses[x,]
      chromo.last.loci.index <- map.info$last.locus.per.chrom   # The last number of each loci for all chromsomes
      num.chromos <- length(chromo.last.loci.index)
      QTLSNPs      <- map.info$QTLSNP.loci # Vector of loci which are snpqtl
      m <- map.info$genetic.map
      
      par1<-match(y[1],colnames(parent.info$genos.3d)) # assigns par1 to be the first parent in crossdesign matrix
      par2<-match(y[2],colnames(parent.info$genos.3d)) # assigns par2 to be the second parent in the crossdesign matrix
      cross.prog<-as.numeric(y[3]) #assigns number of progeny to be the third column for cross "X"
      
      # Create empty matrix to hold gametes
      # dimensions are (total # of loci) x  (# of cross progeny)
      # rownames are the loci names
      gametes1<-matrix(rep(NA,map.info$total.loci.num*cross.prog),nrow=map.info$total.loci.num,ncol=cross.prog,dimnames=list(map.info$map$loci,NULL))
      gametes2<-matrix(rep(NA,map.info$total.loci.num*cross.prog),nrow=map.info$total.loci.num,ncol=cross.prog,dimnames=list(map.info$map$loci,NULL))
      
      par1.alleles <- (parent.info$genos.3d[,par1,])
      par2.alleles <- (parent.info$genos.3d[,par2,])
      chr.ind.r <- vector("list")
      for(each in 1:cross.prog){
        first.pos <- 1
        ch.r <- vector("list")
        for(i in 1:length(chromo.last.loci.index)){
          last.pos <- chromo.last.loci.index[i]
          t <- sample(seq(0,1,.0001),chromo.last.loci.index[1],replace = F)
          l <- which(t < m[first.pos:last.pos,7])
          while(length(l) < 1){t <- sample(seq(0,1,.0001),chromo.last.loci.index[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])}
          ch.r[[i]] <- seq(first.pos,last.pos,1)[l]
          first.pos <- last.pos+ 1
        }
        chr.ind.r[[each]] <- unlist(ch.r)}
      
      for(i in 1:cross.prog) {
        allele <- sample(1:2,1)
        end <- chromo.last.loci.index[num.chromos]
        z <-1
        recombination.spots <- chr.ind.r[[i]]
        for (each in 1:length(recombination.spots)) {
          if (each < length(recombination.spots)){
            gametes1[z:recombination.spots[each],i] <- par1.alleles[z:recombination.spots[each],allele]
          } else { gametes1[z:end,i] <- par1.alleles[z:end,allele] }
          if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
          z <- recombination.spots[each] + 1
        }}
      
      chr.ind.r <- vector("list")
      for(each in 1:cross.prog){
        first.pos <- 1
        ch.r <- vector("list")
        for(i in 1:length(chromo.last.loci.index)){
          last.pos <- chromo.last.loci.index[i]
          t <- sample(seq(0,1,.0001),chromo.last.loci.index[1],replace = F)
          l <- which(t < m[first.pos:last.pos,7])
          while(length(l) < 1){t <- sample(seq(0,1,.0001),chromo.last.loci.index[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])}
          ch.r[[i]] <- seq(first.pos,last.pos,1)[l]
          first.pos <- last.pos+ 1
        }
        chr.ind.r[[each]] <- unlist(ch.r)}
      
      for(i in 1:cross.prog) {
        allele <- sample(1:2,1)
        end <- chromo.last.loci.index[num.chromos]
        z <-1
        recombination.spots <- chr.ind.r[[i]]
        for (each in 1:length(recombination.spots)) {
          if (each < length(recombination.spots)){
            gametes2[z:recombination.spots[each],i] <- par2.alleles[z:recombination.spots[each],allele]
          } else { gametes2[z:end,i] <- par2.alleles[z:end,allele] }
          if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
          z <- recombination.spots[each] + 1
        }}
      array.out <- abind(gametes1,gametes2,along=3)
      
      locus.names <- map.info$genetic.map$loci # The locus names pulled from the mab object
      QTLSNPs <- map.info$QTLSNP.loci     # vector of the loci which are snpqtl
      QTLSNP.num <- array.out[QTLSNPs,,] # genotypes of both alleles pulled from the current generation
      num.QTL <- length(map.info$rQTL.loci) # the number of additive qtl
      
      # Create 2 matrices: One to hold snpqtl values for all parents and the other to hold marker marker values for blup analysis
      num.SNPQTL <- map.info$total.SNPQTL.num # the number of loci which are snpqtl
      QTLSNP.values <-matrix(NA,nrow=num.SNPQTL,ncol=1) # matrix to hold snpqtl values
      
      QTLSNPaa <- which(QTLSNP.num[,1]=="a" & QTLSNP.num[,2]=="a")
      QTLSNPcc <- which(QTLSNP.num[,1]=="c" & QTLSNP.num[,2]=="c")
      QTLSNPac <- which(QTLSNP.num[,1]=="a" & QTLSNP.num[,2]=="c")
      QTLSNPca <- which(QTLSNP.num[,1]=="c"  & QTLSNP.num[,2]=="a")
      if (dom.coeff==0){
        QTLSNP.values[QTLSNPaa,1] <- A*2
        QTLSNP.values[QTLSNPcc,1] <- a*2
        QTLSNP.values[QTLSNPac,1] <- (A+a+dom.coeff)
        QTLSNP.values[QTLSNPca,1] <- (A+a+dom.coeff) } else {
          QTLSNP.values[QTLSNPaa,1] <- A
          QTLSNP.values[QTLSNPcc,1] <- a
          QTLSNP.values[QTLSNPac,1] <- (A*dom.coeff)
          QTLSNP.values[QTLSNPca,1] <- (A*dom.coeff)}
      
      par.QTL.allele1 <- as.numeric(array.out[map.info$rQTL.loci,,1])
      par.QTL.allele2 <- as.numeric(array.out[map.info$rQTL.loci,,2])
      
      
      # Genetic values of progeny
      #genetic.values <- sum(QTLSNP.values[,1]) + finalqtl
      genetic.values <- sum(QTLSNP.values[,1]) + sum(par.QTL.allele1) + sum(par.QTL.allele2)
      genetic.values
    })
  } else {
    g.val <- lapply(1:nrow(OP.crosses),FUN = function(x,crosses=OP.crosses){
      num.chromos <- length(chromo.last.loci.index)
      y= crosses[x,]
      chromo.last.loci.index <- map.info$last.locus.per.chrom   # The last number of each loci for all chromsomes
      QTLSNPs      <- map.info$QTLSNP.loci # Vector of loci which are snpqtl
      m <- map.info$genetic.map
      
      par1<-match(y[1],colnames(parent.info$genos.3d)) # assigns par1 to be the first parent in crossdesign matrix
      par2<-match(y[2],colnames(parent.info$genos.3d)) # assigns par2 to be the second parent in the crossdesign matrix
      cross.prog<-as.numeric(y[3]) #assigns number of progeny to be the third column for cross "X"
      
      # Create empty matrix to hold gametes
      # dimensions are (total # of loci) x  (# of cross progeny)
      # rownames are the loci names
      gametes1<-matrix(rep(NA,map.info$total.loci.num*cross.prog),nrow=map.info$total.loci.num,ncol=cross.prog,dimnames=list(map.info$map$loci,NULL))
      gametes2<-matrix(rep(NA,map.info$total.loci.num*cross.prog),nrow=map.info$total.loci.num,ncol=cross.prog,dimnames=list(map.info$map$loci,NULL))
      
      par1.alleles <- (parent.info$genos.3d[,par1,])
      par2.alleles <- (parent.info$genos.3d[,par2,])
      chr.ind.r <- vector("list")
      for(each in 1:cross.prog){
        first.pos <- 1
        ch.r <- vector("list")
        for(i in 1:length(chromo.last.loci.index)){
          last.pos <- chromo.last.loci.index[i]
          t <- sample(seq(0,1,.0001),chromo.last.loci.index[1],replace = F)
          l <- which(t < m[first.pos:last.pos,7])
          while(length(l) < 1){t <- sample(seq(0,1,.0001),chromo.last.loci.index[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])}
          ch.r[[i]] <- seq(first.pos,last.pos,1)[l]
          first.pos <- last.pos+ 1
        }
        chr.ind.r[[each]] <- unlist(ch.r)}
      
      for(i in 1:cross.prog) {
        allele <- sample(1:2,1)
        end <- chromo.last.loci.index[num.chromos]
        z <-1
        recombination.spots <- chr.ind.r[[i]]
        for (each in 1:length(recombination.spots)) {
          if (each < length(recombination.spots)){
            gametes1[z:recombination.spots[each],i] <- par1.alleles[z:recombination.spots[each],allele]
          } else { gametes1[z:end,i] <- par1.alleles[z:end,allele] }
          if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
          z <- recombination.spots[each] + 1
        }}
      
      chr.ind.r <- vector("list")
      for(each in 1:cross.prog){
        first.pos <- 1
        ch.r <- vector("list")
        for(i in 1:length(chromo.last.loci.index)){
          last.pos <- chromo.last.loci.index[i]
          t <- sample(seq(0,1,.0001),chromo.last.loci.index[1],replace = F)
          l <- which(t < m[first.pos:last.pos,7])
          while(length(l) < 1){t <- sample(seq(0,1,.0001),chromo.last.loci.index[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])}
          
          ch.r[[i]] <- seq(first.pos,last.pos,1)[l]
          first.pos <- last.pos+ 1
        }
        chr.ind.r[[each]] <- unlist(ch.r)}
      
      for(i in 1:cross.prog) {
        allele <- sample(1:2,1)
        end <- chromo.last.loci.index[num.chromos]
        z <-1
        recombination.spots <- chr.ind.r[[i]]
        for (each in 1:length(recombination.spots)) {
          if (each < length(recombination.spots)){
            gametes2[z:recombination.spots[each],i] <- par2.alleles[z:recombination.spots[each],allele]
          } else { gametes2[z:end,i] <- par2.alleles[z:end,allele] }
          if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
          z <- recombination.spots[each] + 1
        }}
      array.out <- abind(gametes1,gametes2,along=3)
      
      locus.names <- map.info$genetic.map$loci # The locus names pulled from the mab object
      QTLSNPs <- map.info$QTLSNP.loci     # vector of the loci which are snpqtl
      QTLSNP.num <- array.out[QTLSNPs,,] # genotypes of both alleles pulled from the current generation
      num.QTL <- length(map.info$rQTL.loci) # the number of additive qtl
      
      # Create 2 matrices: One to hold snpqtl values for all parents and the other to hold marker marker values for blup analysis
      num.SNPQTL <- map.info$total.SNPQTL.num # the number of loci which are snpqtl
      QTLSNP.values <-matrix(NA,nrow=num.SNPQTL,ncol=1) # matrix to hold snpqtl values
      
      QTLSNPaa <- which(QTLSNP.num[,1]=="a" & QTLSNP.num[,2]=="a")
      QTLSNPcc <- which(QTLSNP.num[,1]=="c" & QTLSNP.num[,2]=="c")
      QTLSNPac <- which(QTLSNP.num[,1]=="a" & QTLSNP.num[,2]=="c")
      QTLSNPca <- which(QTLSNP.num[,1]=="c"  & QTLSNP.num[,2]=="a")
      if (dom.coeff==0){
        QTLSNP.values[QTLSNPaa,1] <- A*2
        QTLSNP.values[QTLSNPcc,1] <- a*2
        QTLSNP.values[QTLSNPac,1] <- (A+a+dom.coeff)
        QTLSNP.values[QTLSNPca,1] <- (A+a+dom.coeff) } else {
          QTLSNP.values[QTLSNPaa,1] <- A
          QTLSNP.values[QTLSNPcc,1] <- a
          QTLSNP.values[QTLSNPac,1] <- (A*dom.coeff)
          QTLSNP.values[QTLSNPca,1] <- (A*dom.coeff)}
      
      par.QTL.allele1 <- as.numeric(array.out[map.info$rQTL.loci,,1])
      par.QTL.allele2 <- as.numeric(array.out[map.info$rQTL.loci,,2])
      
      # Genetic values of progeny
      genetic.values <- sum(QTLSNP.values[,1]) + sum(par.QTL.allele1) + sum(par.QTL.allele2)
      genetic.values
    })
  }
  
  g.val <- unlist(g.val)
  total.indiv <- length(g.val)
  if(length(E.var) > 0){ phenos <- g.val + rnorm(total.indiv,mean = 0,sd = E.var) } else {
    phenos <- g.val + rnorm(total.indiv,mean = 0,sd = sqrt(var(g.val)/h2))}
  
  mean.phenos <- vector()
  first <- 1
  last <- num.pars-1
  for(i in 1:num.pars){
    mean.phenos <- c(mean.phenos,mean(phenos[first:last]))
    first <- last+1
    last <- first + num.pars-2} 
  names(parent.phenos$phenos) <- colnames(parent.info$genos.3d)
  names(mean.phenos) <- colnames(parent.info$genos.3d)
  pars <- sort(mean.phenos,decreasing = T)[1:num.select]
  delt.alleles <- parent.info$delt.allele[as.numeric(names(pars))]
  genetic.values <- parent.phenos$genetic.values[as.numeric(names(pars))]
  phenos <- parent.phenos$phenos[as.numeric(names(pars))]
  genos.3d <- parent.info$genos.3d[,as.numeric(names(pars)),]
  marker.matrix <- parents.TGV$markers.matrix[as.numeric(names(pars)),]
  
  out <- list(pars=pars,delt.alleles=delt.alleles,genetic.values=genetic.values,phenos=phenos, genos.3d=genos.3d,
              marker.matrix=marker.matrix)
  return(out)
}