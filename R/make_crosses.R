#' Make crosses
#'
#' This function makes crosses among parents in a format specified by the cross design file
#' @param parent.info Object returned from create_parents()
#' @param map.info Object returned from create_genetic_map()
#' @param cross.design Object returned from create_cross_design()
#' @param run.parallel logical. Set TRUE to run function in parallel
#' @param num.cores The number of cores to use if running in parallel
#' @keywords create crosses among parents
#' @export
#' @examples
#' progeny1 <- make_crosses(parent.info = op.families,map.info = the.map,cross.design = cross.file,run.parallel = T,num.cores = 3)

####Create Make Crosses####
make_crosses <- function(parent.info,map.info,cross.design, run.parallel = F, num.cores = NULL){
  library(parallel); library(abind)
  cross.design <- cross.design$cross.design
  num.crosses <- as.numeric(nrow(cross.design))
  num.chromos <- length(map.info$last.locus.per.chrom)

  if(run.parallel){
    gametes1 <- mclapply(1:num.crosses,function(x){
      chromo.loci.index <- map.info$last.locus.per.chrom   # The last number of each loci for all chromsomes
      QTLSNPs      <- map.info$QTLSNP.loci # Vector of loci which are snpqtl
      m <- map.info$genetic.map

      par1<-match(cross.design[x,1],colnames(parent.info$genos.3d)) # assigns par1 to be the first parent in cross.design matrix
      par2<-match(cross.design[x,2],colnames(parent.info$genos.3d)) # assigns par2 to be the second parent in the cross.design matrix
      cross.prog<-as.numeric(cross.design[x,3]) #assigns number of progeny to be the third column for cross "X"

      # Create empty matrix to hold gametes
      # dimensions are (total # of loci) x  (# of cross progeny)
      # rownames are the loci names
      gametes1<-matrix(rep(NA,map.info$total.loci.num*cross.prog),nrow=map.info$total.loci.num,ncol=cross.prog,dimnames=list(map.info$map$loci,NULL))

      par1.alleles <- (parent.info$genos.3d[,par1,])
      par2.alleles <- (parent.info$genos.3d[,par2,])
      chr.ind.r <- vector("list")
      for(each in 1:cross.prog){
        first.pos <- 1
        ch.r <- vector("list")
        for(i in 1:length(chromo.loci.index)){
          last.pos <- chromo.loci.index[i]
          t <- sample(seq(0,1,.0001),chromo.loci.index[1],replace = F)
          l <- which(t < m[first.pos:last.pos,7])
          while(length(l) < 1){t <- sample(seq(0,1,.0001),chromo.loci.index[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])}
          ch.r[[i]] <- seq(first.pos,last.pos,1)[l]
          first.pos <- last.pos+ 1}
        chr.ind.r[[each]] <- unlist(ch.r)}

      for(i in 1:cross.prog) {
        allele <- sample(1:2,1)
        end <- chromo.loci.index[num.chromos]
        z <-1
        recombination.spots <- chr.ind.r[[i]]
        for (each in 1:length(recombination.spots)) {
          if (each < length(recombination.spots)){
            gametes1[z:recombination.spots[each],i] <- par1.alleles[z:recombination.spots[each],allele]
          } else { gametes1[z:end,i] <- par1.alleles[z:end,allele] }
          if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
          z <- recombination.spots[each] + 1}}

      gametes1
    },mc.cores=num.cores)
    gametes1 <- matrix(unlist(gametes1), ncol = num.crosses*mean(as.numeric(cross.design[,3])))
    gametes2 <- mclapply(1:num.crosses,function(x){
      chromo.loci.index <- map.info$last.locus.per.chrom   # The last number of each loci for all chromsomes
      QTLSNPs      <- map.info$QTLSNP.loci # Vector of loci which are snpqtl
      m <- map.info$genetic.map

      par1<-match(cross.design[x,1],colnames(parent.info$genos.3d)) # assigns par1 to be the first parent in cross.design matrix
      par2<-match(cross.design[x,2],colnames(parent.info$genos.3d)) # assigns par2 to be the second parent in the cross.design matrix
      cross.prog<-as.numeric(cross.design[x,3]) #assigns number of progeny to be the third column for cross "X"

      # Create empty matrix to hold gametes
      # dimensions are (total # of loci) x  (# of cross progeny)
      # rownames are the loci names
      gametes2<-matrix(rep(NA,map.info$total.loci.num*cross.prog),nrow=map.info$total.loci.num,ncol=cross.prog,dimnames=list(map.info$map$loci,NULL))

      par1.alleles <- (parent.info$genos.3d[,par1,])
      par2.alleles <- (parent.info$genos.3d[,par2,])

      chr.ind.r <- vector("list")
      for(each in 1:cross.prog){
        first.pos <- 1
        ch.r <- vector("list")
        for(i in 1:length(chromo.loci.index)){
          last.pos <- chromo.loci.index[i]
          t <- sample(seq(0,1,.0001),chromo.loci.index[1],replace = F)
          l <- which(t < m[first.pos:last.pos,7])
          while(length(l) < 1){t <- sample(seq(0,1,.0001),chromo.loci.index[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])}
          ch.r[[i]] <- seq(first.pos,last.pos,1)[l]
          first.pos <- last.pos+ 1
        }
        chr.ind.r[[each]] <- unlist(ch.r)}

      for(i in 1:cross.prog) {
        allele <- sample(1:2,1)
        end <- chromo.loci.index[num.chromos]
        z <-1
        recombination.spots <- chr.ind.r[[i]]
        for (each in 1:length(recombination.spots)) {
          if (each < length(recombination.spots)){
            gametes2[z:recombination.spots[each],i] <- par2.alleles[z:recombination.spots[each],allele]
          } else { gametes2[z:end,i] <- par2.alleles[z:end,allele] }
          if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
          z <- recombination.spots[each] + 1
        }}

      gametes2
    },mc.cores=num.cores)
    gametes2 <- matrix(unlist(gametes2), ncol = num.crosses*mean(as.numeric(cross.design[,3])))
    genos.3d <- abind(gametes1,gametes2,along = 3)
    out <- list(genos.3d=genos.3d)
  } else {
    gametes1 <- lapply(1:num.crosses,function(x){
      chromo.loci.index <- map.info$last.locus.per.chrom   # The last number of each loci for all chromsomes
      QTLSNPs      <- map.info$QTLSNP.loci # Vector of loci which are snpqtl
      m <- map.info$genetic.map

      par1<-match(cross.design[x,1],colnames(parent.info$genos.3d)) # assigns par1 to be the first parent in cross.design matrix
      par2<-match(cross.design[x,2],colnames(parent.info$genos.3d)) # assigns par2 to be the second parent in the cross.design matrix
      cross.prog<-as.numeric(cross.design[x,3]) #assigns number of progeny to be the third column for cross "X"

      # Create empty matrix to hold gametes
      # dimensions are (total # of loci) x  (# of cross progeny)
      # rownames are the loci names
      gametes1<-matrix(rep(NA,map.info$total.loci.num*cross.prog),nrow=map.info$total.loci.num,ncol=cross.prog,dimnames=list(map.info$map$loci,NULL))

      par1.alleles <- (parent.info$genos.3d[,par1,])
      par2.alleles <- (parent.info$genos.3d[,par2,])
      chr.ind.r <- vector("list")
      for(each in 1:cross.prog){
        first.pos <- 1
        ch.r <- vector("list")
        for(i in 1:length(chromo.loci.index)){
          last.pos <- chromo.loci.index[i]
          t <- sample(seq(0,1,.0001),chromo.loci.index[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])
          while(length(l) < 1){t <- sample(seq(0,1,.0001),chromo.loci.index[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])}
          if(length(l) > 2){ l <- l[1:2]} else{l <- l}
          ch.r[[i]] <- seq(first.pos,last.pos,1)[l]
          first.pos <- last.pos+ 1
        }
        chr.ind.r[[each]] <- unlist(ch.r)}

      for(i in 1:cross.prog) {
        allele <- sample(1:2,1)
        end <- chromo.loci.index[num.chromos]
        z <-1
        recombination.spots <- chr.ind.r[[i]]
        for (each in 1:length(recombination.spots)) {
          if (each < length(recombination.spots)){
            gametes1[z:recombination.spots[each]-1,i] <- par1.alleles[z:recombination.spots[each]-1,allele]
          } else { gametes1[z:end,i] <- par1.alleles[z:end,allele] }
          if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
          z <- recombination.spots[each]
        }}

      gametes1
    })
    gametes1 <- matrix(unlist(gametes1), ncol = num.crosses*mean(as.numeric(cross.design[,3])))
    gametes2 <- lapply(1:num.crosses,function(x){
      chromo.loci.index <- map.info$last.locus.per.chrom   # The last number of each loci for all chromsomes
      QTLSNPs      <- map.info$QTLSNP.loci # Vector of loci which are snpqtl
      m <- map.info$genetic.map

      par1<-match(cross.design[x,1],colnames(parent.info$genos.3d)) # assigns par1 to be the first parent in cross.design matrix
      par2<-match(cross.design[x,2],colnames(parent.info$genos.3d)) # assigns par2 to be the second parent in the cross.design matrix
      cross.prog<-as.numeric(cross.design[x,3]) #assigns number of progeny to be the third column for cross "X"

      # Create empty matrix to hold gametes
      # dimensions are (total # of loci) x  (# of cross progeny)
      # rownames are the loci names
      gametes2<-matrix(rep(NA,map.info$total.loci.num*cross.prog),nrow=map.info$total.loci.num,ncol=cross.prog,dimnames=list(map.info$map$loci,NULL))

      par1.alleles <- (parent.info$genos.3d[,par1,])
      par2.alleles <- (parent.info$genos.3d[,par2,])

      chr.ind.r <- vector("list")
      for(each in 1:cross.prog){
        first.pos <- 1
        ch.r <- vector("list")
        for(i in 1:length(chromo.loci.index)){
          last.pos <- chromo.loci.index[i]
          t <- sample(seq(0,1,.0001),chromo.loci.index[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])
          while(length(l) < 1){t <- sample(seq(0,1,.0001),chromo.loci.index[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])}
          if(length(l) > 2){ l <- l[1:2]} else{l <- l}
          ch.r[[i]] <- seq(first.pos,last.pos,1)[l]
          first.pos <- last.pos+ 1
        }
        chr.ind.r[[each]] <- unlist(ch.r)}

      for(i in 1:cross.prog) {
        allele <- sample(1:2,1)
        end <- chromo.loci.index[num.chromos]
        z <-1
        recombination.spots <- chr.ind.r[[i]]
        for (each in 1:length(recombination.spots)) {
          if (each < length(recombination.spots)){
            gametes2[z:recombination.spots[each]-1,i] <- par2.alleles[z:recombination.spots[each]-1,allele]
          } else { gametes2[z:end,i] <- par2.alleles[z:end,allele] }
          if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
          z <- recombination.spots[each]
        }}

      gametes2
    })
    gametes2 <- matrix(unlist(gametes2), ncol = num.crosses*mean(as.numeric(cross.design[,3])))
    genos.3d <- abind(gametes1,gametes2,along = 3)
    out <- list(genos.3d=genos.3d)
  }
  return(out)}
