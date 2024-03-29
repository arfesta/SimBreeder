#' Make crosses
#'
#' This function makes crosses among parents in a format specified by the cross design file
#' @param parent.info Object returned from create_parents() or extract_selections()
#' @param map.info Object returned from create_genetic_map()
#' @param cross.design Object returned from create_cross_design()
#' @param num.cores The number of cores to use if running in parallel
#' @keywords create crosses among parents
#' @export
#' @examples
#' progeny1 <- make_crosses(parent.info = op.families,map.info = the.map,cross.design = cross.file,run.parallel = T,num.cores = 3)

####Create Make Crosses####
make_crosses <- function(parent.info,map.info,cross.design, num.cores = num.of.cores){
  library(parallel); library(abind)
  cross_design <- cross.design$cross.design
  num.crosses <- as.numeric(nrow(cross_design))
  num.chromos <- length(unique(map.info$chr))
  last.locus.per.chrom <- vector()
  for(each in 1:num.chromos){
    this.chr <- map.info$loci[which(map.info$chr == unique(map.info$chr)[each])]
    last.locus.per.chrom <- c(last.locus.per.chrom,which(map.info$loci == this.chr[length(this.chr)]))
  }
  loci.per.chromo <- unlist(lapply(1:num.chromos,function(x){
    this.chrom <- paste0("chr",x)
    length(which(map.info$chr == this.chrom))
  }))
  chromo.last.loci.index <-  unlist(lapply(1:num.chromos,function(x) {
    this.chrom <- paste0("chr",x)
    max(which(map.info$chr == this.chrom))
  }))
  
  gametes1 <- mclapply(1:num.crosses,function(x){
    chromo.loci.index <- last.locus.per.chrom   # The last number of each loci for all chromsomes
    QTLSNPs      <- which(map.info$types %in% c("snpqtl"))      # A vector of the loci which are snpqtl
    
    par1<-match(cross_design[x,1],unique(colnames(parent.info$genos.3d))) # assigns par1 to be the first parent in cross.design matrix
    par2<-match(cross_design[x,2],unique(colnames(parent.info$genos.3d))) # assigns par2 to be the second parent in the cross.design matrix
    cross.prog<-as.numeric(cross_design[x,3]) #assigns number of progeny to be the third column for cross "X"
    
    # Create empty matrix to hold gametes
    # dimensions are (total # of loci) x  (# of cross progeny)
    # rownames are the loci names
    gametes1<-matrix(rep(NA,nrow(map.info)*cross.prog),nrow=nrow(map.info),ncol=cross.prog,dimnames=list(map.info$loci,NULL))
    
    if(length(dim(parent.info$genos.3d)) == 2) { par1.alleles <- parent.info$genos.3d; par2.alleles <- parent.info$genos.3d} else {
      par1.alleles <- (parent.info$genos.3d[,par1,])
      par2.alleles <- (parent.info$genos.3d[,par2,])
    }
    chr.ind.r <- vector("list")
    for(each in 1:cross.prog){ # For each progeny we are going to do the following
      first.pos <- 1           # Specify the first position to be row 1 on the genetic map
      ch.r <- vector("list")   # Create empty vector list to hold recombination spots
      
      for(i in 1:length(chromo.last.loci.index)){ # Now for each chromosome do the following
        last.pos <- chromo.last.loci.index[i]     # Subset the last position of this chromosome
        l <- which((rbinom(n = (first.pos+1):last.pos,size = 1, prob = map.info$dist[(first.pos+1):last.pos])) == 1)
        while(length(l) < 1){   # If none were less than rec. freq, sample again
          l <- which((rbinom(n = (first.pos+1):last.pos,size = 1, prob = map.info$dist[(first.pos+1):last.pos])) == 1)
        }
        ch.r[[i]] <- seq(first.pos,last.pos,1)[l+1] # Now subest those positions from the loci to find which loci will recombine
        first.pos <- last.pos+ 1  # The first position will now be the last + 1
      }
      
      chr.ind.r[[each]] <- unlist(ch.r)} # Store output of recombination spots
    
    # Assign gametes 1
    for(i in 1:cross.prog) {
      allele <- sample(1:2,num.chromos,replace = T) # sample 12 alleles to start for each chromosome
      rec.spots <- chr.ind.r[[i]]          #subset the list of recombination spots for this indivdual
      first.locus <- 1                     # Always start with chr1
      for(each.chromo in 1:num.chromos){   # For each chromosome do the following
        # Subset the recspots for this indvidiual that are on this chromsome and the starting allele
        chr.rec.spots <- rec.spots[which(rec.spots >= first.locus & rec.spots <= chromo.last.loci.index[each.chromo])]
        starting.allele <- allele[each.chromo]
        for(each.rec.spot in 1:length(chr.rec.spots)){ # Then, for each of those spots do the following
          # Specify gametes to be the first poistion up to the rec spot
          gametes1[first.locus:chr.rec.spots[each.rec.spot],i] <- par1.alleles[first.locus:chr.rec.spots[each.rec.spot],starting.allele]
          # Add one to the rec spot so that next set will get next locus
          first.locus <- chr.rec.spots[each.rec.spot] + 1
          # Add one to the starting allele so it will switch sides
          if (starting.allele==1) { starting.allele <- starting.allele +1 } else { starting.allele <- starting.allele -1}
        }
        if(chr.rec.spots[length(chr.rec.spots)] != chromo.last.loci.index[each.chromo]){
          gametes1[first.locus:chromo.last.loci.index[each.chromo],i] <- par1.alleles[first.locus:chromo.last.loci.index[each.chromo],starting.allele]
          first.locus <- chromo.last.loci.index[each.chromo] + 1
        }
      }
    }
    
    gametes1
  },mc.cores=num.cores)
  gametes1 <- matrix(unlist(gametes1), ncol = num.crosses*mean(as.numeric(cross_design[,3])))
  gametes2 <- mclapply(1:num.crosses,function(x){
    chromo.loci.index <- last.locus.per.chrom   # The last number of each loci for all chromsomes
    QTLSNPs      <- which(map.info$types %in% c("snpqtl"))      # A vector of the loci which are snpqtl
    
    par1<-match(cross_design[x,1],unique(colnames(parent.info$genos.3d))) # assigns par1 to be the first parent in cross_design matrix
    par2<-match(cross_design[x,2],unique(colnames(parent.info$genos.3d))) # assigns par2 to be the second parent in the cross_design matrix
    cross.prog<-as.numeric(cross_design[x,3]) #assigns number of progeny to be the third column for cross "X"
    
    # Create empty matrix to hold gametes
    # dimensions are (total # of loci) x  (# of cross progeny)
    # rownames are the loci names
    gametes2<-matrix(rep(NA,nrow(map.info)*cross.prog),nrow=nrow(map.info),ncol=cross.prog,dimnames=list(map.info$loci,NULL))
    if(length(dim(parent.info$genos.3d)) == 2) { par1.alleles <- parent.info$genos.3d; par2.alleles <- parent.info$genos.3d} else {
      par1.alleles <- (parent.info$genos.3d[,par1,])
      par2.alleles <- (parent.info$genos.3d[,par2,])}
    
    chr.ind.r <- vector("list")
    for(each in 1:cross.prog){ # For each progeny we are going to do the following
      first.pos <- 1           # Specify the first position to be row 1 on the genetic map
      ch.r <- vector("list")   # Create empty vector list to hold recombination spots
      
      for(i in 1:length(chromo.last.loci.index)){ # Now for each chromosome do the following
        last.pos <- chromo.last.loci.index[i]     # Subset the last position of this chromosome
        l <- which((rbinom(n = (first.pos+1):last.pos,size = 1, prob = map.info$dist[(first.pos+1):last.pos])) == 1)
        while(length(l) < 1){   # If none were less than rec. freq, sample again
          l <- which((rbinom(n = (first.pos+1):last.pos,size = 1, prob = map.info$dist[(first.pos+1):last.pos])) == 1)
        }
        ch.r[[i]] <- seq(first.pos,last.pos,1)[l+1] # Now subest those positions from the loci to find which loci will recombine
        first.pos <- last.pos+ 1  # The first position will now be the last + 1
      }
      
      chr.ind.r[[each]] <- unlist(ch.r)} # Store output of recombination spots
    
    # Assign gametes 2
    for(i in 1:cross.prog) {
      allele <- sample(1:2,num.chromos,replace = T) # sample 12 alleles to start for each chromosome
      rec.spots <- chr.ind.r[[i]]          #subset the list of recombination spots for this indivdual
      first.locus <- 1                     # Always start with chr1
      for(each.chromo in 1:num.chromos){   # For each chromosome do the following
        # Subset the recspots for this indvidiual that are on this chromsome and the starting allele
        chr.rec.spots <- rec.spots[which(rec.spots >= first.locus & rec.spots <= chromo.last.loci.index[each.chromo])]
        starting.allele <- allele[each.chromo]
        for(each.rec.spot in 1:length(chr.rec.spots)){ # Then, for each of those spots do the following
          # Specify gametes to be the first poistion up to the rec spot
          gametes2[first.locus:chr.rec.spots[each.rec.spot],i] <- par2.alleles[first.locus:chr.rec.spots[each.rec.spot],starting.allele]
          # Add one to the rec spot so that next set will get next locus
          first.locus <- chr.rec.spots[each.rec.spot] + 1
          # Add one to the starting allele so it will switch sides
          if (starting.allele==1) { starting.allele <- starting.allele +1 } else { starting.allele <- starting.allele -1}
        }
        if(chr.rec.spots[length(chr.rec.spots)] != chromo.last.loci.index[each.chromo]){
          gametes2[first.locus:chromo.last.loci.index[each.chromo],i] <- par2.alleles[first.locus:chromo.last.loci.index[each.chromo],starting.allele]
          first.locus <- chromo.last.loci.index[each.chromo] + 1
        }
      }
    }
    
    gametes2
  },mc.cores=num.cores)
  gametes2 <- matrix(unlist(gametes2), ncol = num.crosses*mean(as.numeric(cross_design[,3])))
  genos.3d <- abind(gametes1,gametes2,along = 3)
  colnames(genos.3d) <- cross.design$progeny.pedigree[,1]
  out <- list(genos.3d=genos.3d)
  return(out)}
