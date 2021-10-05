#' Generate Open Pollinated Families
#'
#' This function simulates the creation of an Open Pollinated design and can make selections of which Parents perform the best
#' @param map.info Object returned from create_genetic_map()
#' @param parent.info Object returned from create_parents()
#' @param parent.phenos Object returned from sim_parents_phenos()
#' @param parents.TGV Object returned from calc_parents_TGV()
#' @param cross.prog The number of crosses to be made between each founder parent
#' @param dom.coeff The dominance cofficient used in estimating genetic values
#' @param A The value assigned to the major allele
#' @param a The value assigned to the minor allele
#' @param h2 Heritability to be used in estimating phenotypes
#' @param n.cores The number of cores which will be used if ran in parallel
#' @keywords Generate OP families
#' @export
#' @examples
#' op.families <- OP_testing( map.info = the.map, parent.info = the.parents, parent.phenos = parent.PHENOS,
#' parents.TGV = parent.TGV, cross.prog = 1, num.select = 64, dom.coeff = 1, A = 1, a = -100, h2 = .3, 
#' run.parallel = T,n.cores = 3)


#Extract OP genotypic values####
OP_testing <- function(map.info,parent.info,parent.phenos,parents.TGV,cross.prog = 1,dom.coeff,A,a,
                       h2, E.var = NULL, n.cores= 2 ) {
  library(abind); library(parallel)

### Generate OP cross design ####
  # Retrieve number of parents
  num.pars <- parent.info$num.parents          
  
  # Create string to hold crosses for p1 and p2
  parent1 <-  rep(1:num.pars,each=(num.pars-1)); parent2 <- vector()
  
  # Fill parent 2 as every other parent except the ith parent
  for(i in 1:num.pars){
    all.parents <- 1:num.pars
    parent2 <- c(parent2,all.parents[-i])
    }
  
  # Create vector to hold number of progeny to create
  progeny <- rep(cross.prog,length(parent2))
  
  # Finish cross design for OP
  OP.crosses <- data.frame(par1=parent1,par2=parent2,num.prog=progeny)

### Data for calculating gval ####    
  # Subset the number of chromosome
  num.chromos = as.numeric(length(unique(map.info$chr)))
  
  # Get index of last loci for each chromosome
  chromo.last.loci.index <-  unlist(lapply(1:num.chromos,function(x) {
     this.chrom <- paste0("chr",x)
     max(which(map.info$chr == this.chrom))
   }))
  
  # Vector of loci which are snpqtl
  QTLSNPs      <- which(map.info$types == "snpqtl")
  
  # determine number of loci per chromosome
  loci.per.chromo <- unlist(lapply(1:num.chromos,function(x){
    this.chrom <- paste0("chr",x)
    length(which(map.info$chr == this.chrom))
  }))
 
### Estimate genetic value for each cross #### 
    g.val <- mclapply(1:nrow(OP.crosses),mc.cores = n.cores,FUN = function(x,crosses=OP.crosses){
      # subset this cross 
      y= crosses[x,]
      
      # subset parent genotypes for this cross 
      par1<-match(y[1],colnames(parent.info$genos.3d)) # assigns par1 to be the first parent in crossdesign matrix
      par2<-match(y[2],colnames(parent.info$genos.3d)) # assigns par2 to be the second parent in the crossdesign matrix
      cross.prog<-as.numeric(y[3]) #assigns number of progeny to be the third column for cross "X"
      
      # Create empty matrix to hold gametes
      # dimensions are (total # of loci) x  (# of cross progeny)
      # rownames are the loci names
      gametes1<-matrix(rep(NA,nrow(map.info)*cross.prog),nrow=nrow(map.info),ncol=cross.prog,dimnames=list(map.info$loci,NULL))
      gametes2<-matrix(rep(NA,nrow(map.info)*cross.prog),nrow=nrow(map.info),ncol=cross.prog,dimnames=list(map.info$loci,NULL))
      
      par1.alleles <- (parent.info$genos.3d[,par1,])
      par2.alleles <- (parent.info$genos.3d[,par2,])
      
      # Create recombination breakpoints for each cross
      chr.ind.r <- vector("list")
      for(each in 1:cross.prog){ # For each progeny we are going to do the following
        first.pos <- 1           # Specify the first position to be row 1 on the genetic map
        ch.r <- vector("list")   # Create empty vector list to hold recombination spots
        
        for(i in 1:length(chromo.last.loci.index)){ # Now for each chromosome do the following
          last.pos <- chromo.last.loci.index[i]     # Subset the last position of this chromosome
          l <- which((rbinom(n = (first.pos+1):last.pos,size = 1, prob = map.info$dist[(first.pos+1):last.pos])) == 1)
          #t <- sample(seq(0,1,.0001),loci.per.chromo[i],replace = F)  # Sample a range of frequencies to test against rec. freqs
          #l <- which(t < map.info$recfreqs[first.pos:last.pos])  # Identify which of the sampled numbers were less than the rec. freq
          while(length(l) < 1){   # If none were less than rec. freq, sample again
            #t <- sample(seq(0,1,.0001),loci.per.chromo[i],replace = F) 
            #l <- which(t < map.info$recfreqs[first.pos:last.pos])
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

        
      # Identify recombination spots for parent 2
      chr.ind.r <- vector("list")
      for(each in 1:cross.prog){ # For each progeny we are going to do the following
        first.pos <- 1           # Specify the first position to be row 1 on the genetic map
        ch.r <- vector("list")   # Create empty vector list to hold recombination spots
        
        for(i in 1:length(chromo.last.loci.index)){ # Now for each chromosome do the following
          last.pos <- chromo.last.loci.index[i]     # Subset the last position of this chromosome
          l <- which((rbinom(n = (first.pos+1):last.pos,size = 1, prob = map.info$dist[(first.pos+1):last.pos])) == 1)
          #t <- sample(seq(0,1,.0001),loci.per.chromo[i],replace = F)  # Sample a range of frequencies to test against rec. freqs
          #l <- which(t < map.info$recfreqs[first.pos:last.pos])  # Identify which of the sampled numbers were less than the rec. freq
          while(length(l) < 1){   # If none were less than rec. freq, sample again
            #t <- sample(seq(0,1,.0001),loci.per.chromo[i],replace = F) 
            #l <- which(t < map.info$recfreqs[first.pos:last.pos])
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
      array.out <- abind(gametes1,gametes2,along=3)
      
      locus.names <- map.info$loci # The locus names pulled from the mab object
      QTLSNPs      <- which(map.info$types == "snpqtl")    # vector of the loci which are snpqtl
      QTLSNP.num <- array.out[QTLSNPs,,] # genotypes of both alleles pulled from the current generation
      rQTL.loci      <- which(map.info$types == "qtl")
      num.QTL <- length(rQTL.loci) # the number of additive qtl
      
      # Create 2 matrices: One to hold snpqtl values for all parents and the other to hold marker marker values for blup analysis
      num.SNPQTL <- length(QTLSNPs) # the number of loci which are snpqtl
      QTLSNP.values <-matrix(NA,nrow=num.SNPQTL,ncol=1) # matrix to hold snpqtl values
      
     # Dominance coefficient *h* of 1 means bad allele dominance, 0 mean good allele dominance
    difference <- A-a
    QTLSNPaa <- A*length(which(QTLSNP.num[,1]=="a" & QTLSNP.num[,2]=="a"))
    QTLSNPcc <- a*length(which(QTLSNP.num[,1]=="c" & QTLSNP.num[,2]=="c"))
    QTLSNPac <- (A-(difference*dom.coeff))  * length(which(QTLSNP.num[,1]=="a" & QTLSNP.num[,2]=="c"|QTLSNP.num[,1]=="c" & QTLSNP.num[,2]=="a"))
  
  QTLSNP.values <- QTLSNPaa+QTLSNPcc+QTLSNPac
      
      par.QTL.allele1 <- as.numeric(array.out[rQTL.loci,,1])
      par.QTL.allele2 <- as.numeric(array.out[rQTL.loci,,2])
      
      
      # Genetic values of progeny
      genetic.value <- sum(QTLSNP.values) + sum(par.QTL.allele1) + sum(par.QTL.allele2)
      genetic.value
    })
  
  g.val <- unlist(g.val)
  total.indiv <- length(g.val)
  if(length(E.var) > 0){ phenos <- g.val + rnorm(total.indiv,mean = 0,sd = E.var) } else {
    phenos <- g.val + rnorm(total.indiv,mean = 0,sd = sqrt(var(g.val)/h2 - var(g.val)))
    }
  
  mean.gval <- vector()
  first <- 1
  last <- num.pars-1
  for(i in 1:num.pars){
    mean.gval <- c(mean.gval,mean(g.val[first:last]))
    first <- last+1
    last <- first + num.pars-2
    } 
  
  
  mean.phenos <- vector()
  first <- 1
  last <- num.pars-1
  for(i in 1:num.pars){
    mean.phenos <- c(mean.phenos,mean(phenos[first:last]))
    first <- last+1
    last <- first + num.pars-2} 
  
  names(parent.phenos$phenos) <- colnames(parent.info$genos.3d)
  names(mean.phenos) <- colnames(parent.info$genos.3d)
  pars <- sort(mean.phenos,decreasing = T)
  mean.gval = mean.gval[as.numeric(names(pars))]
  delt.alleles <- parent.info$delt.allele[as.numeric(names(pars))]
  genetic.values <- parent.phenos$genetic.values[as.numeric(names(pars))]
  phenos <- parent.phenos$phenos[as.numeric(names(pars))]
  genos.3d <- parent.info$genos.3d[,as.numeric(names(pars)),]
  marker.matrix <- parents.TGV$markers.matrix[as.numeric(names(pars)),]
  
  out <- list(delt.alleles=delt.alleles,genetic.values=genetic.values,phenos=phenos, genos.3d=genos.3d, mean.parent.phenos = pars,mean.parent.tgv = mean.gval,
              marker.matrix=marker.matrix)
  return(out)
}