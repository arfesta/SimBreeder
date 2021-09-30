#' Create parents for base population
#'
#' This function creates the base population that can be used for testing further mating and selection strategies
#' @param map.info The object returned from create_map function
#' @param num.parents Number of parents that should be generated
#' @param max.delt.allele The Maximum number of deleterious alleles any single parent can have
#' @param par.markers.unqiue logical. Should each parent have a unique set of markers? Default: FALSE
#' @param heterozygous.markers logical. Should the markers for all parents be heterozygous?  Default: FALSE
#' @param QTL.sd A number providing the standard deviation of random QTL effects generated for the parents. Default: 0.25
#' @keywords create parents
#' @export
#' @examples
#' the.parents <- create_parents(map.info = the.map, num.parents = 96, max.delt.allele = 14)


####Create Founder Parent Population####
create_parents <- function(map.info, num.parents, max.delt.allele, par.markers.unique = F, heterozygous.markers=F, inbred.parents=F, QTL.sd = 0.25){

  total.loci <- as.numeric(nrow(map.info))  # Specifies the total # of loci by pulling from the map object
  locus.names <- map.info$loci  # Specifies locus names from map object
  all.MAFs <- map.info$MAF     # Assign all.MAFs as the minor allele frequencies in the map object
  total.QTL <- length(which(map.info$types %in% c("snpqtl","qtl")))   # Total number of qtl pulled from map object
  QTLSNPs <- which(map.info$types %in% c("snpqtl"))      # A vector of the loci which are snpqtl
  num.SNPQTL <- length(QTLSNPs)  # the number of snpqtl pulled from map object
  num.QTL <- total.QTL - num.SNPQTL  # the number of additive qtl (rQTL)
  QTLoci <-  which(map.info$types %in% c("qtl"))   # A vector of the loci which are additive QTL
  all.markers <- as.numeric(rownames(map.info)[-c(QTLSNPs,QTLoci)])
  num.markers <- length(all.markers)  # the number of potential markers

  #  Create Dimension names for 2 arrays
  par.IDs <- 1:num.parents; allele.IDs <- c("a1","a2")

  # Parents array holds both alleles for all parents
  # QTLSNP array holds the alleles for the loci of all parents under dominance control
  parents <- array(rep(NA,total.loci*num.parents*2), dim=c(total.loci,num.parents,2),dimnames=list(locus.names,par.IDs,allele.IDs))
  QTLSNP.array <- array(0, dim=c(num.SNPQTL,num.parents,2),dimnames=list(QTLSNPs,par.IDs,allele.IDs))

  #  Sample rQTL values from two normal distributions with mean 0 and std deviation half of that specified by QTLsd variable
  #  Use 'ceiling' function for one distribution, 'floor' function for other, then add them together to get centered distribution
  dist1 <- ceiling(rnorm((num.QTL*2), 0, (QTL.sd*(sqrt(2)/2)))); dist2 <- floor(rnorm((num.QTL*2), 0, (QTL.sd*(sqrt(2)/2))))
  QTL.alleles <- dist1 + dist2  # Vector contains the values for the rQTL alleles

  # Assign deleterious alleles
  d.marker.genotypes <- vector(); r.marker.genotypes <- vector(); marker.genotypes <- vector()
  if(par.markers.unique == T){
    if(inbred.parents == T){
      for(i in 1:26){marker.genotypes <- c(marker.genotypes, LETTERS[i],letters[i])}
    } else {
      for(i in 1:26){
        d.marker.genotypes <- c(d.marker.genotypes,paste(LETTERS[i],LETTERS, sep=""))
        r.marker.genotypes <-  c(r.marker.genotypes,paste(letters[i],letters, sep=""))
      }}} else {
      d.marker.genotypes <- rep("a",num.parents)
      r.marker.genotypes <- rep("c",num.parents)}

  #  Create empty vectors to hold either allele1 or allele2  for each parent
  allele1<-rep(NA, total.loci); allele2<-rep(NA, total.loci)

  # For each parent assign them alleles
  while(anyNA(parents)){
    if(max.delt.allele == 0 ){
      for(par in 1:num.parents){
        maf.test <- runif(n = num.SNPQTL,min = 0,max = 1)
        minor1 <- which(map.info$MAFs < maf.test)
        maf.test <- runif(n = num.SNPQTL,min = 0,max = 1)
        minor2 <- which(map.info$MAFs < maf.test)
        major1 <- QTLSNPs[which(!QTLSNPs %in% minor1)]
        major2 <- QTLSNPs[which(!QTLSNPs %in% minor2)]
        allele1[minor1]<-"c"            # minor allele is always "c", major allele is always "a"
        allele1[major1] <-"a"
        allele2[minor2]<-"c"
        allele2[major2]<-"a"

        rand1<-runif(num.markers, min=0, max=1)
        rand2<-runif(num.markers, min=0, max=1)
        minor1<-all.markers[which(all.MAFs[all.markers] > rand1)]
        major1 <- all.markers[which(!all.markers %in% minor1)]
        minor2<-all.markers[which(all.MAFs[all.markers] > rand2)]
        major2 <- all.markers[which(!all.markers %in% minor2)]

        if(heterozygous.markers){
          all.markers.allele1 <- all.markers
          minor1 <- sample(all.markers,size = length(all.markers)/2,replace=F)
          major1 <- all.markers[which(!all.markers %in% minor1)]
          minor2<- major1
          major2 <- minor1
        }
        if(inbred.parents){
          allele1[all.markers] <- marker.genotypes[par]
          allele2[all.markers] <- marker.genotypes[par]
        } else {

          allele1[minor1] <- r.marker.genotypes[par]
          allele1[major1] <- d.marker.genotypes[par]
          allele2[minor2] <- r.marker.genotypes[par]
          allele2[major2] <- d.marker.genotypes[par]}

        parents[,par,1]<-allele1
        parents[,par,2]<-allele2

        parents[QTLoci,par,1] <- sample(QTL.alleles,num.QTL)
        parents[QTLoci,par,2] <- sample(QTL.alleles,num.QTL)
        ## Write a copy of the SNP alleles at the SNPQTL into an array to use in calculating parental values
        QTLSNP.array[,par,1]<-allele1[QTLSNPs]
        QTLSNP.array[,par,2]<-allele2[QTLSNPs]
      }


      } else {
      available.delt1 <- QTLSNPs; available.delt2 <- QTLSNPs
      ratio <- (max.delt.allele/2)/num.SNPQTL

    for(par in 1:num.parents){
      num.delt <- sample(0:(num.SNPQTL*ratio),1) # pick a number of deleterious alleles for allele 1 for this parent
      if(num.delt > length(available.delt1)) {break}
      minor1 <- sample(available.delt1,num.delt)  # pick out of the avaiable deleterious alleles for this parent
      available.delt1 <- available.delt1[which(!available.delt1 %in% minor1 )] # this new list makes it so the next parent can't have any of these alleles at allele 1 positions
      major1 <- QTLSNPs[which(!QTLSNPs %in% minor1)]

      num.delt <- sample(0:(num.SNPQTL*ratio),1)
      if(num.delt > length(available.delt1[which(!available.delt1 %in% minor1)])) {break}
      minor2 <- sample(available.delt1[which(!available.delt1 %in% minor1)],num.delt)
      available.delt1 <- available.delt1[which(!available.delt1 %in% minor2 )]
      major2 <- QTLSNPs[which(!QTLSNPs %in% minor2)]

      allele1[minor1]<-"c"            # minor allele is always "c", major allele is always "a"
      allele1[major1] <-"a"
      allele2[minor2]<-"c"
      allele2[major2]<-"a"

      rand1<-runif(num.markers, min=0, max=1)
      rand2<-runif(num.markers, min=0, max=1)
      minor1<-all.markers[which(all.MAFs[all.markers] > rand1)]
      major1 <- all.markers[which(!all.markers %in% minor1)]
      minor2<-all.markers[which(all.MAFs[all.markers] > rand2)]
      major2 <- all.markers[which(!all.markers %in% minor2)]

       if(heterozygous.markers){
        all.markers.allele1 <- all.markers
        minor1 <- sample(all.markers,size = length(all.markers)/2,replace=F)
        major1 <- all.markers[which(!all.markers %in% minor1)]
        minor2<- major1
        major2 <- minor1
      }
      if(inbred.parents){
        allele1[all.markers] <- marker.genotypes[par]
        allele2[all.markers] <- marker.genotypes[par]
      } else {

      allele1[minor1] <- r.marker.genotypes[par]
      allele1[major1] <- d.marker.genotypes[par]
      allele2[minor2] <- r.marker.genotypes[par]
      allele2[major2] <- d.marker.genotypes[par]}

      parents[,par,1]<-allele1
      parents[,par,2]<-allele2

      parents[QTLoci,par,1] <- sample(QTL.alleles,num.QTL)
      parents[QTLoci,par,2] <- sample(QTL.alleles,num.QTL)
      ## Write a copy of the SNP alleles at the SNPQTL into an array to use in calculating parental values
      QTLSNP.array[,par,1]<-allele1[QTLSNPs]
      QTLSNP.array[,par,2]<-allele2[QTLSNPs]
    }}}

  # Output the snpqtls genotypes to a matrix
  parent.SNP.genos <- matrix(paste(parents[QTLSNPs,,1],parents[QTLSNPs,,2],sep=""),nrow=num.SNPQTL,ncol=num.parents)
  dimnames(parent.SNP.genos) <- list(locus.names[QTLSNPs],par.IDs)

  # Output the marker gentoypes to a matrix
  parent.markers <- matrix(paste(parents[all.markers,,1],parents[all.markers,,2],sep=""),nrow=num.markers,ncol=num.parents)
  dimnames(parent.markers) <- list(locus.names[all.markers],par.IDs)

  result.parent <- sapply(rep(1:num.parents,1),function(x) length(which(parent.SNP.genos[,x]=="ac" | parent.SNP.genos[,x]=="ca" | parent.SNP.genos[,x]=="cc")))
  result.parent <- unlist(result.parent); names(result.parent) <- colnames(parents)
  # If save is true then then we save the parentSNPgenos object and parent.markers object

  parentinfo<-list(genos.3d=parents,delt.allele=result.parent, parent.SNPQTL.matrix=parent.SNP.genos, parent.Marker.matrix=parent.markers, parent.IDs=par.IDs,num.parents=num.parents)
  #cat("The returned object is a list containing a 3-D array of parent genotypes\n")
  return(parentinfo)
  # # # # #    End of function    # # # # # # # #
}
