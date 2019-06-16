#' Cacluate Parents Total Genetic Value
#'
#' This function esimates the total genetic value of parents produced from the create_parents function
#' @param parents The object that is returned from the create_parents function
#' @param map.info The object returned from create_map function
#' @param A Value assigned to the Major SNPQTL allele
#' @param a Value assigned to the Minor SNPQTL allele
#' @param dom.coeff The dominance coeffcient used for SNPQTLs
#' @param usesnpeffect logical.  Set to TRUE if specificying user provided snp effects. Default: FALSE
#' @param snp.effects The vector of SNP effects to use. Only use if usesnpeffect = TRUE.
#' @param prefix Name prefix to add to ouptut if save = TRUE. Default does not write output
#' @param save logical. Saves the output of genetic map (Default: FALSE)
#' @keywords parent genetic value
#' @export
#' @examples
#' parent.TGV <- calc_parents_TGV(parents = the.parents, map.info = the.map, A = 1, a = -100, dom.coeff = 1)

####Calc Parent TGV Genetic Value####
calc_parents_TGV <- function(parents, map.info,
                             A, a, dom.coeff,
                             usesnpeffect = F,snp.effects = NULL,
                             save = F, prefix = NULL) {
  
  locus.names <- map.info$genetic.map$loci # The locus names pulled from the mab object
  QTLSNPs <- map.info$QTLSNP.loci     # vector of the loci which are snpqtl
  QTLSNP.num <- parents$genos.3d[QTLSNPs,,] # genotypes of both alleles pulled from the current generation
  markers <- map.info$available.Markers# a list of all the markers pulled from map object
  num.markers <- length(markers)
  marker.loci <- sort(sample(markers,num.markers,replace=F),decreasing=F)
  num.markers <- length(marker.loci) # length of markers that were selected
  marker.select.genos <- parents$genos.3d[marker.loci,,] # genotypes of the markers pulled from the current generation
  marker.map <- map.info$genetic.map[marker.loci,c(1,6)]
  num.QTL <- length(map.info$rQTL.loci) # the number of additive qtl
  par.IDs <- parents$parent.IDs
  
  
  
  # Create 2 matrices: One to hold snpqtl values for all parents and the other to hold marker marker values for blup analysis
  num.SNPQTL <- map.info$total.SNPQTL.num # the number of loci which are snpqtl
  num.parents <- length(par.IDs) # the number of parents
  QTLSNP.values <-matrix(NA,nrow=num.SNPQTL,ncol=num.parents) # matrix to hold snpqtl values
  marker.matrix <- matrix(NA,nrow=num.markers,ncol=num.parents) # matrix to hold marker values
  capital.genotypes <- vector()
  lowercase.genotypes <- vector()
  for (i in 1:26){
    capital.genotypes <- c(capital.genotypes,paste(LETTERS[i],LETTERS, sep=""))
    lowercase.genotypes <-  c(lowercase.genotypes,paste(letters[i],letters, sep=""))}
  
  for (i in 1:length(par.IDs)){
    QTLSNPaa <- which(QTLSNP.num[,i,1]=="a" & QTLSNP.num[,i,2]=="a")
    QTLSNPcc <- which(QTLSNP.num[,i,1]=="c" & QTLSNP.num[,i,2]=="c")
    QTLSNPac <- which(QTLSNP.num[,i,1]=="a" & QTLSNP.num[,i,2]=="c")
    QTLSNPca <- which(QTLSNP.num[,i,1]=="c"  & QTLSNP.num[,i,2]=="a")
    if (dom.coeff==0){
      QTLSNP.values[QTLSNPaa,i] <- A*2
      QTLSNP.values[QTLSNPcc,i] <- a*2
      QTLSNP.values[QTLSNPac,i] <- (A+a+dom.coeff)
      QTLSNP.values[QTLSNPca,i] <- (A+a+dom.coeff) } else {
        QTLSNP.values[QTLSNPaa,i] <- A
        if(usesnpeffect){
          QTLSNP.values[QTLSNPcc,i] <- snp.effects[a]
        } else{
          QTLSNP.values[QTLSNPcc,i] <- a}
        QTLSNP.values[QTLSNPac,i] <- (A*dom.coeff)
        QTLSNP.values[QTLSNPca,i] <- (A*dom.coeff)
      }
    
    markers.aa <- which(marker.select.genos[,i,1] %in% capital.genotypes & marker.select.genos[,i,2] %in% capital.genotypes)
    markers.cc <- which(marker.select.genos[,i,1] %in% lowercase.genotypes & marker.select.genos[,i,2] %in% lowercase.genotypes)
    markers.ac <-which(marker.select.genos[,i,1] %in% capital.genotypes & marker.select.genos[,i,2] %in% lowercase.genotypes)
    markers.ca <- which(marker.select.genos[,i,1] %in% lowercase.genotypes & marker.select.genos[,i,2] %in% capital.genotypes)
    
    marker.matrix[markers.aa,i] <- "0"
    marker.matrix[markers.cc,i] <- "2"
    marker.matrix[markers.ac,i] <- "1"
    marker.matrix[markers.ca,i] <- "1"
  }
  marker.matrix <- t(marker.matrix)
  colnames(marker.matrix) <- marker.loci
  rownames(marker.matrix) <- par.IDs
  
  # Convert the 'invisible' rQTL genotypes to numeric matrices, merge the alleles to paired values also
  par.QTL.allele1 <- matrix(as.integer(parents$genos.3d[map.info$rQTL,,1]),nrow=num.QTL,ncol=num.parents)
  colnames(par.QTL.allele1) <- c(par.IDs)
  par.QTL.allele2 <- matrix(as.integer(parents$genos.3d[map.info$rQTL,,2]),nrow=num.QTL,ncol=num.parents)
  colnames(par.QTL.allele2) <- c(par.IDs)
  QTL.values     <- matrix(paste(par.QTL.allele1,par.QTL.allele2,sep=","),nrow=num.QTL,ncol=num.parents)
  dimnames(QTL.values)<-list(locus.names[map.info$rQTL],par.IDs)
  
  
  genetic.values <- colSums(QTLSNP.values) + colSums(par.QTL.allele1) + colSums(par.QTL.allele2) 
  #parQTLalleles <- par.QTL.allele2 + rqtldom
  # Genetic values of progeny
  #geneticvals <- colSums(QTLSNP.values) + colSums(final.qtl)
  
  if(save) {
    pedfilename=paste(prefix,".txt", sep="")
    write.table(genetic.values,pedfilename, quote = F, row.names = F, col.names = T, sep="\t")}
  
  TGV <- list(genetic.values=genetic.values, SNP.value.matrix=QTLSNP.values, marker.matrix=marker.matrix, marker.loci=marker.loci, marker.map=marker.map)
  return(TGV)
}
