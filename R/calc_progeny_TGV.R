#' Cacluate Progeny Total Genetic Value
#'
#' This function esimates the total genetic value of parents produced from the create_progeny function
#' @param prog.info The object that is returned from the make_crosses()
#' @param map.info The object returned from create_map()
#' @param cross.design Object returned from create_cross_design()
#' @param A Value assigned to the Major SNPQTL allele
#' @param a Value assigned to the Minor SNPQTL allele
#' @param dom.coeff The dominance coeffcient used for SNPQTLs
#' @param usesnpeffect logical.  Set to TRUE if specificying user provided snp effects. Default: FALSE
#' @param snp.effects The vector of SNP effects to use. Only use if usesnpeffect = TRUE.
#' @param prefix Name prefix to add to ouptut if save = TRUE. Default does not write output
#' @param save logical. Saves the output of genetic map (Default: FALSE)
#' @keywords progeny genetic value
#' @export
#' @examples
#' progeny1.TGV <- calc_progeny_TGV(prog.info = progeny1, map.info = the.map, crossdesign = cross.file,
#' A = 1, a = -100, dom.coeff = 1)

####Create TGV2####
calc_progeny_TGV <- function(prog.info, map.info, cross.design, A,a, dom.coeff,
                             save=F, prefix=NULL){

  QTLSNPs <- map.info$QTLSNP.loci     # vector of the loci which are snpqtl
  QTLSNP.num <- prog.info$genos.3d[QTLSNPs,,] # genotypes of both alleles pulled from the current generation
  markers <- map.info$available.Markers# a list of all the markers pulled from map object
  num.markers <- length(markers) # length of markers that were selected
  markers.select <- sort(sample(markers,num.markers,replace=F),decreasing=F)
  marker.select.genos <- prog.info$genos.3d[markers.select,,] # genotypes of the markers pulled from the current generation
  map.markers <- map.info$genetic.map[markers.select,c(1,6)]
  num.QTL <- length(map.info$rQTL.loci) # the number of additive qtl

  par.IDs <- cross.design$parent.IDs
  length.prog <- length(par.IDs)

  num.SNPQTL <- map.info$total.SNPQTL.num # the number of loci which are snpqtl
  num.parents <- length(par.IDs) # the number of parents
  #QTLSNP.values <-matrix(NA,nrow=num.SNPQTL,ncol=num.parents) # matrix to hold snpqtl values
  marker.values <- matrix(NA,nrow=num.markers,ncol=num.parents) # matrix to hold marker values
  capital.genotypes <- vector()
  lowercase.genotypes <- vector()
  for (i in 1:26){
    capital.genotypes <- c(capital.genotypes,paste(LETTERS[i],LETTERS, sep=""))
    lowercase.genotypes <-  c(lowercase.genotypes,paste(letters[i],letters, sep=""))}

  QTLSNPaa <- sapply(1:length.prog,function(x){
    A*length(which(QTLSNP.num[,x,1]=="a" & QTLSNP.num[,x,2]=="a"))},simplify = T)
  QTLSNPcc <- sapply(1:length.prog,function(x){
    a*length(which(QTLSNP.num[,x,1]=="c" & QTLSNP.num[,x,2]=="c"))},simplify = T)
  QTLSNPac <- sapply(1:length.prog,function(x){
    (A * dom.coeff) * length(which(QTLSNP.num[,x,1]=="a" & QTLSNP.num[,x,2]=="c"))},simplify = T)
  QTLSNPca <- sapply(1:length.prog,function(x){
    (A * dom.coeff) * length(which(QTLSNP.num[,x,1]=="c" & QTLSNP.num[,x,2]=="a"))},simplify = T)
  QTLSNP.values <- QTLSNPaa+QTLSNPcc+QTLSNPac+QTLSNPca

  markers.aa <-   sapply(1:length.prog,function(x){which(marker.select.genos[,x,1] %in% capital.genotypes & marker.select.genos[,x,2] %in% capital.genotypes)},simplify = T)
  markers.cc <-   sapply(1:length.prog,function(x){which(marker.select.genos[,x,1] %in% lowercase.genotypes & marker.select.genos[,x,2] %in% lowercase.genotypes)},simplify = T)
  markers.ac <-   sapply(1:length.prog,function(x){which(marker.select.genos[,x,1] %in% capital.genotypes & marker.select.genos[,x,2] %in% lowercase.genotypes)},simplify = T)
  markers.ca <-   sapply(1:length.prog,function(x){which(marker.select.genos[,x,1] %in% lowercase.genotypes & marker.select.genos[,x,2] %in% capital.genotypes)},simplify = T)
  marker.values <- matrix(NA,nrow=num.markers,ncol=length.prog) # matrix to hold marker values

  for(i in 1:length.prog){
    marker.values[markers.aa[[i]],i] <- "0"
    marker.values[markers.cc[[i]],i] <- "2"
    marker.values[markers.ac[[i]],i] <- "1"
    marker.values[markers.ca[[i]],i] <- "1"}

  marker.values <- t(marker.values)
  colnames(marker.values) <- markers.select
  rownames(marker.values) <- par.IDs
  #removed marker.vals from return (markers.matrix=marker.values)
  par.QTL.allele1 <- matrix(as.integer(prog.info$genos.3d[map.info$rQTL,,1]),nrow=num.QTL,ncol=num.parents)
  colnames(par.QTL.allele1) <- c(par.IDs)
  par.QTL.allele2 <- matrix(as.integer(prog.info$genos.3d[map.info$rQTL,,2]),nrow=num.QTL,ncol=num.parents)
  colnames(par.QTL.allele2) <- c(par.IDs)
  #QTLvalues     <- matrix(paste(par.QTL.allele1,par.QTL.allele2,sep=","),nrow=num.QTL,ncol=num.parents)
  #dimnames(QTLvalues)<-list(map.info$[map.info$rQTL],par.IDs)
  genetic.values <- QTLSNP.values + colSums( matrix(as.integer(prog.info$genos.3d[map.info$rQTL,,1]),nrow=num.QTL,ncol=length.prog)) + colSums(matrix(as.integer(prog.info$genos.3d[map.info$rQTL,,2]),nrow=num.QTL,ncol=length.prog))
  names(genetic.values) <- par.IDs

  TGV <- list(genetic.values=genetic.values, SNP.value.matrix=QTLSNP.values,markers.matrix=marker.values, marker.loci=markers.select, marker.map=map.markers)
  return(TGV)}
