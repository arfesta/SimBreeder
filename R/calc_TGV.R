#' Cacluate Total Genetic Value
#'
#' This function esimates the total genetic value of parents produced from the create_progeny function
#' @param geno.info The object that is returned from create parents or make_crosses
#' @param map.info The object returned from create_map()
#' @param cross.design Object returned from create_cross_design()
#' @param A Value assigned to the Major SNPQTL allele
#' @param a Value assigned to the Minor SNPQTL allele
#' @param dom.coeff The dominance coeffcient used for SNPQTLs
#' @param founder logical. Only should be set to true for founder population.
#' @param prefix Name prefix to add to ouptut if save = TRUE. Default does not write output
#' @param save logical. Saves the output of genetic map (Default: FALSE)
#' @keywords progeny genetic value
#' @export
#' @examples
#' progeny1.TGV <- calc_TGV(geno.info = progeny1, map.info = the.map, crossdesign = cross.file,
#' A = 1, a = -100, dom.coeff = 1)

####Create TGV2####
calc_TGV <- function(geno.info, map.info, cross.design=NULL, A,a, dom.coeff, founder=F, rqtl.dom=NULL,
                             save=F, prefix=NULL){

  QTLSNPs <- which(map.info$types %in% c("snpqtl"))      # A vector of the loci which are snpqtl
  QTLSNP.num <- geno.info$genos.3d[QTLSNPs,,] # genotypes of both alleles pulled from the current generation
  markers <- which(map.info$types %in% c("m"))# a list of all the markers pulled from map object
  num.markers <- length(markers) # length of markers that were selected
  markers.select <- sort(sample(markers,num.markers,replace=F),decreasing=F)
  marker.select.genos <- geno.info$genos.3d[markers.select,,] # genotypes of the markers pulled from the current generation
  map.markers <- map.info[markers.select,c("chr","pos")]
  num.QTL <- length(which(map.info$types %in% c("qtl"))  ) # the number of additive qtl
  rQTL <- which(map.info$types %in% c("qtl"))  # the number of additive qtl

  if(founder == F) {par.IDs <- cross.design$parent.IDs} else {  par.IDs <- geno.info$parent.IDs}
  length.prog <- length(par.IDs)

  num.SNPQTL <- length(QTLSNPs)
  num.parents <- length(par.IDs) # the number of parents
  #QTLSNP.values <-matrix(NA,nrow=num.SNPQTL,ncol=num.parents) # matrix to hold snpqtl values
  marker.values <- matrix(NA,nrow=num.markers,ncol=num.parents) # matrix to hold marker values
  capital.genotypes <- vector()
  lowercase.genotypes <- vector()
  for (i in 1:26){
    capital.genotypes <- c(capital.genotypes,paste(LETTERS[i],LETTERS, sep=""))
    lowercase.genotypes <-  c(lowercase.genotypes,paste(letters[i],letters, sep=""))}

  # Dominance coefficient *h* of 1 means bad allele dominance, 0 mean good allele dominance
    difference <- A-a
    QTLSNPaa <- sapply(1:length.prog,function(x){
      A*length(which(QTLSNP.num[,x,1]=="a" & QTLSNP.num[,x,2]=="a"))},simplify = T)
    QTLSNPcc <- sapply(1:length.prog,function(x){
      a*length(which(QTLSNP.num[,x,1]=="c" & QTLSNP.num[,x,2]=="c"))},simplify = T)
    QTLSNPac <- sapply(1:length.prog,function(x){
      (A-(difference*dom.coeff))  * length(which(QTLSNP.num[,x,1]=="a" & QTLSNP.num[,x,2]=="c"|QTLSNP.num[,x,1]=="c" & QTLSNP.num[,x,2]=="a"))},simplify = T)
  
  QTLSNP.values <- QTLSNPaa+QTLSNPcc+QTLSNPac

  markers.aa <-   lapply(1:length.prog,function(x){which(marker.select.genos[,x,1] %in% c(LETTERS,capital.genotypes) & marker.select.genos[,x,2] %in% c(LETTERS,capital.genotypes))})
  markers.cc <-   lapply(1:length.prog,function(x){which(marker.select.genos[,x,1] %in% c(letters,lowercase.genotypes) & marker.select.genos[,x,2] %in% c(letters,lowercase.genotypes))})
  markers.ac <-   lapply(1:length.prog,function(x){which(marker.select.genos[,x,1] %in% c(LETTERS,capital.genotypes) & marker.select.genos[,x,2] %in% c(letters,lowercase.genotypes))})
  markers.ca <-   lapply(1:length.prog,function(x){which(marker.select.genos[,x,1] %in% c(letters,lowercase.genotypes) & marker.select.genos[,x,2] %in% c(LETTERS,capital.genotypes))})
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
  par.QTL.allele1 <- matrix(as.integer(geno.info$genos.3d[rQTL,,1]),nrow=num.QTL,ncol=num.parents)
  colnames(par.QTL.allele1) <- c(par.IDs)
  par.QTL.allele2 <- matrix(as.integer(geno.info$genos.3d[rQTL,,2]),nrow=num.QTL,ncol=num.parents)
  colnames(par.QTL.allele2) <- c(par.IDs)
  #QTLvalues     <- matrix(paste(par.QTL.allele1,par.QTL.allele2,sep=","),nrow=num.QTL,ncol=num.parents)
  #dimnames(QTLvalues)<-list(map.info$[map.info$rQTL],par.IDs)
  
  if(length(rqtl.dom) > 0){
    equal.rqtls <- which(par.QTL.allele1 == par.QTL.allele2)
    rqt1.better.than.rqtl2 <- which(par.QTL.allele1 > par.QTL.allele2)
    rqt2.better.than.rqtl1 <- which(par.QTL.allele2 > par.QTL.allele1)
    fill.vector <- rep(0,length(par.QTL.allele1))
    fill.vector[equal.rqtls] <- unlist(par.QTL.allele1)[c(equal.rqtls)]
    fill.vector[rqt1.better.than.rqtl2] <- unlist(par.QTL.allele1)[c(rqt1.better.than.rqtl2)] - unlist(par.QTL.allele2)[c(rqt1.better.than.rqtl2)]
    fill.vector[rqt2.better.than.rqtl1] <- unlist(par.QTL.allele2)[c(rqt2.better.than.rqtl1)] - unlist(par.QTL.allele1)[c(rqt2.better.than.rqtl1)] 
    
    rqtl.dom.coes <- c(1:length(rqtl.dom))
    rqtl.dom.coes <- rqtl.dom.coes[which(rqtl.dom.coes %in% as.numeric(names(table(ceiling(rqt1.better.than.rqtl2/ncol(par.QTL.allele1))))))]
    rqtl.dom.number <- table(ceiling(rqt1.better.than.rqtl2/ncol(par.QTL.allele1)))
    rqtl.dom.1 <- vector()
    for(each.rqtl in 1:length(rqtl.dom.coes)){
      rqtl.dom.1 <- c(rqtl.dom.1,rep(rqtl.dom.coes[each.rqtl],rqtl.dom.number[each.rqtl]))
    }

    fill.vector[rqt1.better.than.rqtl2] <- unlist(par.QTL.allele1)[c(rqt1.better.than.rqtl2)] -(fill.vector[rqt1.better.than.rqtl2]*rqtl.dom[rqtl.dom.1])
  
    rqtl.dom.coes <- c(1:length(rqtl.dom))
    rqtl.dom.coes <- rqtl.dom.coes[which(rqtl.dom.coes %in% as.numeric(names(table(ceiling(rqt2.better.than.rqtl1/ncol(par.QTL.allele1))))))]
    rqtl.dom.number <- table(ceiling(rqt2.better.than.rqtl1/ncol(par.QTL.allele1)))
    rqtl.dom.1 <- vector()
    for(each.rqtl in 1:length(rqtl.dom.coes)){
      rqtl.dom.1 <- c(rqtl.dom.1,rep(rqtl.dom.coes[each.rqtl],rqtl.dom.number[each.rqtl]))
    }
    fill.vector[rqt2.better.than.rqtl1] <- unlist(par.QTL.allele2)[c(rqt2.better.than.rqtl1)] -(fill.vector[rqt2.better.than.rqtl1]*rqtl.dom[rqtl.dom.1])
    
    fill.vector <- colSums(matrix(fill.vector,nrow=num.QTL,ncol=length.prog))
   genetic.values <-  QTLSNP.values + fill.vector
    } else {
  genetic.values <- QTLSNP.values + colSums( matrix(as.integer(geno.info$genos.3d[rQTL,,1]),nrow=num.QTL,ncol=length.prog)) + colSums(matrix(as.integer(geno.info$genos.3d[rQTL,,2]),nrow=num.QTL,ncol=length.prog))
    }
  names(genetic.values) <- par.IDs

  TGV <- list(genetic.values=genetic.values, SNP.value.matrix=QTLSNP.values,markers.matrix=marker.values, marker.loci=markers.select, marker.map=map.markers)
  return(TGV)
  }
