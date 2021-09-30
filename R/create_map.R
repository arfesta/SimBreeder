#' Simulate Genetic Map
#'
#' This function simulates the creation of a genetic map which mimics the genetic architecture of user defined inputs
#' @param num.chromos The number of chromosomes to simulate
#' @param map.length The genetic map length in centi-morgans
#' @param num.markers Number of markers that can be used for genomic selection or estimated relatedness
#' @param total.qtl The total number of QTL to place on the genetic map
#' @param num.snpqtl The number of total.qtl which will be SNPs
#' @param map.dist Type of distance function to use in calculating recombination frequences.  Options are: "haldane" & "kosambi".  (Default = "haldane)
#' @param distribute.loci Default is NULL so that loci are randomly distributed to linkage groups. Options: "even" or "list"
#' @param marker.distribution Changes the assignment of markers to linkage groups from random to equally spaced ("equally-spaced")
#' @param snp.qtl.maf A vector of 2 numbers indiciating the min and max range of SNP QTL MAFs
#' @param marker.maf Default is Null. Unless specified markers MAFs will be drawn from beta distribution
#' @param chromosome.size.range The chromosome size range from the mean length of chromosomes (Default: 0.2)
#' @param seed logical. Set seed so that genetic map is the same every time (Default: FALSE)
#' @param save logical. Saves the output of genetic map (Default: FALSE)
#' @keywords cat
#' @export
#' @examples
#' create_genetic_map(num.chromos = 12, map.length = 1800, num.markers = 120, total.loci = signif((120+2500)/12,digits = 2)*12, total.qtl = 2500, num.snpqtl = 1960)

#Create Map Function####
create_genetic_map <- function(num.chromos, map.length, num.markers, total.qtl, num.snpqtl,
                               map.dist = "haldane", distribute.loci = NULL, marker.distribution =NULL,
                               total.loci = NULL, marker.maf=NULL, snp.qtl.maf= NULL, chromosome.size.range=.2, seed=F, save=F) {

  #Determine average length of a chromosome
  avg.chromo.size <- map.length / num.chromos

  #Use mean chromosome size length and the chromosomse size range variable to specify chromosome sizes
  all.chromo.sizes <- floor(runif(num.chromos, min=(avg.chromo.size - chromosome.size.range * avg.chromo.size),
                                  max=(avg.chromo.size + chromosome.size.range * avg.chromo.size)))

  #Determine the number of intervals that will be calculated between loci and the number of loci per chromosome
 total.loci <- sum(num.markers,total.qtl)

  num.intervals <- total.loci - 1 #Specifying the number of intervals that will need to be calcuated
  if(is.null(distribute.loci)) { loci.per.chromo <- round(rep(total.loci / num.chromos, num.chromos),digits = 0)
  } else if(distribute.loci == "list") {
    loci.per.chromo <- unlist(loci.per.chromo)
  } else if(distribute.loci == "even") {
    loci.per.chromo <- rep(total.loci / num.chromos, num.chromos)#Each chromosome will get an equal number of loci
    }

  #Create Map and Recombination Frequencies####
  # Define map function used to create recombination frequencies:
  calc_rec_freq <- function(x,mapdist=map.dist,interval=intervals) {
    loci.intervals <- interval[[x]]
    if(map.dist == "haldane") { out <- (1 - (exp(-2 * loci.intervals))) / 2 }
    if(map.dist == "kosambi") { out <- ((exp(2 * loci.intervals) - exp(-2 * loci.intervals)) / (exp(2 * loci.intervals) + exp(-2 * loci.intervals))) / 2 }
    out
  }

  # Set up genetic map data frame
  list1 <- lapply(1:num.chromos,function(x) rep(paste("chr",x,sep=""),loci.per.chromo[x])) # List "chr#" character set for all loci
  list2 <- lapply(1:num.chromos,function(x) 1:loci.per.chromo[x]) # List locus# for all chromosomes
  intervals <- lapply(1:num.chromos,function(x) c(round(runif((loci.per.chromo[x]-1),min=round(0.2*all.chromo.sizes[x]/loci.per.chromo[x],1),
                                                              max=round(1.8*all.chromo.sizes[x]/loci.per.chromo[x],1)),2)/100,5))
  rec.freqs <- lapply(1:num.chromos,calc_rec_freq)
  #Recfreqs for each of the loci are determined by using the interval distance (for each locus) in Morgans and the haldane mapping function:
  #r=(1/2)*(1-e^(2*M))

  #Since intervals are in Morgans, we multiply by 100 to turn it back into cM for positions (cM is unit for chromosome lengths)
  positions <- lapply(1:num.chromos,function(x) {
    all.except.last <- c(round(100*cumsum(intervals[[x]][-length(intervals[[x]])]),digits=8))
    last <- abs(all.chromo.sizes[x]/100 + all.except.last[length(all.except.last)])
    c(all.except.last,last)})


  ## paste together locus names from list2 and list1
  locus.names<- paste(unlist(list1),unlist(list2),sep="_")

  #Create map that has locusnames, distance between loci, empty vector of types (SNPQTL,rQTL, or Marker), MAFs, position on chromosome, & rec freqs
  map <- data.frame(chr=unlist(list1),loci=locus.names, dist= (unlist(intervals)),types= (rep("NA",total.loci)), MAF=rep("NA",total.loci),
                    pos=as.numeric(unlist(positions)),recfreqs=as.numeric(unlist(rec.freqs)),stringsAsFactors = F)

  #Now sample from the map to specify SNPQTL & rQTL, the remainder are potential markers that can be used
  all.loci <- 1:total.loci  # vector that contains 1 through the number of all loci

  #Sample 10 markers for each chromosome uniformly distributed####
  #Specify the last loci for each chromosome
  chromo.loci.index <- vector("list")
  smp <- num.markers/num.chromos
  for (i in 1:length(loci.per.chromo)) {
    if (i==1) {
      chromo.loci.index[[1]] <- loci.per.chromo[i]
    } else {
      chromo.loci.index[[i]] <- loci.per.chromo[[i]] + chromo.loci.index[[i-1]]
    }
  }
  chromo.loci.index <- unlist(chromo.loci.index)
  last.pos.chrs <-  sapply(1:num.chromos,function(x){round(map$pos[chromo.loci.index[x]-1],0)})
  markers <- vector(); start <- 1
  if(is.null(marker.distribution)) { markers <- sample(x = 1:nrow(map),size = num.markers,replace = F)
  } else if(marker.distribution == "equally-spaced"){
    for(i in 1:num.chromos){
    equal.dist.per.chr <- last.pos.chrs[i]/smp
    map.posit <- map$pos[start:chromo.loci.index[i]]
    marker <- sapply(1:smp,function(x){
      the.pos <- equal.dist.per.chr*x
      which.min(abs(round(map.posit - the.pos,1)))})
    if(i > 1){marker <- marker + chromo.loci.index[i-1]}
    markers <- c(markers,marker)
    start <- chromo.loci.index[i] + 1} }

  map$types[markers] <- "m"
  # Specify in map data frame that all loci which are not qtl or snpqtl are markers
  MAFs <- sample(rbeta(200000,.4,.4),total.loci,replace=F)
  MAFs[which(MAFs > .5)] <- 1- MAFs[which(MAFs > .5)]
  #MAFs <- sample(runif(totalloci,.2,.49),size = totalloci)
  map$MAF <- sample(MAFs,length(map$MAF),replace=F) # Assign minor allele frequencies to markers and qtl

  SNPQTL <- sample(all.loci[-markers],num.snpqtl,replace=FALSE)
  if(!is.null(snp.qtl.maf)){ SNPQTL.MAFs <- runif(num.snpqtl,min=snp.qtl.maf[1],max=snp.qtl.maf[2]) }
  map$types[SNPQTL] <- "snpqtl"   # Specify in the map data frame that these loci are snpqtl
  map$MAF[SNPQTL] <- SNPQTL.MAFs   # Assign these loci the specificed minor allele frequencies generated by user


  num.QTL <- total.qtl - num.snpqtl  # vector that contains the number of QTL which will have slightly dominant effects
  rQTL <- sample(all.loci[-c(SNPQTL,markers)],num.QTL,replace=F)
  map$types[rQTL] <- "qtl"        # Specify in the map data frame that these loci are qtl

  if(!is.null(marker.maf)){
    map$MAF[which(map$types == "m")] <- .49}

  #Now save map output:
  datevalue <- date()
  datevector <- unlist(strsplit(datevalue,"\\s"))
  timevector <- unlist(strsplit(datevector[4],":"))
  newtime <- paste(timevector[1],timevector[2],sep="h")
  newdate <- paste(datevector[3],datevector[2],datevector[5],sep="")
  namestem <- paste(newdate,newtime,sep="_")

  if(save){
    mapname<-paste(namestem,"_map.txt")
    write.table(map,file=mapname, quote=F, row.names=F, col.names=T, sep="\t")}

  #mapinfo<-list(genetic.map=map, total.loci.num=total.loci,total.qtl.num=total.qtl, total.SNPQTL.num=num.snpqtl,
  #              QTLSNP.loci=sort(SNPQTL),date.time=namestem, rQTL.loci=sort(rQTL), available.Markers=markers,last.locus.per.chrom=chromo.loci.index)
  #cat("The returned object is a list containing:\n")
  #cat("Genetic map = a data.frame of 7 columns: chr, locus names, intervals, types (marker, snp, or qtl), MAFs, position, and recombination fractions;\n")
  return(map)
  # End of function #
}

