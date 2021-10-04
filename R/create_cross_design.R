#' Simulate the creation of a cross design object
#'
#' This function creates the cross design object which can be generated from parents or the previous generation selections
#' @param parent.info The parent object from which contains individuals to make crosses from. This should come from either the create_parents() or extract_selections()
#' @param prog.percross Number of progeny that should be made per cross
#' @param generation Number that specifies which generation
#' @param mating.design Specify the type of mating design to be used in making crosses: "AssortiveMating" (Default), "Self", "RandomMating","SP","MP", or "cross.file.input". (See details for description)
#' @param coancest.thresh logical. Should a coancestry threshold be used as to minimize making related crosses?
#' @param use.par.marker.effects logical.  Set to TRUE if using the parents marker effects
#' @param use.op.par.phenos logical. Set to TRUE if passing parents created from OP_crossing()
#' @param cross.file User provided cross.file; Should be used if mating.design = "cross.design.input"
#' @param save logical. Saves the output of genetic map (Default: FALSE)
#' @param parent.phenos Object returned from sim_parent_phenos()
#' @keywords create cross design
#' @export
#' @examples
#' cross.file  <- create_cross_design(mapinfo = the.map,parent.info = op.families,prog.percross = 60,
#' generation = 1, use.op.par.phenos = T)


####Create Cross Design & Pedigree's####
create_cross_design <- function(parent.info, prog.percross, generation,
                                mating.design=mt.design, coancest.thresh=T,
                                use.par.marker.effects = F, use.op.par.phenos = F,
                                cross.file= NULL, save=F, parent.phenos = NULL){
  library(reshape2)
  
  coancest=NULL
  
  if(length(dim(parent.info$genos.3d)) == 2){num.parents = 1; parent.names <- unique(parent.info$select.ped.ids)} else {
    num.parents=ncol(parent.info$genos.3d[,,1])
    parent.names <- colnames(parent.info$genos.3d)}
  
  if(mating.design == "Self"){
    par1 <- parent.names
    par2 <- parent.names
    prog.par <- rep(prog.percross,length(par1))
    cross.design <- cbind(par1,par2,prog.par)}
  
  #Must be divisible by 4
  if(mating.design == "AssortiveMating"){
    total <- num.parents/4
    if(generation == 1){
      if(use.par.marker.effects | use.op.par.phenos) {sorted.ebvs <- parent.info$mean.parent.phenos
      } else {
        if(length(parent.info$phenos) > 0) {sorted.ebvs <- parent.info$selection.phenos
        } else {sorted.ebvs <- sort(parent.phenos$phenos,decreasing=T)}}
      
      relationship.matrix <- as.data.frame(diag(1,nrow=num.parents,ncol=num.parents))
      colnames(relationship.matrix) <- names(sorted.ebvs); rownames(relationship.matrix) <- names(sorted.ebvs)
      
      top.25 <- names(sorted.ebvs)[1:total]
      top.50 <- names(sorted.ebvs)[(total+1):(total*2)]
      bottom.50 <- names(sorted.ebvs)[(total*2+1):(total*3)]
      bottom.25 <- names(sorted.ebvs)[(total*3+1):(total*4)]
    } else {
      sorted.ebvs <- sort(parent.info$select.EBVs,decreasing = T)
      relationship.matrix <- parent.info$relmat[which(colnames(parent.info$relmat) %in% names(sorted.ebvs)),which(colnames(parent.info$relmat) %in% names(sorted.ebvs))]
      top.25 <- names(sorted.ebvs)[1:total]
      top.50 <- names(sorted.ebvs)[(total+1):(total*2)]
      bottom.50 <- names(sorted.ebvs)[(total*2+1):(total*3)]
      bottom.25 <- names(sorted.ebvs)[(total*3+1):(total*4)]
    }
    
    all <- vector("list")
    for(i in 1:total){all[[i]] <- top.25[-i]}
    
    if(coancest.thresh){
      E <- upper.tri(relationship.matrix); relationship.matrix[E] <- NA
      coancestry <- seq(.01,1,.005)
      for(i in 1:length(coancestry)){
        new <- subset(reshape2::melt(as.matrix(t(relationship.matrix))), value < coancestry[i])
        matches1 <- match(new[,1],names(sorted.ebvs))
        matches2 <- match(new[,2],names(sorted.ebvs))
        a <- sorted.ebvs[matches1]
        b <- sorted.ebvs[matches2]
        new$mean.bv <- (a+b)/2
        new <- new[order(new$mean.bv,decreasing=T),]
        coancest <- coancestry[i]
        cross.design <- matrix(NA,nrow=1,ncol=2)
        for(each in 1:length(top.25)){
          par <- top.25[each]
          options1 <- which(new[,1] %in% par)
          options2 <- which(new[,2] %in% par)
          options1.other <- options1[which(as.character(new[options1,2]) %in% all[[each]])]
          options2.other <- options2[which(new[options2,1] %in% all[[each]])]
          options <- sort(c(options1.other,options2.other))
          this.cross <- new[options[1:2],c(1:2)] ; colnames(this.cross) <- c("V1","V2")
          cross.design <- rbind(cross.design,this.cross)
          if(each == 1){cross.design <- cross.design[-1,]}
          new <- new[-c(which(new[,1] %in% this.cross[,1] & new[,2] %in% this.cross[,2] )),]
        }
        
        top50.samp <-top.50
        for(each in 1:length(top.25)){
          par <- top.25[each]
          options1 <-   which(new[,1] %in% par)
          options2 <-  which( new[,2] %in% par )
          options1.other <- options1[which(as.character(new[options1,2]) %in% top50.samp)]
          options2.other <- options2[which(new[options2,1] %in% top50.samp)]
          options <- sort(c(options1.other,options2.other))
          this.cross <- new[options[1],c(1:2)] ; colnames(this.cross) <- c("V1","V2")
          cross.design <- rbind(cross.design,this.cross)
          new <- new[-c(which(new[,1] %in% this.cross[,1] & new[,2] %in% this.cross[,2] )),]
          top50.samp <- top50.samp[-c(which(top50.samp %in% this.cross[,1] | top50.samp %in% this.cross[,2]))]
        }
        
        parent1.2.top25 <- sample(top.25,length(top.25),replace=F)
        parent2.2.bottom50 <- sample(bottom.50,length(top.25),replace=F)
        parent1.3.top50 <- sample(top.50,length(top.25),replace=F)
        parent2.3.bottom25 <- sample(bottom.25,length(top.25),replace=F)
        parent1.4.top50 <- sample(top.50,length(top.25),replace=F)
        parent2.4.bottom50 <- sample(bottom.50,length(top.25),replace=F)
        
        all.other.crosses1 <- c(parent1.2.top25,parent1.3.top50,parent1.4.top50)
        all.other.crosses2 <- c(parent2.2.bottom50,parent2.3.bottom25,parent2.4.bottom50)
        all.other.crosses <- cbind(all.other.crosses1,all.other.crosses2)
        
        colnames(all.other.crosses) <- c("V1","V2")
        c <- rbind(cross.design,all.other.crosses)
        if(anyNA(c)){cross.design <- NA} else {cross.design <- c}
        if(anyNA(cross.design)){} else { break}
      }
      cross.design[,3] <- rep(prog.percross,nrow(cross.design))
      coancest <- coancestry[i]
    } else {
      parent2 <- vector()
      for(each in 1:length(top.25)){
        par1 <- each
        par2 <- sample(all[[each]],2,replace=F)
        parent2 <- c(parent2,par2)
        for(i in 1:length(par2)){
          all[[par2[i]]] <- all[[par2[i]]][-each]
        }}
      parent1.top25 <- rep(top.25,each=2)
      parent2.top25 <- parent2
      parent1.1.top25 <- sample(top.25,length(top.25),replace=F)
      parent2.1.top50 <- sample(top.50,length(top.25),replace=F)
      parent1.2.top25 <- sample(top.25,length(top.25),replace=F)
      parent2.2.bottom50 <- sample(bottom.50,length(top.25),replace=F)
      parent1.3.top50 <- sample(top.50,length(top.25),replace=F)
      parent2.3.bottom25 <- sample(bottom.25,length(top.25),replace=F)
      parent1.4.top50 <- sample(top.50,length(top.25),replace=F)
      parent2.4.bottom50 <- sample(bottom.50,length(top.25),replace=F)
      
      pars1 <- c(parent1.top25,parent1.1.top25,parent1.2.top25,parent1.3.top50,parent1.4.top50)
      
      pars2 <- c(parent2.top25,parent2.1.top50,parent2.2.bottom50,parent2.3.bottom25,parent2.4.bottom50)
      
      
      cross.design <- cbind(par1=pars1,
                            par2=pars2,
                            prog=rep(prog.percross,length(pars1)))}}
  
  if(mating.design == "RandomMating"){
    cross.design= matrix(sample(colnames(parent.info$genos.3d),replace=FALSE, size=num.parents),nrow = num.parents/2,ncol=2)
    cross.design <- cbind(cross.design,as.numeric(rep(prog.percross,length(cross.design[,1]))))}
  
  if(mating.design == "SP" ){
    if(coancest.thresh){
      par1 <- vector(); par2 <- vector()
      relationship.matrix <- parent.info$relmat[which(colnames(parent.info$relmat) %in% parent.info$select.ped.ids),which(colnames(parent.info$relmat) %in% parent.info$select.ped.ids)]
      E <- upper.tri(relationship.matrix); relationship.matrix[E] <- NA
      coancestry <- seq(.01,1,.005)
      for(i in 1:length(coancestry)){
        new <- subset(melt(as.matrix(t(relationship.matrix))), value < coancestry[i])
        matches1 <- match(new[,1],names(parent.info$select.EBVs))
        matches2 <- match(new[,2],names(parent.info$select.EBVs))
        a <- parent.info$select.EBVs[matches1]
        b <- parent.info$select.EBVs[matches2]
        new$mean.bv <- (a+b)/2
        new <- new[order(new$mean.bv,decreasing=T),]
        coancest <- coancestry[i]
        cross.design <- matrix(NA,nrow=(num.parents/2),ncol=3)
        for(i in 1:(num.parents/2)) {
          cross.design[i,1] <-  new[1,1]
          cross.design[i,2] <-  new[1,2]
          cross.design[i,3] <- prog.percross
          new <- new[which(new[,1] != cross.design[i,1] & new[,1] != cross.design[i,2]),]
          new <- new[which(new[,2] != cross.design[i,1] & new[,2] != cross.design[i,2]),]}
        if (!anyNA(cross.design)) {break}
      }} else {
        if(generation==1){
          if(use.par.marker.effects | use.op.par.phenos){sorted.ebvs <- parent.info$pars
          }} else {
            sorted.ebvs <- sort(parent.info$select.EBVs,decreasing = T)}
        n.s.ebvs <- length(sorted.ebvs)
        parent1 <- names(sorted.ebvs)[c(seq(from = 1,to = n.s.ebvs,by = 2))]
        parent2 <- names(sorted.ebvs)[c(seq(2,n.s.ebvs,2))]
        cross.design <- cbind(par1=parent1,par2=parent2,prog=rep(prog.percross,length(parent1)))}}
  
  if(mating.design == "MP") {
    if(coancest.thresh){
      par1 <- vector()
      par2 <- vector()
      relationship.matrix <- parent.info$relmat[which(colnames(parent.info$relmat) %in% parent.info$select.ped.ids),which(colnames(parent.info$relmat) %in% parent.info$select.ped.ids)]
      D <- relationship.matrix
      coancestry <- seq(.01,1,.005)
      for(i in 1:length(coancestry)){
        D[D < coancestry[i]] <- 0
        D[D > coancestry[i]] <- 1
        D <- 1-D
        minimum <- colSums(D, na.rm=T)
        coancest <- coancestry[i]
        D <- relationship.matrix
        #if(min(minimum) >= 90) break }
        sorted.minimum <- sort(minimum, decreasing=F)
        relationship.matrix <- relationship.matrix[,names(sorted.minimum)]
        E <- upper.tri(relationship.matrix)
        relationship.matrix[E] <- NA
        new <- subset(melt(as.matrix(t(relationship.matrix))), value < coancest)
        matches1 <- match(new[,1],names(parent.info$select.EBVs))
        matches2 <- match(new[,2],names(parent.info$select.EBVs))
        a <- parent.info$select.EBVs[matches1]
        b <- parent.info$select.EBVs[matches2]
        new$mean.bv <- (a+b)/2
        sib1 <- sample(colnames(relationship.matrix)[seq(1,num.parents,2)],num.parents/2)
        sib2 <- colnames(relationship.matrix)[seq(2,num.parents,2)]
        test <- new[order(-new[,4]),]
        par1 <- test[1:(length(parent.info$selections)/2),1]
        par2 <- test[1:(length(parent.info$selections)/2),2]
        progeny <- rep(prog.percross, length(par1))
        cross.design=matrix(data=c(par1,par2,progeny),nrow=length(progeny),ncol=3)
        if (!anyNA(cross.design)) {break}}} else {
          if(generation==1){
            if(use.par.marker.effects | use.op.par.phenos){sorted.ebvs <- parent.info$pars
            }} else {
              sorted.ebvs <- sort(parent.info$select.EBVs,decreasing = T)}
          n.s.ebvs <- as.numeric(length(sorted.ebvs))
          parent1 <- rep(names(sorted.ebvs[1]),n.s.ebvs-1)
          parent2 <- names(sorted.ebvs[-1])
          parent3 <- rep(names(sorted.ebvs[2]),(n.s.ebvs/2)+1)
          parent4 <- names(sorted.ebvs[3:((n.s.ebvs/2)+3)])
          cross.design <- cbind(par1=c(parent1,parent3),par2=c(parent2,parent4),prog=rep(prog.percross,length(c(parent1,parent3))))
          
        }}
  
  if(mating.design == "cross.file.input") {
    colClasses = c("integer", "integer", "integer")
    cross.design = as.matrix(read.table(cross.file,header=F,stringsAsFactors = F))}
  
  num.crosses <- nrow(cross.design)
  total.progeny <- sum(as.numeric(cross.design[,3]))
  par1.list<-vector(); par2.list<-vector()
  par1.id<-vector(); par2.id<-vector()
  the.generation<- rep(generation,total.progeny)
  
  if (generation==1){
    total.indiv   <- total.progeny+length(unique(c(cross.design[,1],cross.design[,2])))
    cumul.total <- length(unique(c(cross.design[,1],cross.design[,2]))) + total.progeny
    indiv.IDs<-c(seq(max(as.numeric(parent.names))+1,max(as.numeric(parent.names))+(total.progeny)))
    for (m in 1:num.crosses){
      cross.prog<-as.numeric(cross.design[m,3])
      par1<-as.character(cross.design[m,1])
      par2<-as.character(cross.design[m,2])
      par1.list<-c(par1.list,rep(par1,cross.prog))
      par2.list<-c(par2.list,rep(par2,cross.prog))
      }
    
    uni <- length(unique(c(cross.design[,1],cross.design[,2])))
    parent.pedigree <- cbind("ID"=rep(unique(c(cross.design[,1],cross.design[,2])),1),"Par1" = rep(0,uni), "Par2"= rep(0,uni),"gener"=rep(0,uni))
    progeny.pedigree<-cbind("ID"=indiv.IDs[1:length(indiv.IDs)], "Par1"=par1.list[1:length(par1.list)],"Par2"=par2.list[1:length(par1.list)],"gener"=the.generation[1:length(the.generation)])
    full.pedigree <- rbind(parent.pedigree,progeny.pedigree)
    selection.pedigree <- full.pedigree
    
  } else {
    indiv.IDs<- c("indiv.IDs",seq(max(as.numeric(parent.info$fullped[,1]))+1,max(as.numeric(parent.info$fullped[,1]))+total.progeny))
    cumul.total<- parent.info$cumulative.total + total.progeny
    vector2 <- parent.info$select.ped.ids
    for (m in 1:num.crosses){
      cross.prog<-as.numeric(cross.design[m,3])
      par1<-c(cross.design[m,1])
      par2<-c(cross.design[m,2])
      #par1.id<- (vector2[par1])
      #par2.id<- (vector2[par2])
      par1.list<-c(par1.list,rep(as.vector(par1),cross.prog))
      par2.list<-c(par2.list,rep(as.vector(par2),cross.prog))}
    progeny.pedigree<-cbind("ID"=indiv.IDs[2:length(indiv.IDs)], "Par1"=par1.list,"Par2"=par2.list,"gener"=the.generation)
    full.pedigree <- rbind(parent.info$fullped,progeny.pedigree)
    selection.pedigree <- rbind(parent.info$ped,progeny.pedigree)
  }
  
  #if(save) {
  #  pedfilename=paste(rep.num,mapinfo$date.time,prefix,"-pedigree.txt",sep="")
  #  write.table(progeny.pedigree,pedfilename, quote = F, row.names = F, col.names = T, sep=" ")}
  
  ###Imp that genos.3d is first in list item, and total progeny is next list item bc in createtgv routine it expects genos.3d as first element and total prog as second
  cross.design <-list(cross.design=cross.design, total.progeny.number=total.progeny, progeny.pedigree=progeny.pedigree, full.pedigree= full.pedigree,
                      num.parents=parent.info$num.parents, coancestry.threshold=coancest, selection.ped=selection.pedigree,
                      num.crosses=num.crosses,cross.prog=cross.prog,cumul.total=cumul.total, parent.IDs=progeny.pedigree[,1])
  return(cross.design)}
