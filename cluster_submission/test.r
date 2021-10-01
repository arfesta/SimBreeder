#!/usr/bin/env Rscript
args = commandARrgs(trailingOnly=TRUE)

# Run this on Domino ####
source('~/SimBreeder/R/sim_phenos.R')
source('~/SimBreeder/R/create_cross_design.R')
source('~/SimBreeder/R/make_crosses.R')
source('~/SimBreeder/R/calc_TGV.R')
source('~/SimBreeder/R/extract_selections.new3.R')

## Run Values ####
NumParents=64
prog.per.cross = 30
af.selection = "GBLUP"
rel.mat.cross = "markers"
num.sel.af = 64
num.of.cores = 1

#Values for alleles and trait variation
Major.value =1
Minor.value =-100
Dominance.Coeff=0
indiv.tree.h2 = .3

  load.pop <- paste0("~/SimBreeder/base_population/base_population_op_test_",args[1],".RData")
load(load.pop)

###Create matrices to hold outputs####
simulations=1; generations=10
genetic.gain.mine <- matrix(data=NA,nrow=simulations,ncol=generations)
genotypic.variance <- matrix(data=NA,nrow=simulations,ncol=generations)
phenotypic.gain.mine <- matrix(data=NA,nrow=simulations,ncol=generations)
mean.select.inbreeding<- matrix(data=NA,nrow=simulations,ncol=generations)
mean.pop.inbreeding<- matrix(data=NA,nrow=simulations,ncol=generations)
bulmer.effect <- matrix(data=NA,nrow=simulations,ncol=generations)
coancest.threshold <- matrix(data=NA,nrow=simulations,ncol=generations)
delt.alleles <- matrix(data=NA,nrow=simulations,ncol=generations)
max.delt.alleles <- matrix(data=NA,nrow=simulations,ncol=generations)
min.delt.alleles <- matrix(data=NA,nrow=simulations,ncol=generations)


###Subset the open pollinated parent output to use as input for 1st Generation####
op.test$mean.parent.phenos <- op.test$mean.parent.phenos[1:NumParents]
op.test$mean.parent.tgv <- op.test$mean.parent.tgv[1:NumParents]
op.test$delt.alleles <- op.test$delt.alleles[1:NumParents]
op.test$genetic.values <- op.test$genetic.values[1:NumParents]
op.test$phenos <- op.test$phenos[1:NumParents]
op.test$genos.3d <- op.test$genos.3d[,1:NumParents,]
op.test$marker.matrix <- op.test$marker.matrix[1:NumParents,]

#### Create 1st gen progeny#######
cross.design <- create_cross_design(parent.info = op.test,prog.percross = prog.per.cross,generation  = 1,use.op.par.phenos = T)
progeny1 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = base_pop_data$parents,num.cores = num.of.cores)
progeny1.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny1,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
set.seed(NULL); progeny1.phenos <- sim_phenos(TGV.object = progeny1.TGV,h2 = indiv.tree.h2)
progeny1.extractions <- extract_selections(relationship.matrix.type = rel.mat.cross, map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = base_pop_data$parents.TGV,
                                           past.phenos = base_pop_data$parents.phenos,parent.info = op.test,progeny.info = progeny1,
                                           progeny.TGV = progeny1.TGV,progeny.phenos = progeny1.phenos,among.family.selection = af.selection,
                                           num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)                                     

genetic.gain.mine[1,1] <- mean(progeny1.extractions$select.genval) - mean(op.test$genetic.values)
genotypic.variance[1,1] <- var(progeny1.extractions$select.genval)
phenotypic.gain.mine[1,1] <-mean(progeny1.extractions$selection.phenos)-mean(op.test$phenos)
mean.select.inbreeding[1,1] <- mean(progeny1.extractions$select.inbred.level)
mean.pop.inbreeding[1,1]<- mean(progeny1.extractions$prog.inbred.level)
bulmer.effect[1,1] <- progeny1.extractions$bulmer.effect
coancest.threshold[1,1] <- 0
delt.alleles[1,1] <- sum(progeny1.extractions$delt.alleles)/NumParents
max.delt.alleles[1,1] <- max(progeny1.extractions$delt.allele)
min.delt.alleles[1,1] <- min(progeny1.extractions$delt.allele)
rm(progeny1); gc(verbose = F,full = T)

newList <- list("g"=genetic.gain.mine,"gv"=genotypic.variance, "pg"=phenotypic.gain.mine,"si"=mean.select.inbreeding,
                "pi"=mean.pop.inbreeding,"be"=bulmer.effect, "ct"=coancest.threshold, "maxd"=max.delt.alleles,"mind"=min.delt.alleles,
                "da"=delt.alleles)

the.name <- paste("GBLUP.nodom.64.OP.pop_",args[1], sep="")
assign(x=the.name,value=newList)
save(list=the.name,file=paste("~/",the.name,".RData",sep=""))
