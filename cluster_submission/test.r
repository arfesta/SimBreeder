#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Run this on Domino ####
source('/share/rosswhet/arfesta/SimBreeder/R/sim_phenos.R')
source('/share/rosswhet/arfesta/SimBreeder/R/create_cross_design.R')
source('/share/rosswhet/arfesta/SimBreeder/R/make_crosses.R')
source('/share/rosswhet/arfesta/SimBreeder/R/calc_TGV.R')
source('/share/rosswhet/arfesta/SimBreeder/R/extract_selections.new3.R')

## Run Values ####
NumParents=64
prog.per.cross = 60
af.selection = "GBLUP"
rel.mat.cross = "markers"
num.sel.af = 64
num.of.cores = 1

#Values for alleles and trait variation
Major.value =1
Minor.value =-100
Dominance.Coeff=0
indiv.tree.h2 = .3

load("/share/rosswhet/arfesta/SimBreeder/cluster_submission/all_seeds.RData")
the_seeds <- all_seeds[args[1]]

  load.pop <- paste0("/share/rosswhet/arfesta/SimBreeder/base_population/base_population_op_test_",args[1],".RData")
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
set.seed(the_seeds[1]); progeny1.phenos <- sim_phenos(TGV.object = progeny1.TGV,h2 = indiv.tree.h2)
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

####Generation 2#####
cross.design <- create_cross_design(parent.info = progeny1.extractions,prog.percross = prog.per.cross,generation  = 2)
progeny2 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny1.extractions,num.cores = num.of.cores)
progeny2.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny2,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
set.seed(the_seeds[2]); progeny2.phenos <- sim_phenos(TGV.object = progeny2.TGV,h2 = indiv.tree.h2)
progeny2.extractions <- extract_selections(relationship.matrix.type = rel.mat.cross, map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny1.TGV,
                                           past.phenos = progeny1.phenos,parent.info = progeny1.extractions,progeny.info = progeny2,
                                           progeny.TGV = progeny2.TGV,progeny.phenos = progeny2.phenos,among.family.selection = af.selection,
                                           num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)


genetic.gain.mine[1,2] <- mean(progeny2.extractions$select.genval) - mean(op.test$genetic.values)
genotypic.variance[1,2]<- var(progeny2.extractions$select.genval)
phenotypic.gain.mine[1,2] <-mean(progeny2.extractions$selection.phenos) - mean(op.test$phenos)
mean.select.inbreeding[1,2]<- mean(progeny2.extractions$select.inbred.level)
mean.pop.inbreeding[1,2]<- mean(progeny2.extractions$prog.inbred.level)
bulmer.effect[1,2]<- progeny2.extractions$bulmer.effect
coancest.threshold[1,2] <- cross.design$coancestry.threshold
delt.alleles[1,2] <- sum(progeny2.extractions$delt.alleles)/NumParents
max.delt.alleles[1,2] <- max(progeny2.extractions$delt.allele)
min.delt.alleles[1,2] <- min(progeny2.extractions$delt.allele)

rm(progeny1.TGV,progeny1.extractions, progeny1.phenos,progeny2);gc(verbose = F,full = T)


####Generations 3####
cross.design <- create_cross_design(parent.info = progeny2.extractions,prog.percross = prog.per.cross,generation  = 3)
progeny3 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny2.extractions,num.cores = num.of.cores)
progeny3.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny3,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
set.seed(the_seeds[3]); progeny3.phenos <- sim_phenos(TGV.object = progeny3.TGV,h2 = indiv.tree.h2)
progeny3.extractions <- extract_selections(relationship.matrix.type = rel.mat.cross, map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny2.TGV,
                                           past.phenos = progeny2.phenos,parent.info = progeny2.extractions,progeny.info = progeny3,
                                           progeny.TGV = progeny3.TGV,progeny.phenos = progeny3.phenos,among.family.selection = af.selection,
                                           num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)


genetic.gain.mine[1,3] <- mean(progeny3.extractions$select.genval) - mean(op.test$genetic.values)
genotypic.variance[1,3]<- var(progeny3.extractions$select.genval)
phenotypic.gain.mine[1,3] <-mean(progeny3.extractions$selection.phenos) - mean(op.test$phenos)
mean.select.inbreeding[1,3]<- mean(progeny3.extractions$select.inbred.level)
mean.pop.inbreeding[1,3]<- mean(progeny3.extractions$prog.inbred.level)
bulmer.effect[1,3]<- progeny3.extractions$bulmer.effect
coancest.threshold[1,3] <- cross.design$coancestry.threshold
delt.alleles[1,3] <- sum(progeny3.extractions$delt.alleles)/NumParents
max.delt.alleles[1,3] <- max(progeny3.extractions$delt.allele)
min.delt.alleles[1,3] <- min(progeny3.extractions$delt.allele)

rm(progeny2.TGV,progeny2.extractions, progeny2.phenos,progeny3);gc(verbose = F,full = T)

####Generations 4####
cross.design <- create_cross_design(parent.info = progeny3.extractions,prog.percross = prog.per.cross,generation  = 4)
progeny4 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny3.extractions,num.cores = num.of.cores)
progeny4.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny4,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
set.seed(the_seeds[4]); progeny4.phenos <- sim_phenos(TGV.object = progeny4.TGV,h2 = indiv.tree.h2)
progeny4.extractions <- extract_selections(relationship.matrix.type = rel.mat.cross, map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny3.TGV,
                                           past.phenos = progeny3.phenos,parent.info = progeny3.extractions,progeny.info = progeny4,
                                           progeny.TGV = progeny4.TGV,progeny.phenos = progeny4.phenos,among.family.selection = af.selection,
                                           num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)

genetic.gain.mine[1,4] <- mean(progeny4.extractions$select.genval) - mean(op.test$genetic.values)
genotypic.variance[1,4]<- var(progeny4.extractions$select.genval)
phenotypic.gain.mine[1,4] <-mean(progeny4.extractions$selection.phenos) - mean(op.test$phenos)
mean.select.inbreeding[1,4]<- mean(progeny4.extractions$select.inbred.level)
mean.pop.inbreeding[1,4]<- mean(progeny4.extractions$prog.inbred.level)
bulmer.effect[1,4]<- progeny4.extractions$bulmer.effect
coancest.threshold[1,4] <- cross.design$coancestry.threshold
delt.alleles[1,4] <- sum(progeny4.extractions$delt.alleles)/NumParents
max.delt.alleles[1,4] <- max(progeny4.extractions$delt.allele)
min.delt.alleles[1,4] <- min(progeny4.extractions$delt.allele)

rm(progeny3.TGV,progeny3.extractions, progeny3.phenos,progeny4);gc(verbose = F,full = T)


####Generations 5####
cross.design <- create_cross_design(parent.info = progeny4.extractions,prog.percross = prog.per.cross,generation  = 5)
progeny5 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny4.extractions,num.cores = num.of.cores)
progeny5.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny5,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
set.seed(the_seeds[5]); progeny5.phenos <- sim_phenos(TGV.object = progeny5.TGV,h2 = indiv.tree.h2)
progeny5.extractions <- extract_selections(relationship.matrix.type = rel.mat.cross, map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny4.TGV,
                                           past.phenos = progeny4.phenos,parent.info = progeny4.extractions,progeny.info = progeny5,
                                           progeny.TGV = progeny5.TGV,progeny.phenos = progeny5.phenos,among.family.selection = af.selection,
                                           num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)

genetic.gain.mine[1,5] <- mean(progeny5.extractions$select.genval) - mean(op.test$genetic.values)
genotypic.variance[1,5]<- var(progeny5.extractions$select.genval)
phenotypic.gain.mine[1,5] <-mean(progeny5.extractions$selection.phenos) - mean(op.test$phenos)
mean.select.inbreeding[1,5]<- mean(progeny5.extractions$select.inbred.level)
mean.pop.inbreeding[1,5]<- mean(progeny5.extractions$prog.inbred.level)
bulmer.effect[1,5]<- progeny5.extractions$bulmer.effect
coancest.threshold[1,5] <- cross.design$coancestry.threshold
delt.alleles[1,5] <- sum(progeny5.extractions$delt.alleles)/NumParents
max.delt.alleles[1,5] <- max(progeny5.extractions$delt.allele)
min.delt.alleles[1,5] <- min(progeny5.extractions$delt.allele)

rm(progeny4.TGV,progeny4.extractions, progeny4.phenos,progeny5);gc(verbose = F,full = T)

####Generations 6####
cross.design <- create_cross_design(parent.info = progeny5.extractions,prog.percross = prog.per.cross,generation  = 6)
progeny6 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny5.extractions,num.cores = num.of.cores)
progeny6.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny6,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
set.seed(the_seeds[6]); progeny6.phenos <- sim_phenos(TGV.object = progeny6.TGV,h2 = indiv.tree.h2)
progeny6.extractions <- extract_selections(relationship.matrix.type = rel.mat.cross, map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny5.TGV,
                                           past.phenos = progeny5.phenos,parent.info = progeny5.extractions,progeny.info = progeny6,
                                           progeny.TGV = progeny6.TGV,progeny.phenos = progeny6.phenos,among.family.selection = af.selection,
                                           num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)

genetic.gain.mine[1,6] <- mean(progeny6.extractions$select.genval) - mean(op.test$genetic.values)
genotypic.variance[1,6]<- var(progeny6.extractions$select.genval)
phenotypic.gain.mine[1,6] <-mean(progeny6.extractions$selection.phenos) - mean(op.test$phenos)
mean.select.inbreeding[1,6]<- mean(progeny6.extractions$select.inbred.level)
mean.pop.inbreeding[1,6]<- mean(progeny6.extractions$prog.inbred.level)
bulmer.effect[1,6]<- progeny6.extractions$bulmer.effect
coancest.threshold[1,6] <- cross.design$coancestry.threshold
delt.alleles[1,6] <- sum(progeny6.extractions$delt.alleles)/NumParents
max.delt.alleles[1,6] <- max(progeny6.extractions$delt.allele)
min.delt.alleles[1,6] <- min(progeny6.extractions$delt.allele)

rm(progeny5.TGV,progeny5.extractions, progeny5.phenos,progeny6);gc(verbose = F,full = T)

####Generations 7####
cross.design <- create_cross_design(parent.info = progeny6.extractions,prog.percross = prog.per.cross,generation  = 7)
progeny7 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny6.extractions,num.cores = num.of.cores)
progeny7.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny7,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
set.seed(the_seeds[7]); progeny7.phenos <- sim_phenos(TGV.object = progeny7.TGV,h2 = indiv.tree.h2)
progeny7.extractions <- extract_selections(relationship.matrix.type = rel.mat.cross, map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny6.TGV,
                                           past.phenos = progeny6.phenos,parent.info = progeny6.extractions,progeny.info = progeny7,
                                           progeny.TGV = progeny7.TGV,progeny.phenos = progeny7.phenos,among.family.selection = af.selection,
                                           num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)

genetic.gain.mine[1,7] <- mean(progeny7.extractions$select.genval) - mean(op.test$genetic.values)
genotypic.variance[1,7]<- var(progeny7.extractions$select.genval)
phenotypic.gain.mine[1,7] <-mean(progeny7.extractions$selection.phenos) - mean(op.test$phenos)
mean.select.inbreeding[1,7]<- mean(progeny7.extractions$select.inbred.level)
mean.pop.inbreeding[1,7]<- mean(progeny7.extractions$prog.inbred.level)
bulmer.effect[1,7]<- progeny7.extractions$bulmer.effect
coancest.threshold[1,7] <- cross.design$coancestry.threshold
delt.alleles[1,7] <- sum(progeny7.extractions$delt.alleles)/NumParents
max.delt.alleles[1,7] <- max(progeny7.extractions$delt.allele)
min.delt.alleles[1,7] <- min(progeny7.extractions$delt.allele)

rm(progeny6.TGV,progeny6.extractions, progeny6.phenos,progeny7);gc(verbose = F,full = T)

####Generations 8####
cross.design <- create_cross_design(parent.info = progeny7.extractions,prog.percross = prog.per.cross,generation  = 8)
progeny8 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny7.extractions,num.cores = num.of.cores)
progeny8.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny8,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
set.seed(the_seeds[8]); progeny8.phenos <- sim_phenos(TGV.object = progeny8.TGV,h2 = indiv.tree.h2)
progeny8.extractions <- extract_selections(relationship.matrix.type = rel.mat.cross, map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny7.TGV,
                                           past.phenos = progeny7.phenos,parent.info = progeny7.extractions,progeny.info = progeny8,
                                           progeny.TGV = progeny8.TGV,progeny.phenos = progeny8.phenos,among.family.selection = af.selection,
                                           num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)

genetic.gain.mine[1,8] <- mean(progeny8.extractions$select.genval) - mean(op.test$genetic.values)
genotypic.variance[1,8]<- var(progeny8.extractions$select.genval)
phenotypic.gain.mine[1,8] <-mean(progeny8.extractions$selection.phenos) - mean(op.test$phenos)
mean.select.inbreeding[1,8]<- mean(progeny8.extractions$select.inbred.level)
mean.pop.inbreeding[1,8]<- mean(progeny8.extractions$prog.inbred.level)
bulmer.effect[1,8]<- progeny8.extractions$bulmer.effect
coancest.threshold[1,8] <- cross.design$coancestry.threshold
delt.alleles[1,8] <- sum(progeny8.extractions$delt.alleles)/NumParents
max.delt.alleles[1,8] <- max(progeny8.extractions$delt.allele)
min.delt.alleles[1,8] <- min(progeny8.extractions$delt.allele)

rm(progeny7.TGV,progeny7.extractions, progeny7.phenos,progeny8);gc(verbose = F,full = T)

####Generations 9####
cross.design <- create_cross_design(parent.info = progeny8.extractions,prog.percross = prog.per.cross,generation  = 9)
progeny9 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny8.extractions,num.cores = num.of.cores)
progeny9.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny9,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
set.seed(the_seeds[9]); progeny9.phenos <- sim_phenos(TGV.object = progeny9.TGV,h2 = indiv.tree.h2)
progeny9.extractions <- extract_selections(relationship.matrix.type = rel.mat.cross, map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny8.TGV,
                                           past.phenos = progeny8.phenos,parent.info = progeny8.extractions,progeny.info = progeny9,
                                           progeny.TGV = progeny9.TGV,progeny.phenos = progeny9.phenos,among.family.selection = af.selection,
                                           num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)

genetic.gain.mine[1,9] <- mean(progeny9.extractions$select.genval) - mean(op.test$genetic.values)
genotypic.variance[1,9]<- var(progeny9.extractions$select.genval)
phenotypic.gain.mine[1,9] <-mean(progeny9.extractions$selection.phenos) - mean(op.test$phenos)
mean.select.inbreeding[1,9]<- mean(progeny9.extractions$select.inbred.level)
mean.pop.inbreeding[1,9]<- mean(progeny9.extractions$prog.inbred.level)
bulmer.effect[1,9]<- progeny9.extractions$bulmer.effect
coancest.threshold[1,9] <- cross.design$coancestry.threshold
delt.alleles[1,9] <- sum(progeny9.extractions$delt.alleles)/NumParents
max.delt.alleles[1,9] <- max(progeny9.extractions$delt.allele)
min.delt.alleles[1,9] <- min(progeny9.extractions$delt.allele)

rm(progeny8.TGV,progeny8.extractions, progeny8.phenos,progeny9);gc(verbose = F,full = T)

####Generations 10####
cross.design <- create_cross_design(parent.info = progeny9.extractions,prog.percross = prog.per.cross,generation  = 10)
progeny10 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny9.extractions,num.cores = num.of.cores)
progeny10.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny10,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
set.seed(the_seeds[10]); progeny10.phenos <- sim_phenos(TGV.object = progeny10.TGV,h2 = indiv.tree.h2)
progeny10.extractions <- extract_selections(relationship.matrix.type = rel.mat.cross, map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny9.TGV,
                                           past.phenos = progeny9.phenos,parent.info = progeny9.extractions,progeny.info = progeny10,
                                           progeny.TGV = progeny10.TGV,progeny.phenos = progeny10.phenos,among.family.selection = af.selection,
                                           num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)

genetic.gain.mine[1,10] <- mean(progeny10.extractions$select.genval) - mean(op.test$genetic.values)
genotypic.variance[1,10]<- var(progeny10.extractions$select.genval)
phenotypic.gain.mine[1,10] <-mean(progeny10.extractions$selection.phenos) - mean(op.test$phenos)
mean.select.inbreeding[1,10]<- mean(progeny10.extractions$select.inbred.level)
mean.pop.inbreeding[1,10]<- mean(progeny10.extractions$prog.inbred.level)
bulmer.effect[1,10]<- progeny10.extractions$bulmer.effect
coancest.threshold[1,10] <- cross.design$coancestry.threshold
delt.alleles[1,10] <- sum(progeny10.extractions$delt.alleles)/NumParents
max.delt.alleles[1,10] <- max(progeny10.extractions$delt.allele)
min.delt.alleles[1,10] <- min(progeny10.extractions$delt.allele)

rm(progeny9.TGV,progeny9.extractions, progeny9.phenos,progeny10);gc(verbose = F,full = T)

newList <- list("g"=genetic.gain.mine,"gv"=genotypic.variance, "pg"=phenotypic.gain.mine,"si"=mean.select.inbreeding,
                "pi"=mean.pop.inbreeding,"be"=bulmer.effect, "ct"=coancest.threshold, "maxd"=max.delt.alleles,"mind"=min.delt.alleles,
                "da"=delt.alleles)

the.name <- paste("GBLUP.nodom.64.OP.pop_",args[1], sep="")
assign(x=the.name,value=newList)
save(list=the.name,file=paste("/share/rosswhet/arfesta/results/",the.name,".RData",sep=""))
