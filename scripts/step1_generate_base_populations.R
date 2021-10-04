# Run this on Domino ####
source('/mnt/simulator/create_map.R')
source('/mnt/simulator/create_parents.R')
source('/mnt/simulator/calc_TGV.R')
source('/mnt/simulator/sim_phenos.R')
library(parallel);set.seed(448234)
seed.vector <- sample(1235:4536345,size = 300,replace =F) 
rm(.Random.seed)
mclapply(1:300,function(pop){
  set.seed(seed.vector[pop])
  genetic.map <- create_genetic_map(num.chromos = 12,map.length = 1800,num.markers = 120,total.qtl = 2640,num.snpqtl = 1960,distribute.loci = "even",marker.distribution = "equally-spaced",snp.qtl.maf = c(0.01,0.02))
  parents <- create_parents(map.info = genetic.map,num.parents = 280,max.delt.allele = 14,par.markers.unique = T,heterozygous.markers = T)
  parents.TGV <- calc_TGV(geno.info = parents,map.info = genetic.map,A = 1,a = -100,dom.coeff = 0,founder = T)
  parents.phenos <- sim_phenos(TGV.object = parents.TGV,h2 = .3)
  founder.h2 <- var(parents.phenos$genetic.values)/var(parents.phenos$phenos)
  E.sd <- sqrt((var(parents.phenos$genetic.values)/founder.h2) - var(parents.phenos$genetic.values))
  
  this.pop <- paste0("/mnt/data/base_populations/","base_population_",pop,".RData")
  base_pop_data <- list(
   "genetic.map" = genetic.map,
    "parents" = parents,
    "parents.TGV"=parents.TGV,
    "parents.phenos"=parents.phenos,
    "founder.h2"=founder.h2,
    "E.sd"=E.sd)
  save(base_pop_data,file = this.pop,compress = T)
  },mc.cores=30)

