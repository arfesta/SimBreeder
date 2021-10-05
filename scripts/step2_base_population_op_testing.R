# Now for each base population conduct OP testing ####
source('/mnt/simulator/OP_testing.R')
load("all_seeds.RData")
lapply(1:300,function(pop){
  s <- unlist(all_seeds[pop])[1]+10
  load.pop.name <- paste0("/mnt/data_v2/base_populations/base_population_",pop,".RData")
  load(load.pop.name)
  op.test <- OP_testing(map.info=base_pop_data$genetic.map,parent.info=base_pop_data$parents,
                      parent.phenos=base_pop_data$parents.phenos,parents.TGV=base_pop_data$parents.TGV,
                      cross.prog = 1,dom.coeff = 0,A =1,a = -100,h2 = .3,n.cores= 30,use.seed = T,starting_seed = s)
  
  save.pop.name <- paste0("/mnt/data_v2/base_population_op_testing/base_population_op_test_",pop,".RData")
  save(base_pop_data,op.test,file=save.pop.name,compress=T)
  })
