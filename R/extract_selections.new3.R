extract_selections <- function(map.info, cross.design, past.tgv, 
                               past.phenos, parent.info, progeny.info, progeny.TGV, progeny.phenos, 
                               selection.strategy, among.family.selection, within.family.selection = NULL, 
                               num.selections.within.family = 1, num.selections.among.family = NULL, 
                               relationship.matrix.type = "pedigree", prefix = NULL, rep.num = NULL, 
                               reduced = F, weighted = F,num.cores = NULL) {
  
  library(MatrixModels)
  library(parallel)
  library(pedigreemm)
  library(data.table)
  library(pedigree)
  library(tibble)
  library(tidyr)
  library(Rfast)
  
  {
    marker.loci <- which(map.info$types == "m")
    snpqtl.loci <- which(map.info$types == "snpqtl")
    current.cross.design <- cross.design$cross.design
    prog.percross <- as.numeric(current.cross.design[1,3])
    prog.pedigree <- cross.design$progeny.pedigree
    generation <- as.numeric(prog.pedigree[1, 4])
    numparents <- cross.design$num.parents
    full.ped <- cross.design$full.pedigree
    selection.ped <- cross.design$selection.ped
    if (generation == 1) {
      if (reduced) { parent.markers <- parent.info$genos.3d[marker.loci, (names(parent.info$mean.parent.phenos)), ]}  else {parent.markers <- past.tgv$markers.matrix}
    } else { parent.markers <- parent.info$all.markers }
    prog.markers <- progeny.TGV$markers.matrix
    prog.phenos <- progeny.phenos$phenos
    prog.genetic.values <- progeny.TGV$genetic.values
    map.markers <- progeny.TGV$marker.map
    colnames(map.markers) <- c("chr", "pos")
    map.markers$chr <- as.character(map.markers$chr)
    if (generation == 1) {
      ped <- pedigree(label = full.ped[, 1], sire = full.ped[,2], dam = full.ped[, 3])
      if (reduced) {
        all.phenos <- c(parent.info$phenos, progeny.phenos)
        all.genetic.vals <- c(parent.info$genetic.values, 
                              prog.genetic.values)
      }
      else {
        all.phenos <- c(past.phenos$phenos, progeny.phenos)
        all.genetic.vals <- c(past.phenos$genetic.values, 
                              prog.genetic.values)
      }
    } else {
      ped <- pedigree(label = selection.ped[, 1], sire = selection.ped[,2], dam = selection.ped[, 3])
      all.phenos <- c(parent.info$all.phenos, prog.phenos)
      all.genetic.vals <- c(past.phenos$genetic.values, prog.genetic.values)
    }
  }
  ### Selection Method ####
  if (among.family.selection == "ABLUP") {
    
    the.data <- as.matrix(getAInv(ped))
    the.data = the.data[i = match(names(prog.phenos), 
                                  colnames(the.data)), j = match(names(prog.phenos), 
                                                                 colnames(the.data))]
    
    s <- ncol(the.data)/8
    e <- s - 1
    first.seq <- seq(1, ncol(the.data), ncol(the.data)/8)
    
    {
      my.env <- new.env()
      my.env$l1 <- the.data[i = c(first.seq[1]:(first.seq[1] + 
                                                  e), first.seq[2]:(first.seq[2] + e)), j = c(first.seq[1]:(first.seq[1] + 
                                                                                                              e), first.seq[2]:(first.seq[2] + e)), drop = T]
      my.env$l2 <- the.data[i = c(first.seq[1]:(first.seq[1] + 
                                                  e), first.seq[3]:(first.seq[3] + e)), j = c(first.seq[1]:(first.seq[1] + 
                                                                                                              e), first.seq[3]:(first.seq[3] + e)), drop = T]
      my.env$l3 <- the.data[i = c(first.seq[1]:(first.seq[1] + 
                                                  e), first.seq[4]:(first.seq[4] + e)), j = c(first.seq[1]:(first.seq[1] + 
                                                                                                              e), first.seq[4]:(first.seq[4] + e)), drop = T]
      my.env$l4 <- the.data[i = c(first.seq[1]:(first.seq[1] + 
                                                  e), first.seq[5]:(first.seq[5] + e)), j = c(first.seq[1]:(first.seq[1] + 
                                                                                                              e), first.seq[5]:(first.seq[5] + e)), drop = T]
      my.env$l5 <- the.data[i = c(first.seq[1]:(first.seq[1] + 
                                                  e), first.seq[6]:(first.seq[6] + e)), j = c(first.seq[1]:(first.seq[1] + 
                                                                                                              e), first.seq[6]:(first.seq[6] + e)), drop = T]
      my.env$l6 <- the.data[i = c(first.seq[1]:(first.seq[1] + 
                                                  e), first.seq[7]:(first.seq[7] + e)), j = c(first.seq[1]:(first.seq[1] + 
                                                                                                              e), first.seq[7]:(first.seq[7] + e)), drop = T]
      my.env$l7 <- the.data[i = c(first.seq[1]:(first.seq[1] + 
                                                  e), first.seq[8]:(first.seq[8] + e)), j = c(first.seq[1]:(first.seq[1] + 
                                                                                                              e), first.seq[8]:(first.seq[8] + e)), drop = T]
      my.env$l8 <- the.data[i = c(first.seq[2]:(first.seq[2] + 
                                                  e), first.seq[3]:(first.seq[3] + e)), j = c(first.seq[2]:(first.seq[2] + 
                                                                                                              e), first.seq[3]:(first.seq[3] + e)), drop = T]
      my.env$l9 <- the.data[i = c(first.seq[2]:(first.seq[2] + 
                                                  e), first.seq[4]:(first.seq[4] + e)), j = c(first.seq[2]:(first.seq[2] + 
                                                                                                              e), first.seq[4]:(first.seq[4] + e)), drop = T]
      my.env$l10 <- the.data[i = c(first.seq[2]:(first.seq[2] + 
                                                   e), first.seq[5]:(first.seq[5] + e)), j = c(first.seq[2]:(first.seq[2] + 
                                                                                                               e), first.seq[5]:(first.seq[5] + e)), drop = T]
      my.env$l11 <- the.data[i = c(first.seq[2]:(first.seq[2] + 
                                                   e), first.seq[6]:(first.seq[6] + e)), j = c(first.seq[2]:(first.seq[2] + 
                                                                                                               e), first.seq[6]:(first.seq[6] + e)), drop = T]
      my.env$l12 <- the.data[i = c(first.seq[2]:(first.seq[2] + 
                                                   e), first.seq[7]:(first.seq[7] + e)), j = c(first.seq[2]:(first.seq[2] + 
                                                                                                               e), first.seq[7]:(first.seq[7] + e)), drop = T]
      my.env$l13 <- the.data[i = c(first.seq[2]:(first.seq[2] + 
                                                   e), first.seq[8]:(first.seq[8] + e)), j = c(first.seq[2]:(first.seq[2] + 
                                                                                                               e), first.seq[8]:(first.seq[8] + e)), drop = T]
      my.env$l14 <- the.data[i = c(first.seq[3]:(first.seq[3] + 
                                                   e), first.seq[4]:(first.seq[4] + e)), j = c(first.seq[3]:(first.seq[3] + 
                                                                                                               e), first.seq[4]:(first.seq[4] + e)), drop = T]
      my.env$l15 <- the.data[i = c(first.seq[3]:(first.seq[3] + 
                                                   e), first.seq[5]:(first.seq[5] + e)), j = c(first.seq[3]:(first.seq[3] + 
                                                                                                               e), first.seq[5]:(first.seq[5] + e)), drop = T]
      my.env$l16 <- the.data[i = c(first.seq[3]:(first.seq[3] + 
                                                   e), first.seq[6]:(first.seq[6] + e)), j = c(first.seq[3]:(first.seq[3] + 
                                                                                                               e), first.seq[6]:(first.seq[6] + e)), drop = T]
      my.env$l17 <- the.data[i = c(first.seq[3]:(first.seq[3] + 
                                                   e), first.seq[7]:(first.seq[7] + e)), j = c(first.seq[3]:(first.seq[3] + 
                                                                                                               e), first.seq[7]:(first.seq[7] + e)), drop = T]
      my.env$l18 <- the.data[i = c(first.seq[3]:(first.seq[3] + 
                                                   e), first.seq[8]:(first.seq[8] + e)), j = c(first.seq[3]:(first.seq[3] + 
                                                                                                               e), first.seq[8]:(first.seq[8] + e)), drop = T]
      my.env$l19 <- the.data[i = c(first.seq[4]:(first.seq[4] + 
                                                   e), first.seq[5]:(first.seq[5] + e)), j = c(first.seq[4]:(first.seq[4] + 
                                                                                                               e), first.seq[5]:(first.seq[5] + e)), drop = T]
      my.env$l20 <- the.data[i = c(first.seq[4]:(first.seq[4] + 
                                                   e), first.seq[6]:(first.seq[6] + e)), j = c(first.seq[4]:(first.seq[4] + 
                                                                                                               e), first.seq[6]:(first.seq[6] + e)), drop = T]
      my.env$l21 <- the.data[i = c(first.seq[4]:(first.seq[4] + 
                                                   e), first.seq[7]:(first.seq[7] + e)), j = c(first.seq[4]:(first.seq[4] + 
                                                                                                               e), first.seq[7]:(first.seq[7] + e)), drop = T]
      my.env$l22 <- the.data[i = c(first.seq[4]:(first.seq[4] + 
                                                   e), first.seq[8]:(first.seq[8] + e)), j = c(first.seq[4]:(first.seq[4] + 
                                                                                                               e), first.seq[8]:(first.seq[8] + e)), drop = T]
      my.env$l23 <- the.data[i = c(first.seq[5]:(first.seq[5] + 
                                                   e), first.seq[6]:(first.seq[6] + e)), j = c(first.seq[5]:(first.seq[5] + 
                                                                                                               e), first.seq[6]:(first.seq[6] + e)), drop = T]
      my.env$l24 <- the.data[i = c(first.seq[5]:(first.seq[5] + 
                                                   e), first.seq[7]:(first.seq[7] + e)), j = c(first.seq[5]:(first.seq[5] + 
                                                                                                               e), first.seq[7]:(first.seq[7] + e)), drop = T]
      my.env$l25 <- the.data[i = c(first.seq[5]:(first.seq[5] + 
                                                   e), first.seq[8]:(first.seq[8] + e)), j = c(first.seq[5]:(first.seq[5] + 
                                                                                                               e), first.seq[8]:(first.seq[8] + e)), drop = T]
      my.env$l26 <- the.data[i = c(first.seq[6]:(first.seq[6] + 
                                                   e), first.seq[7]:(first.seq[7] + e)), j = c(first.seq[6]:(first.seq[6] + 
                                                                                                               e), first.seq[7]:(first.seq[7] + e)), drop = T]
      my.env$l27 <- the.data[i = c(first.seq[6]:(first.seq[6] + 
                                                   e), first.seq[8]:(first.seq[8] + e)), j = c(first.seq[6]:(first.seq[6] + 
                                                                                                               e), first.seq[8]:(first.seq[8] + e)), drop = T]
      my.env$l28 <- the.data[i = c(first.seq[7]:(first.seq[7] + 
                                                   e), first.seq[8]:(first.seq[8] + e)), j = c(first.seq[7]:(first.seq[7] + 
                                                                                                               e), first.seq[8]:(first.seq[8] + e)), drop = T]
    }
    gc()
    {
      s.phenos1 <- prog.phenos[c(first.seq[1]:(first.seq[1] + 
                                                 e), first.seq[2]:(first.seq[2] + e))]
      s.phenos2 <- prog.phenos[c(first.seq[1]:(first.seq[1] + 
                                                 e), first.seq[3]:(first.seq[3] + e))]
      s.phenos3 <- prog.phenos[c(first.seq[1]:(first.seq[1] + 
                                                 e), first.seq[4]:(first.seq[4] + e))]
      s.phenos4 <- prog.phenos[c(first.seq[1]:(first.seq[1] + 
                                                 e), first.seq[5]:(first.seq[5] + e))]
      s.phenos5 <- prog.phenos[c(first.seq[1]:(first.seq[1] + 
                                                 e), first.seq[6]:(first.seq[6] + e))]
      s.phenos6 <- prog.phenos[c(first.seq[1]:(first.seq[1] + 
                                                 e), first.seq[7]:(first.seq[7] + e))]
      s.phenos7 <- prog.phenos[c(first.seq[1]:(first.seq[1] + 
                                                 e), first.seq[8]:(first.seq[8] + e))]
      s.phenos8 <- prog.phenos[c(first.seq[2]:(first.seq[2] + 
                                                 e), first.seq[3]:(first.seq[3] + e))]
      s.phenos9 <- prog.phenos[c(first.seq[2]:(first.seq[2] + 
                                                 e), first.seq[4]:(first.seq[4] + e))]
      s.phenos10 <- prog.phenos[c(first.seq[2]:(first.seq[2] + 
                                                  e), first.seq[5]:(first.seq[5] + e))]
      s.phenos11 <- prog.phenos[c(first.seq[2]:(first.seq[2] + 
                                                  e), first.seq[6]:(first.seq[6] + e))]
      s.phenos12 <- prog.phenos[c(first.seq[2]:(first.seq[2] + 
                                                  e), first.seq[7]:(first.seq[7] + e))]
      s.phenos13 <- prog.phenos[c(first.seq[2]:(first.seq[2] + 
                                                  e), first.seq[8]:(first.seq[8] + e))]
      s.phenos14 <- prog.phenos[c(first.seq[3]:(first.seq[3] + 
                                                  e), first.seq[4]:(first.seq[4] + e))]
      s.phenos15 <- prog.phenos[c(first.seq[3]:(first.seq[3] + 
                                                  e), first.seq[5]:(first.seq[5] + e))]
      s.phenos16 <- prog.phenos[c(first.seq[3]:(first.seq[3] + 
                                                  e), first.seq[6]:(first.seq[6] + e))]
      s.phenos17 <- prog.phenos[c(first.seq[3]:(first.seq[3] + 
                                                  e), first.seq[7]:(first.seq[7] + e))]
      s.phenos18 <- prog.phenos[c(first.seq[3]:(first.seq[3] + 
                                                  e), first.seq[8]:(first.seq[8] + e))]
      s.phenos19 <- prog.phenos[c(first.seq[4]:(first.seq[4] + 
                                                  e), first.seq[5]:(first.seq[5] + e))]
      s.phenos20 <- prog.phenos[c(first.seq[4]:(first.seq[4] + 
                                                  e), first.seq[6]:(first.seq[6] + e))]
      s.phenos21 <- prog.phenos[c(first.seq[4]:(first.seq[4] + 
                                                  e), first.seq[7]:(first.seq[7] + e))]
      s.phenos22 <- prog.phenos[c(first.seq[4]:(first.seq[4] + 
                                                  e), first.seq[8]:(first.seq[8] + e))]
      s.phenos23 <- prog.phenos[c(first.seq[5]:(first.seq[5] + 
                                                  e), first.seq[6]:(first.seq[6] + e))]
      s.phenos24 <- prog.phenos[c(first.seq[5]:(first.seq[5] + 
                                                  e), first.seq[7]:(first.seq[7] + e))]
      s.phenos25 <- prog.phenos[c(first.seq[5]:(first.seq[5] + 
                                                  e), first.seq[8]:(first.seq[8] + e))]
      s.phenos26 <- prog.phenos[c(first.seq[6]:(first.seq[6] + 
                                                  e), first.seq[7]:(first.seq[7] + e))]
      s.phenos27 <- prog.phenos[c(first.seq[6]:(first.seq[6] + 
                                                  e), first.seq[8]:(first.seq[8] + e))]
      s.phenos28 <- prog.phenos[c(first.seq[7]:(first.seq[7] + 
                                                  e), first.seq[8]:(first.seq[8] + e))]
      pheno.list <- list(s.phenos1, s.phenos2, s.phenos3, 
                         s.phenos4, s.phenos5, s.phenos6, s.phenos7, s.phenos8, 
                         s.phenos9, s.phenos10, s.phenos11, s.phenos12, 
                         s.phenos13, s.phenos14, s.phenos15, s.phenos16, 
                         s.phenos17, s.phenos18, s.phenos19, s.phenos20, 
                         s.phenos21, s.phenos22, s.phenos23, s.phenos24, 
                         s.phenos25, s.phenos26, s.phenos27, s.phenos28)
    }
    apply_my_fun <- function(idx, env) {
      eval(parse(text = paste0("result <- env$", ls(env)[idx])))
      the.data <- as.matrix(result)
      n.col <- ncol(the.data)
      h.2 <- var(prog.genetic.values)/var(prog.phenos)
      lambda <- (1 - h.2)/h.2
      I <- diag(n.col)
      s.phenos <- pheno.list[[idx]]
      sol <- solve(rbind(cbind(n.col, t(rep(1, n.col))), 
                         cbind(rep(1, n.col), (I + (lambda * the.data)))), 
                   matrix(c(sum(s.phenos), c(as.vector(s.phenos)))))
      sol <- sol[-1, 1]
      sol
    }
    index <- 1:28
    names(index) <- 1:28
    the.blups <- mclapply(index, apply_my_fun, env = my.env, 
                          mc.cores = num.cores)
    rm(my.env)
    gc()
    {
      SG1 <- apply(rbind(the.blups[[1]][1:(length(s.phenos1)/2)], 
                         the.blups[[2]][1:(length(s.phenos1)/2)], the.blups[[3]][1:(length(s.phenos1)/2)], 
                         the.blups[[4]][1:(length(s.phenos1)/2)], the.blups[[5]][1:(length(s.phenos1)/2)], 
                         the.blups[[6]][1:(length(s.phenos1)/2)], the.blups[[7]][1:(length(s.phenos1)/2)]), 
                   2, mean)
      SG2 <- apply(rbind(the.blups[[1]][((length(s.phenos1)/2) + 
                                           1):length(s.phenos1)], the.blups[[8]][1:(length(s.phenos1)/2)], 
                         the.blups[[9]][1:(length(s.phenos1)/2)], the.blups[[10]][1:(length(s.phenos1)/2)], 
                         the.blups[[11]][1:(length(s.phenos1)/2)], the.blups[[12]][1:(length(s.phenos1)/2)], 
                         the.blups[[13]][1:(length(s.phenos1)/2)]), 2, 
                   mean)
      SG3 <- apply(rbind(the.blups[[2]][((length(s.phenos1)/2) + 
                                           1):length(s.phenos1)], the.blups[[8]][((length(s.phenos1)/2) + 
                                                                                    1):length(s.phenos1)], the.blups[[14]][1:(length(s.phenos1)/2)], 
                         the.blups[[15]][1:(length(s.phenos1)/2)], the.blups[[16]][1:(length(s.phenos1)/2)], 
                         the.blups[[17]][1:(length(s.phenos1)/2)], the.blups[[18]][1:(length(s.phenos1)/2)]), 
                   2, mean)
      SG4 <- apply(rbind(the.blups[[3]][((length(s.phenos1)/2) + 
                                           1):length(s.phenos1)], the.blups[[9]][((length(s.phenos1)/2) + 
                                                                                    1):length(s.phenos1)], the.blups[[14]][((length(s.phenos1)/2) + 
                                                                                                                              1):length(s.phenos1)], the.blups[[19]][1:(length(s.phenos1)/2)], 
                         the.blups[[20]][1:(length(s.phenos1)/2)], the.blups[[21]][1:(length(s.phenos1)/2)], 
                         the.blups[[22]][1:(length(s.phenos1)/2)]), 2, 
                   mean)
      SG5 <- apply(rbind(the.blups[[4]][((length(s.phenos1)/2) + 
                                           1):length(s.phenos1)], the.blups[[10]][((length(s.phenos1)/2) + 
                                                                                     1):length(s.phenos1)], the.blups[[15]][((length(s.phenos1)/2) + 
                                                                                                                               1):length(s.phenos1)], the.blups[[19]][((length(s.phenos1)/2) + 
                                                                                                                                                                         1):length(s.phenos1)], the.blups[[23]][1:(length(s.phenos1)/2)], 
                         the.blups[[24]][1:(length(s.phenos1)/2)], the.blups[[25]][1:(length(s.phenos1)/2)]), 
                   2, mean)
      SG6 <- apply(rbind(the.blups[[5]][((length(s.phenos1)/2) + 
                                           1):length(s.phenos1)], the.blups[[11]][((length(s.phenos1)/2) + 
                                                                                     1):length(s.phenos1)], the.blups[[16]][((length(s.phenos1)/2) + 
                                                                                                                               1):length(s.phenos1)], the.blups[[20]][((length(s.phenos1)/2) + 
                                                                                                                                                                         1):length(s.phenos1)], the.blups[[23]][((length(s.phenos1)/2) + 
                                                                                                                                                                                                                   1):length(s.phenos1)], the.blups[[26]][1:(length(s.phenos1)/2)], 
                         the.blups[[27]][1:(length(s.phenos1)/2)]), 2, 
                   mean)
      SG7 <- apply(rbind(the.blups[[6]][((length(s.phenos1)/2) + 
                                           1):length(s.phenos1)], the.blups[[12]][((length(s.phenos1)/2) + 
                                                                                     1):length(s.phenos1)], the.blups[[17]][((length(s.phenos1)/2) + 
                                                                                                                               1):length(s.phenos1)], the.blups[[21]][((length(s.phenos1)/2) + 
                                                                                                                                                                         1):length(s.phenos1)], the.blups[[24]][((length(s.phenos1)/2) + 
                                                                                                                                                                                                                   1):length(s.phenos1)], the.blups[[26]][((length(s.phenos1)/2) + 
                                                                                                                                                                                                                                                             1):length(s.phenos1)], the.blups[[28]][1:(length(s.phenos1)/2)]), 
                   2, mean)
      SG8 <- apply(rbind(the.blups[[7]][((length(s.phenos1)/2) + 
                                           1):length(s.phenos1)], the.blups[[13]][((length(s.phenos1)/2) + 
                                                                                     1):length(s.phenos1)], the.blups[[18]][((length(s.phenos1)/2) + 
                                                                                                                               1):length(s.phenos1)], the.blups[[22]][((length(s.phenos1)/2) + 
                                                                                                                                                                         1):length(s.phenos1)], the.blups[[25]][((length(s.phenos1)/2) + 
                                                                                                                                                                                                                   1):length(s.phenos1)], the.blups[[27]][((length(s.phenos1)/2) + 
                                                                                                                                                                                                                                                             1):length(s.phenos1)], the.blups[[28]][((length(s.phenos1)/2) + 
                                                                                                                                                                                                                                                                                                       1):length(s.phenos1)]), 2, mean)
    }
    progeny.blups <- c(SG1, SG2, SG3, SG4, SG5, SG6, SG7, 
                       SG8)
    names(progeny.blups) <- names(prog.phenos)
    first <- 1
    last <- prog.percross
    mean.progeny.blups <- vector()
    for (each in 1:nrow(current.cross.design)) {
      mean.progeny.blups <- c(mean.progeny.blups, mean(progeny.blups[first:last]))
      first <- last + 1
      last <- first + prog.percross - 1
    }
    sorted.mean.progeny.blups <- sort(mean.progeny.blups, 
                                      decreasing = T)
    top.280.families <- match(sorted.mean.progeny.blups, 
                              mean.progeny.blups)
    the.selections <- vector()
    first.in.family <- 1
    for (family in 1:length(current.cross.design[, 3])) {
      num.offspring <- as.numeric(current.cross.design[family, 
                                                       3])
      last.in.family <- num.offspring + first.in.family - 
        1
      temp <- (progeny.blups[first.in.family:last.in.family])
      sorted <- sort(temp, decreasing = TRUE)
      best.one <- which(temp == sorted[1:num.selections.within.family])
      selected <- best.one + first.in.family - 1
      the.selections <- c(the.selections, selected)
      first.in.family <- last.in.family + 1
    }
    if (length(num.selections.among.family) > 0) {
      the.selections <- the.selections[top.280.families[1:num.selections.among.family]]
    } else {
      the.selections <- the.selections[top.280.families]
    }
  }
  if (among.family.selection == "GBLUP") {
    prog.markers <- progeny.info$genos.3d[marker.loci, , ]
    gc(full=T,verbose = F,reset = T)
    pt <- proc.time()
    all.m <- rep(0,ncol(prog.markers)*ncol(prog.markers))
    for(all.markers in 1:nrow(prog.markers)){
      unique.markers <-  unique(c(prog.markers[all.markers,,1],prog.markers[all.markers,,2]))
      
      for(each.marker in 1:length(unique.markers)){
        these.in.both <- which(prog.markers[all.markers,,1] %in% unique.markers[each.marker] & prog.markers[all.markers,,2] %in% unique.markers[each.marker])
        if(length(these.in.both) > 1){
          test_in_both <- comb_n(n = these.in.both,k = 2)
          these.in.one <-  which(prog.markers[all.markers,,1] %in% unique.markers[each.marker] | prog.markers[all.markers,,2] %in% unique.markers[each.marker] )
          test_in_one <- comb_n(n = these.in.one,k = 2)
          test_in_one <- test_in_one[,-c(which(test_in_one[1,] %in% these.in.both))]
          idx <- ((test_in_both[1,]-1)*ncol(prog.markers)) + test_in_both[2,]
          all.m[idx] <- all.m[idx] +2
          idx <- ((test_in_one[1,]-1)*ncol(prog.markers)) + test_in_one[2,]
          all.m[idx] <- all.m[idx] +1
        } else{
          these.in.one <-  which(prog.markers[all.markers,,1] %in% unique.markers[each.marker] | prog.markers[all.markers,,2] %in% unique.markers[each.marker] )
          if(length(these.in.one) > 1){
            test_in_one <- comb_n(n = these.in.one,k = 2)
            idx <- ((test_in_one[1,]-1)*ncol(prog.markers)) + test_in_one[2,]
            all.m[idx] <- all.m[idx] +1
          } else {}
        }
      }
    }
    proc.time() - pt
    
    all.m <- all.m/(nrow(prog.markers)*2)
    g.mat <- (matrix(all.m, nrow = ncol(prog.markers), ncol = ncol(prog.markers),byrow = F))
    
    rm(all.m); gc(full=T,reset=T)
    diag(g.mat) <- 1
  
    #pt <- proc.time()
    cc <-  Rfast::transpose(g.mat)[upper.tri( Rfast::transpose(g.mat), diag = F)]
    g.mat <- Rfast::upper_tri.assign(x= g.mat,v = cc, diag = F)
    #proc.time()-pt
    A <- as_tibble(as.matrix(getA(ped)))
    A = A[match(names(prog.phenos), colnames(A)), match(names(prog.phenos),  colnames(A))]
    the.data <- as.matrix(g.mat)*.99 + as.matrix(A) * 0.01
    rm(g.mat,A); gc(full=T,reset = T)
  pt <- proc.time()
    the.data <- solve(the.data)
proc.time() - pt
    n.col <- ncol(the.data)
    h.2 <- var(prog.genetic.values)/var(prog.phenos)
    lambda <- (1 - h.2)/h.2
    I <- diag(n.col)
    DD <- rbind(cbind(n.col, t(rep(1, n.col))), cbind(rep(1, n.col), (I + (lambda * the.data))))
    CC <- matrix(c(sum(prog.phenos),  c(as.vector(prog.phenos))))
pt <- proc.time()
sol <- solve(DD,CC)
proc.time() - pt
    rm(I,the.data,DD,CC); gc(full=T,reset=T)
    sol <- sol[-1, 1]
    gprogeny.blups1 <- unlist(sol)
    names(gprogeny.blups1) <- names(prog.phenos)
    first <- 1
    last <- prog.percross
    mean.progeny.blups <- vector()
    for (each in 1:nrow(current.cross.design)) {
      mean.progeny.blups <- c(mean.progeny.blups, mean(gprogeny.blups1[first:last]))
      first <- last + 1
      last <- first + prog.percross - 1
    }
    sorted.mean.progeny.blups <- sort(mean.progeny.blups, decreasing = T)
    top.280.families <- match(sorted.mean.progeny.blups,  mean.progeny.blups)
    the.selections <- vector()
    first.in.family <- 1
    for (family in 1:length(current.cross.design[, 3])) {
      num.offspring <- as.numeric(current.cross.design[family,  3])
      last.in.family <- num.offspring + first.in.family -  1
      temp <- (gprogeny.blups1[first.in.family:last.in.family])
      sorted <- sort(temp, decreasing = TRUE)
      best.one <- which(temp == sorted[1:num.selections.within.family])
      selected <- best.one + first.in.family - 1
      the.selections <- c(the.selections, selected)
      first.in.family <- last.in.family + 1
    }
    if (length(num.selections.among.family) > 0) {
      the.selections <- the.selections[top.280.families[1:num.selections.among.family]]
    } else {the.selections <- the.selections[top.280.families] }
  }
  if (among.family.selection == "Phenotype") {
    # If using Phenotypes as the among family selection, the progeny blups will be the phenos
    progeny.blups <- prog.phenos
    names(progeny.blups) <- names(prog.phenos) # Names of the progeny.blups will be the names of the progeny phenos
    
    first <- 1
    last <- prog.percross
    mean.progeny.blups <- vector()
    for (each in 1:nrow(current.cross.design)) {
      mean.progeny.blups <- c(mean.progeny.blups, mean(progeny.blups[first:last]))
      first <- last + 1
      last <- first + prog.percross - 1
    }
    sorted.mean.progeny.blups <- sort(mean.progeny.blups,decreasing = T)
    top.280.families <- match(sorted.mean.progeny.blups,mean.progeny.blups)
    the.selections <- vector()
    first.in.family <- 1
    for (family in 1:length(current.cross.design[, 3])) {
      num.offspring <- as.numeric(current.cross.design[family,3])
      last.in.family <- num.offspring + first.in.family -  1
      temp <- (progeny.blups[first.in.family:last.in.family])
      sorted <- sort(temp, decreasing = TRUE)
      best.one <- which(names(temp) == names(sorted[1:num.selections.within.family]))
      selected <- best.one + first.in.family - 1
      the.selections <- c(the.selections, selected)
      first.in.family <- last.in.family + 1
    }
    if (length(num.selections.among.family) > 0) {
      the.selections <- the.selections[top.280.families[1:num.selections.among.family]]
    } else {
      the.selections <- the.selections[top.280.families]
    }
  }
  
  {
    selection.phenos <- prog.phenos[the.selections]
    new.pars.genval <- prog.genetic.values[the.selections]
    if (among.family.selection == "GBLUP") {selection.EBVs <- gprogeny.blups1[the.selections]  
    } else {
      selection.EBVs <- progeny.blups[the.selections]
    }
    sorted.top192 <- sort(selection.phenos, decreasing = T)
    if (generation == 1) {
      if (reduced) {
        all.phenos <- c(parent.info$phenos, selection.phenos)
        all.genetic.vals <- c(parent.info$genetic.values, 
                              new.pars.genval)
      } else {
        all.phenos <- c(past.phenos$phenos, selection.phenos)
        all.genetic.vals <- c(past.phenos$genetic.values, 
                              new.pars.genval)
      }
    } else {
      all.phenos <- c(parent.info$all.phenos, selection.phenos)
      all.genetic.vals <- c(parent.info$all.genetic.vals, new.pars.genval)
    }
    ped <- full.ped[match(names(all.phenos), full.ped[, 1]), ]
    colnames(progeny.info$genos.3d) <- names(prog.phenos)
    new.parent.genos <- progeny.info$genos.3d[, the.selections,]
    numselections <- dim(new.parent.genos)[2]
    select.ids <- as.numeric(names(new.pars.genval)) - numparents
    select.ped.ids <- as.numeric(names(new.pars.genval))
    if (generation == 1) {
      if(length(dim(new.parent.genos)) < 3){
        all.markers1 <- cbind(parent.info$genos.3d[marker.loci,,1], new.parent.genos[marker.loci, 1])
        colnames(all.markers1) <- c(colnames(parent.info$genos.3d[marker.loci,,1]),select.ped.ids)
        all.markers2 <- cbind(parent.info$genos.3d[marker.loci,,2], new.parent.genos[marker.loci, 2])
        colnames(all.markers2) <- c(colnames(parent.info$genos.3d[marker.loci,,2]),select.ped.ids)
        all.markers <- abind(all.markers1, all.markers2, along = 3)
      } else {
        colnames(new.parent.genos) <- select.ped.ids
        all.markers1 <- cbind(parent.info$genos.3d[marker.loci,,1], new.parent.genos[marker.loci, ,1])
        all.markers2 <- cbind(parent.info$genos.3d[marker.loci,, 2], new.parent.genos[marker.loci, ,2])
      }
      all.markers <- abind(all.markers1, all.markers2, along = 3)
    } else {
      colnames(new.parent.genos) <- select.ped.ids
      all.markers1 <- cbind(parent.info$genos.3d[marker.loci,,1], 
                            new.parent.genos[marker.loci, , 1])
      all.markers2 <- cbind(parent.info$genos.3d[marker.loci,,2], 
                            new.parent.genos[marker.loci, , 2])
      all.markers <- abind(all.markers1, all.markers2, along = 3)
    }
    f.ped <- pedigree(label = full.ped[, 1], sire = full.ped[,2], dam = full.ped[, 3])
    pedigree.inbreeding <- inbreeding(f.ped)
    names(pedigree.inbreeding) <- full.ped[, 1]
    progeny.inbreeding <- pedigree.inbreeding[prog.pedigree[,1]]
    selections.inbreeding <- progeny.inbreeding[names(the.selections)]
    progeny.gen.var <- var(prog.genetic.values[names(the.selections)])
    bulmer.effect <- var(prog.genetic.values[names(the.selections)]) -  var(prog.genetic.values)
    if(length(the.selections) == 1){
      out <- length(which(new.parent.genos[snpqtl.loci, 1] == "c"))
      out2 <- length(which(new.parent.genos[snpqtl.loci, 2] == "c"))
      result <- out + out2
    } else { 
      result <- sapply(rep(1:length(the.selections), 1), function(x) {
        out <- length(which(new.parent.genos[snpqtl.loci, , 1][, x] == "c"))
        out2 <- length(which(new.parent.genos[snpqtl.loci, , 2][, x] == "c"))
        outer <- out + out2
      })
      result <- unlist(result) }
  }
  
  
  if (relationship.matrix.type == "pedigree") {
    if (generation == 1) {
      ff <- cbind(label = f.ped@label, sire = f.ped@sire, dam = f.ped@dam)
      ff <- data.frame(ff,stringsAsFactors = F)
      l <- which(ff$label %in% names(selection.phenos))
      rmv <- rep(F,nrow(ff))
      rmv[l] <- T
      A <- makeA(ff,which = rmv)
      A <- read.table("A.txt")
      A <- spread(A,V1,V3)
      A <- A[,-1]
      colnames(A) <- ff$label[l]; rownames(A) <- ff$label[l]
      A[lower.tri(A, diag = F)] <-  t(A)[lower.tri( t(A), diag = F)]
    }  else {
      selection.ped <- cross.design$selection.ped
      ff <- cbind(id = selection.ped[, 1], sire = selection.ped[,2], dam = selection.ped[, 3])
      ff <- data.frame(ff,stringsAsFactors = F)
      l <- which(ff$id %in% names(selection.phenos))
      rmv <- rep(F,nrow(ff))
      rmv[l] <- T
      A <- makeA(ff,which = rmv)
      A <- read.table("A.txt")
      A <- spread(A,V1,V3)
      A <- A[,-1]
      colnames(A) <- ff$id[l]; rownames(A) <- ff$id[l]
      A[lower.tri(A, diag = F)] <-  t(A)[lower.tri( t(A), diag = F)]
    }
    l <- match(names(selection.phenos), rownames(A))
    rel.mat <- A[l, l]
  }
  if (relationship.matrix.type == "markers") {
    allele1 <- all.markers[, , 1]
    allele2 <- all.markers[, , 2]
    allele1 <- allele1[, which(colnames(allele1) %in% names(selection.phenos))]
    allele2 <- allele2[, which(colnames(allele2) %in% names(selection.phenos))]
    all.alleles <- abind(allele1, allele2, along = 3)
    
    
    
    pt <- proc.time()
    all.m <- rep(0,ncol(prog.markers)*ncol(prog.markers))
    for(all.markers in 1:nrow(prog.markers)){
      unique.markers <-  unique(c(prog.markers[all.markers,,1],prog.markers[all.markers,,2]))
      
      for(each.marker in 1:length(unique.markers)){
        these.in.both <- which(prog.markers[all.markers,,1] %in% unique.markers[each.marker] & prog.markers[all.markers,,2] %in% unique.markers[each.marker])
        if(length(these.in.both) > 1){
          test_in_both <- comb_n(n = these.in.both,k = 2)
          these.in.one <-  which(prog.markers[all.markers,,1] %in% unique.markers[each.marker] | prog.markers[all.markers,,2] %in% unique.markers[each.marker] )
          test_in_one <- comb_n(n = these.in.one,k = 2)
          test_in_one <- test_in_one[,-c(which(test_in_one[1,] %in% these.in.both))]
          idx <- ((test_in_both[1,]-1)*ncol(prog.markers)) + test_in_both[2,]
          all.m[idx] <- all.m[idx] +2
          idx <- ((test_in_one[1,]-1)*ncol(prog.markers)) + test_in_one[2,]
          all.m[idx] <- all.m[idx] +1
        } else{
          these.in.one <-  which(prog.markers[all.markers,,1] %in% unique.markers[each.marker] | prog.markers[all.markers,,2] %in% unique.markers[each.marker] )
          if(length(these.in.one) > 1){
            test_in_one <- comb_n(n = these.in.one,k = 2)
            idx <- ((test_in_one[1,]-1)*ncol(prog.markers)) + test_in_one[2,]
            all.m[idx] <- all.m[idx] +1
          } else {}
        }
      }
    }
    proc.time() - pt
    all.m <- all.m/240
    g.mat <- matrix(all.m, nrow = ncol(all.alleles), ncol = ncol(all.alleles),byrow = F)
    rm(all.m);gc(full=T,reset=T)
    
    diag(g.mat) <- 1
    cc <-  Rfast::transpose(g.mat)[upper.tri( Rfast::transpose(g.mat), diag = F)]
    g.mat <- Rfast::upper_tri.assign(x= g.mat,v = cc, diag = F)
    rownames(g.mat) <- colnames(allele1)
    colnames(g.mat) <- colnames(allele1)
    l <- match(names(selection.phenos), rownames(g.mat))
    rel.mat <- g.mat[l, l]
    rm(g.mat);gc(full=T)
    
  }
  extraction.info <- list(relmat = rel.mat, delt.alleles = result, 
                          selections = the.selections, bulmer.effect = bulmer.effect, 
                          select.EBVs = selection.EBVs, all.genetic.vals = all.genetic.vals, 
                          selection.phenos = selection.phenos, ped = ped, prog.inbred.level = progeny.inbreeding, 
                          select.inbred.level = selections.inbreeding, genos.3d = new.parent.genos, 
                          num.parents = numselections, select.genval = new.pars.genval, 
                          fullped = full.ped, par.ids = select.ids, select.ped.ids = select.ped.ids, 
                          all.markers = all.markers, all.phenos = all.phenos, cumulative.total = cross.design$cumul.total)
  cat("The returned object is a list containing a matrix of phenotypic data with\n")
  cat("the specified heritability, a vector of unscaled true genetic values,\n")
  cat("to the same total variance as the phenotypic values, and a vector of\n")
  cat("selected individuals, four per family, with the highest phenotype value.\n")
  cat("Both phenotypes and true genetic values include values for parents first followed by progeny.\n")
  return(extraction.info)
}
