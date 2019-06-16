#' Extract selections from current generation
#'
#' This function extracts selections from the current generation
#' @param num.cores Number of cores used to create marker based relatedness matrix
#' @param map.info Object returned from create_genetic_map()
#' @param cross.design Object returned from create_cross_design()
#' @param past.tgv Object returned from previous generation
#' @param past.phenos Object returned from previous generation
#' @param parent.info Object returned from previous generation;
#' @param progeny.info Object returned from make_crosses()
#' @param progeny.TGV Object returned from calc_progeny_TGV()
#' @param progeny.phenos Object returned from sim_progeny_phenos()
#' @param selection.strategy Type of selection strategy to use. Currently either: "among.family" or "within.family". Default: "within.familiy"
#' @param among.family.selection Type of selection strategy to use for among family selection: "ABLUP", "GBLUP", "Phenotype"
#' @param within.family.selection Type of selection strategy to use for within family selection: "ABLUP", "GBLUP", "Phenotype"
#' @param num.selections.within.family The number of selections that should be made within each family
#' @param num.selections.among.family The number of selections that should be made among all families
#' @param relationship.matrix.type The type of relationship matrix to estimate for selections: "markers" or "pedigree"
#' @param prefix A prefix to add to the file if writing the file out
#' @param rep.num The rep number of the current generation
#' @param reduced  logical. Set to TRUE if selecting less individuals than are in the total population
#' @param weighted logical. Should relationship matrix be weighted by A matrix prior to BLUP estimation.
#' @param save logical. Saves the output of genetic map (Default: FALSE)
#' @keywords extract selections
#' @export
#' @examples
#' progeny1.extractions <- extract_selections(num.cores = 3,map.info = the.map,cross.design = cross.file,progeny.info = progeny1,
#' parent.info = op.families,progeny.TGV = progeny1.TGV,progeny.phenos = progeny1.PHENOS,num.selections.among.family = 64,reduced = T,
#' relationship.matrix.type = "markers",among.family.selection = "GBLUP",within.family.selection = "NO",selection.strategy = "within.family",
#' past.tgv = parent.TGV,past.phenos = parent.PHENOS)


####Extract the.selections####
extract_selections <- function(num.cores=NULL, map.info, cross.design, past.tgv, past.phenos, parent.info,
                               progeny.info, progeny.TGV, progeny.phenos, selection.strategy,
                               among.family.selection = NULL, within.family.selection = NULL, num.selections.within.family = 1,
                               num.selections.among.family = NULL, relationship.matrix.type = "pedigree",
                               prefix = NULL, rep.num = NULL, reduced = F, weighted=F) {

  library(MatrixModels); library(parallel); library(pedigreemm); library(pedigree)

  # Extracting cross.design, number of parents, complete(full) & progeny pedigree from progeny.info object
  current.cross.design <- cross.design$cross.design   #Get cross design file for current generation
  prog.percross <- as.numeric(current.cross.design[1,3]) #Determine the number of progeny per cross
  prog.pedigree <- cross.design$progeny.pedigree  #Extract the full progeny pedigree for this generation
  generation <- as.numeric(prog.pedigree[1,4])  #The current generation number is the last column in the prog pedigree
  numparents <- cross.design$numparents  #Get number of parents from parent generation
  full.ped <- cross.design$full.pedigree
  selection.ped <- cross.design$selection.ped

  # Extracting parent markers, progeny tgv & markers, as well as creating marker map from progeny.TGV object
  if (generation==1){
    if(reduced){
      parent.markers <- parent.info$genos.3d[map.info$available.Markers,(names(parent.info$pars)),]
    } else { parent.markers <- past.tgv$markers.matrix}} else {parent.markers <- parent.info$all.markers}
  prog.markers <- progeny.TGV$markers.matrix
  progeny.phenos <- progeny.phenos$phenos
  prog.genetic.values <- progeny.TGV$genetic.values
  map.markers <- progeny.TGV$marker.map
  colnames(map.markers) <- c("chr","pos")
  map.markers$chr <- as.character(map.markers$chr)

  # Generate pedigree using synbreed if 1st generation otherwise just use full.ped object
  # Create 2 objects which hold all markers and all phenotypes respectively

  if (generation==1) {
    ped <- pedigree::pedigree(label = full.ped[,1],sire = full.ped[,2],dam = full.ped[,3])
    #all.markers <- rbind(parent.markers,)
    if(reduced){  all.phenos <- c(parent.info$phenos,progeny.phenos)
    all.genetic.vals <- c(parent.info$genetic.values,prog.genetic.values)
    } else { all.phenos <- c(past.phenos$phenos,progeny.phenos)
    all.genetic.vals <- c(past.phenos$genetic.values,prog.genetic.values)}
  } else {
    ped <- pedigree::pedigree(label = selection.ped[,1],sire = selection.ped[,2],dam = selection.ped[,3])
    #all.markers <- rbind(parent.markers,prog.markers)
    all.phenos <- c(parent.info$all.phenos,progeny.phenos)
    all.genetic.vals <- c(past.phenos$genetic.values,prog.genetic.values)
  }

  if (among.family.selection == "ABLUP"){
    the.data <- as.big.matrix(as.matrix(pedigreemm::getAInv(ped)))
    the.data = the.data[i=match(names(progeny.phenos),colnames(the.data)),j=match(names(progeny.phenos),colnames(the.data))]
    s <-  ncol(the.data)/8
    e <- s-1
    first.seq <- seq(1,ncol(the.data),ncol(the.data)/8)
    {
      my.env <- new.env()
      my.env$l1 <-  the.data[i=c(first.seq[1]:(first.seq[1]+e),first.seq[2]:(first.seq[2]+e)),
                             j=c(first.seq[1]:(first.seq[1]+e),first.seq[2]:(first.seq[2]+e)),drop=T]

      my.env$l2 <- the.data[i=c(first.seq[1]:(first.seq[1]+e),first.seq[3]:(first.seq[3]+e)),
                            j=c(first.seq[1]:(first.seq[1]+e),first.seq[3]:(first.seq[3]+e)),drop=T]

      my.env$l3 <- the.data[i=c(first.seq[1]:(first.seq[1]+e),first.seq[4]:(first.seq[4]+e)),
                            j=c(first.seq[1]:(first.seq[1]+e),first.seq[4]:(first.seq[4]+e)),drop=T]

      my.env$l4 <- the.data[i=c(first.seq[1]:(first.seq[1]+e),first.seq[5]:(first.seq[5]+e)),
                            j=c(first.seq[1]:(first.seq[1]+e),first.seq[5]:(first.seq[5]+e)),drop=T]

      my.env$l5 <- the.data[i=c(first.seq[1]:(first.seq[1]+e),first.seq[6]:(first.seq[6]+e)),
                            j=c(first.seq[1]:(first.seq[1]+e),first.seq[6]:(first.seq[6]+e)),drop=T]

      my.env$l6 <- the.data[i=c(first.seq[1]:(first.seq[1]+e),first.seq[7]:(first.seq[7]+e)),
                            j=c(first.seq[1]:(first.seq[1]+e),first.seq[7]:(first.seq[7]+e)),drop=T]

      my.env$l7 <- the.data[i=c(first.seq[1]:(first.seq[1]+e),first.seq[8]:(first.seq[8]+e)),
                            j=c(first.seq[1]:(first.seq[1]+e),first.seq[8]:(first.seq[8]+e)),drop=T]

      my.env$l8 <- the.data[i=c(first.seq[2]:(first.seq[2]+e),first.seq[3]:(first.seq[3]+e)),
                            j=c(first.seq[2]:(first.seq[2]+e),first.seq[3]:(first.seq[3]+e)),drop=T]

      my.env$l9 <- the.data[i=c(first.seq[2]:(first.seq[2]+e),first.seq[4]:(first.seq[4]+e)),
                            j=c(first.seq[2]:(first.seq[2]+e),first.seq[4]:(first.seq[4]+e)),drop=T]

      my.env$l10 <- the.data[i=c(first.seq[2]:(first.seq[2]+e),first.seq[5]:(first.seq[5]+e)),
                             j=c(first.seq[2]:(first.seq[2]+e),first.seq[5]:(first.seq[5]+e)),drop=T]

      my.env$l11 <- the.data[i=c(first.seq[2]:(first.seq[2]+e),first.seq[6]:(first.seq[6]+e)),
                             j=c(first.seq[2]:(first.seq[2]+e),first.seq[6]:(first.seq[6]+e)),drop=T]

      my.env$l12 <- the.data[i=c(first.seq[2]:(first.seq[2]+e),first.seq[7]:(first.seq[7]+e)),
                             j=c(first.seq[2]:(first.seq[2]+e),first.seq[7]:(first.seq[7]+e)),drop=T]

      my.env$l13 <- the.data[i=c(first.seq[2]:(first.seq[2]+e),first.seq[8]:(first.seq[8]+e)),
                             j=c(first.seq[2]:(first.seq[2]+e),first.seq[8]:(first.seq[8]+e)),drop=T]

      my.env$l14 <- the.data[i=c(first.seq[3]:(first.seq[3]+e),first.seq[4]:(first.seq[4]+e)),
                             j=c(first.seq[3]:(first.seq[3]+e),first.seq[4]:(first.seq[4]+e)),drop=T]

      my.env$l15 <- the.data[i=c(first.seq[3]:(first.seq[3]+e),first.seq[5]:(first.seq[5]+e)),
                             j=c(first.seq[3]:(first.seq[3]+e),first.seq[5]:(first.seq[5]+e)),drop=T]

      my.env$l16 <- the.data[i=c(first.seq[3]:(first.seq[3]+e),first.seq[6]:(first.seq[6]+e)),
                             j=c(first.seq[3]:(first.seq[3]+e),first.seq[6]:(first.seq[6]+e)),drop=T]

      my.env$l17 <- the.data[i=c(first.seq[3]:(first.seq[3]+e),first.seq[7]:(first.seq[7]+e)),
                             j=c(first.seq[3]:(first.seq[3]+e),first.seq[7]:(first.seq[7]+e)),drop=T]

      my.env$l18 <- the.data[i=c(first.seq[3]:(first.seq[3]+e),first.seq[8]:(first.seq[8]+e)),
                             j=c(first.seq[3]:(first.seq[3]+e),first.seq[8]:(first.seq[8]+e)),drop=T]

      my.env$l19 <- the.data[i=c(first.seq[4]:(first.seq[4]+e),first.seq[5]:(first.seq[5]+e)),
                             j=c(first.seq[4]:(first.seq[4]+e),first.seq[5]:(first.seq[5]+e)),drop=T]

      my.env$l20 <- the.data[i=c(first.seq[4]:(first.seq[4]+e),first.seq[6]:(first.seq[6]+e)),
                             j=c(first.seq[4]:(first.seq[4]+e),first.seq[6]:(first.seq[6]+e)),drop=T]

      my.env$l21<- the.data[i=c(first.seq[4]:(first.seq[4]+e),first.seq[7]:(first.seq[7]+e)),
                            j=c(first.seq[4]:(first.seq[4]+e),first.seq[7]:(first.seq[7]+e)),drop=T]

      my.env$l22<- the.data[i=c(first.seq[4]:(first.seq[4]+e),first.seq[8]:(first.seq[8]+e)),
                            j=c(first.seq[4]:(first.seq[4]+e),first.seq[8]:(first.seq[8]+e)),drop=T]

      my.env$l23 <- the.data[i=c(first.seq[5]:(first.seq[5]+e),first.seq[6]:(first.seq[6]+e)),
                             j=c(first.seq[5]:(first.seq[5]+e),first.seq[6]:(first.seq[6]+e)),drop=T]

      my.env$l24 <- the.data[i=c(first.seq[5]:(first.seq[5]+e),first.seq[7]:(first.seq[7]+e)),

                             j=c(first.seq[5]:(first.seq[5]+e),first.seq[7]:(first.seq[7]+e)),drop=T]

      my.env$l25 <- the.data[i=c(first.seq[5]:(first.seq[5]+e),first.seq[8]:(first.seq[8]+e)),
                             j=c(first.seq[5]:(first.seq[5]+e),first.seq[8]:(first.seq[8]+e)),drop=T]

      my.env$l26 <- the.data[i=c(first.seq[6]:(first.seq[6]+e),first.seq[7]:(first.seq[7]+e)),
                             j=c(first.seq[6]:(first.seq[6]+e),first.seq[7]:(first.seq[7]+e)),drop=T]

      my.env$l27 <- the.data[i=c(first.seq[6]:(first.seq[6]+e),first.seq[8]:(first.seq[8]+e)),
                             j=c(first.seq[6]:(first.seq[6]+e),first.seq[8]:(first.seq[8]+e)),drop=T]

      my.env$l28 <- the.data[i=c(first.seq[7]:(first.seq[7]+e),first.seq[8]:(first.seq[8]+e)),
                             j=c(first.seq[7]:(first.seq[7]+e),first.seq[8]:(first.seq[8]+e)),drop=T]
    }

    gc()

    {
      s.phenos1 <- progeny.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[2]:(first.seq[2]+e))]
      s.phenos2 <- progeny.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[3]:(first.seq[3]+e))]
      s.phenos3 <- progeny.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[4]:(first.seq[4]+e))]
      s.phenos4 <- progeny.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[5]:(first.seq[5]+e))]
      s.phenos5 <- progeny.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos6 <- progeny.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos7 <- progeny.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[8]:(first.seq[8]+e))]

      s.phenos8 <- progeny.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[3]:(first.seq[3]+e))]
      s.phenos9 <- progeny.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[4]:(first.seq[4]+e))]
      s.phenos10 <- progeny.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[5]:(first.seq[5]+e))]
      s.phenos11 <- progeny.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos12 <- progeny.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos13 <- progeny.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[8]:(first.seq[8]+e))]

      s.phenos14 <- progeny.phenos[c(first.seq[3]:(first.seq[3]+e),first.seq[4]:(first.seq[4]+e))]
      s.phenos15 <- progeny.phenos[c(first.seq[3]:(first.seq[3]+e),first.seq[5]:(first.seq[5]+e))]
      s.phenos16 <- progeny.phenos[c(first.seq[3]:(first.seq[3]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos17 <- progeny.phenos[c(first.seq[3]:(first.seq[3]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos18 <- progeny.phenos[c(first.seq[3]:(first.seq[3]+e),first.seq[8]:(first.seq[8]+e))]

      s.phenos19 <- progeny.phenos[c(first.seq[4]:(first.seq[4]+e),first.seq[5]:(first.seq[5]+e))]
      s.phenos20 <- progeny.phenos[c(first.seq[4]:(first.seq[4]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos21 <- progeny.phenos[c(first.seq[4]:(first.seq[4]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos22 <- progeny.phenos[c(first.seq[4]:(first.seq[4]+e),first.seq[8]:(first.seq[8]+e))]

      s.phenos23 <- progeny.phenos[c(first.seq[5]:(first.seq[5]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos24 <- progeny.phenos[c(first.seq[5]:(first.seq[5]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos25 <- progeny.phenos[c(first.seq[5]:(first.seq[5]+e),first.seq[8]:(first.seq[8]+e))]

      s.phenos26 <- progeny.phenos[c(first.seq[6]:(first.seq[6]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos27 <- progeny.phenos[c(first.seq[6]:(first.seq[6]+e),first.seq[8]:(first.seq[8]+e))]

      s.phenos28 <- progeny.phenos[c(first.seq[7]:(first.seq[7]+e),first.seq[8]:(first.seq[8]+e))]

      pheno.list <- list(s.phenos1,s.phenos2,s.phenos3,s.phenos4,s.phenos5,s.phenos6,s.phenos7,s.phenos8,s.phenos9,
                         s.phenos10,s.phenos11,s.phenos12,s.phenos13,s.phenos14,s.phenos15,s.phenos16,s.phenos17,
                         s.phenos18,s.phenos19,s.phenos20,s.phenos21,s.phenos22,s.phenos23,s.phenos24,s.phenos25,s.phenos26,s.phenos27,s.phenos28)
    }


    apply_my_fun <- function(idx, env) {
      eval(parse(text = paste0("result <- env$", ls(env)[idx])))
      the.data <- as.matrix(result)
      #the.data <- solve(the.data)
      n.col <- ncol(the.data)
      h.2 <- var(prog.genetic.values)/var(progeny.phenos)
      lambda <- (1-h.2)/h.2
      I <- diag(n.col)
      s.phenos <- pheno.list[[idx]]
      sol <- solve(rBind(cBind(n.col, t(rep(1,n.col))),
                         cBind(rep(1,n.col), (I + (lambda*the.data)))),
                   matrix(c(sum(s.phenos),
                            c(as.vector(s.phenos)))))
      sol <- sol[-1,1]
      sol
    }

    index <- 1:28
    names(index) <- 1:28

    the.blups <- mclapply(index,apply_my_fun,env = my.env, mc.cores=28)
    rm(my.env);gc()

    {
      SG1 <- apply(rbind(the.blups[[1]][1:(length(s.phenos1)/2)],
                         the.blups[[2]][1:(length(s.phenos1)/2)],
                         the.blups[[3]][1:(length(s.phenos1)/2)],
                         the.blups[[4]][1:(length(s.phenos1)/2)],
                         the.blups[[5]][1:(length(s.phenos1)/2)],
                         the.blups[[6]][1:(length(s.phenos1)/2)],
                         the.blups[[7]][1:(length(s.phenos1)/2)]),2,mean)

      SG2 <- apply(rbind(the.blups[[1]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[8]][1:(length(s.phenos1)/2)],
                         the.blups[[9]][1:(length(s.phenos1)/2)],
                         the.blups[[10]][1:(length(s.phenos1)/2)],
                         the.blups[[11]][1:(length(s.phenos1)/2)],
                         the.blups[[12]][1:(length(s.phenos1)/2)],
                         the.blups[[13]][1:(length(s.phenos1)/2)]),2,mean)

      SG3 <- apply(rbind(the.blups[[2]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[8]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[14]][1:(length(s.phenos1)/2)],
                         the.blups[[15]][1:(length(s.phenos1)/2)],
                         the.blups[[16]][1:(length(s.phenos1)/2)],
                         the.blups[[17]][1:(length(s.phenos1)/2)],
                         the.blups[[18]][1:(length(s.phenos1)/2)]),2,mean)

      SG4 <- apply(rbind(the.blups[[3]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[9]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[14]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[19]][1:(length(s.phenos1)/2)],
                         the.blups[[20]][1:(length(s.phenos1)/2)],
                         the.blups[[21]][1:(length(s.phenos1)/2)],
                         the.blups[[22]][1:(length(s.phenos1)/2)]),2,mean)

      SG5 <- apply(rbind(the.blups[[4]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[10]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[15]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[19]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[23]][1:(length(s.phenos1)/2)],
                         the.blups[[24]][1:(length(s.phenos1)/2)],
                         the.blups[[25]][1:(length(s.phenos1)/2)]),2,mean)

      SG6 <- apply(rbind(the.blups[[5]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[11]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[16]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[20]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[23]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[26]][1:(length(s.phenos1)/2)],
                         the.blups[[27]][1:(length(s.phenos1)/2)]),2,mean)

      SG7 <- apply(rbind(the.blups[[6]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[12]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[17]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[21]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[24]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[26]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[28]][1:(length(s.phenos1)/2)]),2,mean)

      SG8 <- apply(rbind(the.blups[[7]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[13]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[18]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[22]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[25]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[27]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[28]][((length(s.phenos1)/2)+1):length(s.phenos1)]),2,mean)
    }
    progeny.blups <- c(SG1,SG2,SG3,SG4,SG5,SG6,SG7,SG8); names(progeny.blups) <- names(progeny.phenos)
    first <- 1
    last <- prog.percross
    mean.progeny.blups <- vector()
    for(each in 1:nrow(current.cross.design)){
      mean.progeny.blups <- c(mean.progeny.blups,mean(progeny.blups[first:last]))
      first <- last+1
      last <- first+prog.percross-1
    }
    sorted.mean.progeny.blups <- sort(mean.progeny.blups,decreasing=T)
    top.280.families <- match(sorted.mean.progeny.blups, mean.progeny.blups)

    the.selections <- vector()
    first.in.family <- 1
    for(family in 1:length(cross.design[,3])){
      num.offspring <- as.numeric(cross.design[family,3])
      last.in.family <- num.offspring + first.in.family - 1
      temp <- (progeny.blups[first.in.family:last.in.family]) ; sorted <-sort(temp,decreasing=TRUE)
      best.one <- which(temp==sorted[1:num.selections.within.family]) ; selected <- best.one + first.in.family - 1 ;the.selections <- c(the.selections, selected)
      first.in.family <- last.in.family + 1 }
    if(length(num.selections.among.family) > 0) {
      the.selections <- the.selections[top.280.families[1:num.selections.among.family]]} else {
        the.selections <- the.selections[top.280.families]
      }

  }
  if(among.family.selection == "GBLUP"){
    m <- map.info$available.Markers
    #if(generation==1){parentmarkers <- parent.info$genos.3d[m,,]} else {parentmarkers <- parent.info$all.markers}
    prog.markers <- progeny.info$genos.3d[m,,]
    prog.1 <- prog.markers[,,1]; colnames(prog.1) <- names(progeny.phenos)
    prog.2 <- prog.markers[,,2]; colnames(prog.2) <- names(progeny.phenos)
    #allele1 <- cbind(parentmarkers[,,1],prog.1)
    #allele2 <- cbind(parentmarkers[,,2],prog.2)
    A <- (as.matrix(pedigreemm::getA(ped)))
    A = A[match(names(progeny.phenos),colnames(A)),match(names(progeny.phenos),colnames(A))]
    #this.order <- match(colnames(allele1),rownames(A))
    #A <- A[c(this.order),c(this.order)]
    #rm(prog.markers,parentmarkers)


    prog.list <- lapply(1:ncol(prog.1),function(x){
      cbind(prog.1[,x],prog.2[,x])
    })
    gc()
    pt <- proc.time()
    genetic.relationships <- mclapply(1:length(prog.list),function(each.prog){
      if(each.prog == length(prog.list)){
        1
      } else {

        first.prog <- prog.list[[each.prog]]

        answer = rep(0,dim(prog.markers[,each.prog:ncol(prog.markers),1])[1] * dim(prog.markers[,each.prog:ncol(prog.markers),1])[2])


        a1.matches.b1.b2 = which(first.prog[,1] == prog.markers[,each.prog:ncol(prog.markers),1] & first.prog[,1] == prog.markers[,each.prog:ncol(prog.markers),2])
        a1.matches.b1.but.not.b2 = which(first.prog[,1] == prog.markers[,each.prog:ncol(prog.markers),1] & first.prog[,1] != prog.markers[,each.prog:ncol(prog.markers),2])
        a1.matches.b2.but.not.b1 = which(first.prog[,1] == prog.markers[,each.prog:ncol(prog.markers),2] & first.prog[,1] != prog.markers[,each.prog:ncol(prog.markers),1])

        a2.matches.b1.b2 = which(first.prog[,2] == prog.markers[,each.prog:ncol(prog.markers),1] & first.prog[,2] == prog.markers[,each.prog:ncol(prog.markers),2])
        a2.matches.b1.but.not.b2 = which(first.prog[,2] == prog.markers[,each.prog:ncol(prog.markers),1] & first.prog[,2] != prog.markers[,each.prog:ncol(prog.markers),2])
        a2.matches.b2.but.not.b1 = which(first.prog[,2] == prog.markers[,each.prog:ncol(prog.markers),2] & first.prog[,2] != prog.markers[,each.prog:ncol(prog.markers),1])

        #These individuals are AA/AA compared to AA/AA
        homozygous.identical = a1.matches.b1.b2[which(a1.matches.b1.b2 %in% a2.matches.b1.b2)]
        if(length(homozygous.identical) > 0 ){answer[homozygous.identical] = 2}

        #These individuals are AA/BB compared to AA/BB
        heterozygous.identical = a1.matches.b1.but.not.b2[which(a1.matches.b1.but.not.b2 %in% a2.matches.b2.but.not.b1)]
        if(length(heterozygous.identical) > 0){answer[heterozygous.identical] = 2}

        #These indivdiuals are AA/BB compared to BB/AA
        other.het.identical = a1.matches.b2.but.not.b1[which(a1.matches.b2.but.not.b1 %in% a2.matches.b1.but.not.b2)]
        if(length(other.het.identical) > 0){answer[other.het.identical] = 2}

        #These individuals are AA/AA compared to AA/BB
        both.a.match.first.other = a1.matches.b1.but.not.b2[which(a1.matches.b1.but.not.b2 %in% a2.matches.b1.but.not.b2)]
        if(length(both.a.match.first.other) > 0){answer[both.a.match.first.other] = 1}

        #These indivdiuals are AA/AA compared to BB/AA
        both.a.match.second.other = a1.matches.b2.but.not.b1[which(a1.matches.b2.but.not.b1 %in% a2.matches.b2.but.not.b1)]
        if(length(both.a.match.second.other) > 0){answer[both.a.match.second.other] = 1}



        #These individuals are AA/BB compared to AA/AA
        if(length(which(a1.matches.b1.b2 %in% a2.matches.b1.b2)) > 0) {
          first.a.matches.both = a1.matches.b1.b2[-c(which(a1.matches.b1.b2 %in% a2.matches.b1.b2))]} else {first.a.matches.both = a1.matches.b1.b2}
        if(length(first.a.matches.both) > 0){answer[first.a.matches.both] = 1}

        #These indivdiuals are AA/BB compared to BB/BB
        if(length(which(a2.matches.b1.b2 %in% a1.matches.b1.b2)) > 0) {
          second.a.matches.both = a2.matches.b1.b2[-c(which(a2.matches.b1.b2 %in% a1.matches.b1.b2))]} else {second.a.matches.both = a2.matches.b1.b2}
        if(length(second.a.matches.both) > 0){answer[second.a.matches.both] = 1}



        #These indivdiuals are AA/BB compared to AA/CC
        only.first.allele.match = a1.matches.b1.but.not.b2[-c(which(a1.matches.b1.but.not.b2 %in% a2.matches.b2.but.not.b1))]
        if(length(only.first.allele.match) > 0){answer[only.first.allele.match] = 1}

        #These individuals are AA/BB compared to CC/AA
        only.first.allele.match2 = a1.matches.b2.but.not.b1[-c(which(a1.matches.b2.but.not.b1 %in% a2.matches.b1.but.not.b2))]
        if(length(only.first.allele.match2) > 0){answer[only.first.allele.match2] = 1}

        #These indivdiuals are AA/BB compared to BB/CC
        only.second.allele.match = a2.matches.b1.but.not.b2[-c(which(a2.matches.b1.but.not.b2 %in% a1.matches.b2.but.not.b1))]
        if(length(only.second.allele.match) > 0){answer[only.second.allele.match] = 1}

        #These indivdiuals are AA/BB compared to CC/BB
        only.second.alle.match2 = a2.matches.b2.but.not.b1[-c(which(a2.matches.b2.but.not.b1 %in% a1.matches.b1.but.not.b2))]
        if(length(only.second.alle.match2) > 0){answer[only.second.alle.match2] = 1}

        the.ans = matrix(answer,nrow = 120,ncol=dim(prog.markers[,each.prog:ncol(prog.markers),1])[2])


        ll <- colSums(the.ans)
        ll/240
      }},mc.cores=60)
    #proc.time() -pt
    #232 seconds
    g.mat <- matrix(0,nrow=length(prog.list),ncol=length(prog.list))
    g.mat[lower.tri(g.mat,diag=T)] <- unlist(genetic.relationships)
    t.gmat2 <- t(g.mat)
    g.mat[upper.tri(g.mat,diag=F)] <- t.gmat2[upper.tri(t.gmat2,diag=F)]
    proc.time() -pt
    rm(t.gmat2); gc()
    the.data <- g.mat*.99 + as.matrix(A)*.01
    the.data <- solve(the.data)
    n.col <- ncol(the.data)
    h.2 <- var(prog.genetic.values)/var(progeny.phenos)
    lambda <- (1-h.2)/h.2
    I <- diag(n.col)
    #s.phenos <- pheno.list[[idx]]
    sol <- solve(rBind(cBind(n.col, t(rep(1,n.col))),
                       cBind(rep(1,n.col), (I + (lambda*the.data)))),
                 matrix(c(sum(progeny.phenos),
                          c(as.vector(progeny.phenos)))))
    sol <- sol[-1,1]
    gprogeny.blups1 <- unlist(sol)
    names(gprogeny.blups1) <- names(progeny.phenos)

    first <- 1
    last <- prog.percross
    mean.progeny.blups <- vector()
    for(each in 1:nrow(current.cross.design)){
      mean.progeny.blups <- c(mean.progeny.blups,mean(gprogeny.blups1[first:last]))
      first <- last+1
      last <- first+prog.percross-1
    }
    sorted.mean.progeny.blups <- sort(mean.progeny.blups,decreasing=T)
    top.280.families <- match(sorted.mean.progeny.blups, mean.progeny.blups)


    the.selections <- vector()
    first.in.family <- 1
    for(family in 1:length(cross.design[,3])){
      num.offspring <- as.numeric(cross.design[family,3])
      last.in.family <- num.offspring + first.in.family - 1
      temp <- (gprogeny.blups1[first.in.family:last.in.family]) ; sorted <-sort(temp,decreasing=TRUE)
      best.one <- which(temp==sorted[1:num.selections.within.family]) ; selected <- best.one + first.in.family - 1 ;the.selections <- c(the.selections, selected)
      first.in.family <- last.in.family + 1 }

    if(length(num.selections.among.family) > 0) {
      the.selections <- the.selections[top.280.families[1:num.selections.among.family]]} else {
        the.selections <- the.selections[top.280.families]
      }
  }
  if (among.family.selection == "Phenotype"){
    progeny.blups <- progeny.phenos; names(progeny.blups) <- names(progeny.phenos)
    first <- 1
    last <- prog.percross
    mean.progeny.blups <- vector()
    for(each in 1:nrow(current.cross.design)){
      mean.progeny.blups <- c(mean.progeny.blups,mean(progeny.blups[first:last]))
      first <- last+1
      last <- first+prog.percross-1
    }
    sorted.mean.progeny.blups <- sort(mean.progeny.blups,decreasing=T)
    top.280.families <- match(sorted.mean.progeny.blups, mean.progeny.blups)

    the.selections <- vector()
    first.in.family <- 1
    for(family in 1:length(current.cross.design[,3])){
      num.offspring <- as.numeric(current.cross.design[family,3])
      last.in.family <- num.offspring + first.in.family - 1
      temp <- (progeny.blups[first.in.family:last.in.family]) ; sorted <-sort(temp,decreasing=TRUE)
      best.one <- which(temp==sorted[1:num.selections.within.family]) ; selected <- best.one + first.in.family - 1 ;the.selections <- c(the.selections, selected)
      first.in.family <- last.in.family + 1 }
    if(length(num.selections.among.family) > 0) {
      the.selections <- the.selections[top.280.families[1:num.selections.among.family]]} else {
        the.selections <- the.selections[top.280.families]
      }
  }
  # Extract phenotypes of new selections
  selection.phenos<- progeny.phenos[the.selections]
  new.pars.genval <- prog.genetic.values[the.selections]
  if (among.family.selection == "GBLUP"){selection.EBVs <- gprogeny.blups1[the.selections]} else {selection.EBVs <- progeny.blups[the.selections]}
  sorted.top192 <- sort(selection.phenos,decreasing=T)
  #all.markers <- rbind(parent.markers,prog.markers[the.selections,])

  if(generation==1){
    if(reduced){  all.phenos <- c(parent.info$phenos,selection.phenos)
    all.genetic.vals <- c(parent.info$genetic.values,new.pars.genval)
    } else {all.phenos <- c(past.phenos$phenos,selection.phenos)
    all.genetic.vals <- c(past.phenos$genetic.values,new.pars.genval)}
  } else { all.phenos <- c(parent.info$all.phenos,selection.phenos)
  all.genetic.vals <- c(parent.info$all.genetic.vals,new.pars.genval)}
  ped <- full.ped[match(names(all.phenos),full.ped[,1]),]

  # Extract genotypes and parent ids of new selections
  colnames(progeny.info$genos.3d) <- names(progeny.phenos)
  new.parent.genos <- progeny.info$genos.3d[,names(the.selections),]
  numselections <- dim(new.parent.genos)[2]
  select.ids <- as.numeric(names(new.pars.genval))-numparents
  select.ped.ids <- as.numeric(names(new.pars.genval))
  if (selection.strategy == "self.test") {colnames(new.parent.genos) <- select.ped.ids} else {colnames(new.parent.genos) <- select.ped.ids}
  if(generation==1){
    all.markers1 <- cbind(parent.info$genos.3d[map.info$available.Markers,,1],new.parent.genos[map.info$available.Markers,,1])
    all.markers2 <- cbind(parent.info$genos.3d[map.info$available.Markers,,2],new.parent.genos[map.info$available.Markers,,2])
    all.markers <- abind(all.markers1,all.markers2,along=3)} else {
      all.markers1 <- cbind(parent.info$all.markers[,,1],new.parent.genos[map.info$available.Markers,,1])
      all.markers2 <- cbind(parent.info$all.markers[,,2],new.parent.genos[map.info$available.Markers,,2])
      all.markers <- abind(all.markers1,all.markers2,along=3)
    }
  # Calculate:
  #       Inbreding level of progeny/selections
  #       Genetic Variance of progeny and selections (Bulmer Effect)
  f.ped <- pedigree::pedigree(label = full.ped[,1],sire = full.ped[,2],dam = full.ped[,3])
  pedigree.inbreeding <- pedigreemm::inbreeding(f.ped)
  names(pedigree.inbreeding) <- full.ped[,1]
  progeny.inbreeding <- pedigree.inbreeding[prog.pedigree[,1]]
  selections.inbreeding <- progeny.inbreeding[names(the.selections)]
  progeny.gen.var <- var(prog.genetic.values[names(the.selections)])
  bulmer.effect <- var(prog.genetic.values[names(the.selections)]) - var(prog.genetic.values)

  #Calculate the number of deletrious alleles for each selection
  result <- sapply(rep(1:length(the.selections),1),function(x) {out <-length(which(new.parent.genos[map.info$QTLSNP.loci,,1][,x]=="c")); out2 <- length(which(new.parent.genos[map.info$QTLSNP.loci,,2][,x]=="c"))
  outer <- out+out2})
  result <- unlist(result)

  if(relationship.matrix.type == "pedigree"){
    if(generation == 1){A <- as.matrix(pedigreemm::getA(f.ped))
    } else {
      selection.ped <- cross.design$selection.ped
      s.ped <- pedigree::pedigree(label = selection.ped[,1],sire = selection.ped[,2],dam = selection.ped[,3])
      A <- as.matrix(pedigreemm::getA(s.ped))}
    l <- match(names(selection.phenos), rownames(A))
    rel.mat <- A[l,l]
  }
  if(relationship.matrix.type == "markers"){
    #m <- map.info$available.Markers
    #progeny1.markers <- new.parent.genos[m,,]
    allele1 <- all.markers[,,1]
    allele2 <- all.markers[,,2]
    allele1 <- allele1[,which(colnames(allele1) %in% names(selection.phenos))]
    allele2 <- allele2[,which(colnames(allele2) %in% names(selection.phenos))]
    all.alleles <- abind(allele1,allele2,along=3)
    prog.list <- lapply(1:ncol(allele1),function(x){
      cbind(allele1[,x],allele2[,x])
    })

    genetic.relationships <- mclapply(1:length(prog.list),function(each.prog){
      if(each.prog == length(prog.list)){
        1
      } else {

        first.prog <- prog.list[[each.prog]]

        answer = rep(0,dim(all.alleles[,each.prog:ncol(all.alleles),1])[1] * dim(all.alleles[,each.prog:ncol(all.alleles),1])[2])


        a1.matches.b1.b2 = which(first.prog[,1] == all.alleles[,each.prog:ncol(all.alleles),1] & first.prog[,1] == all.alleles[,each.prog:ncol(all.alleles),2])
        a1.matches.b1.but.not.b2 = which(first.prog[,1] == all.alleles[,each.prog:ncol(all.alleles),1] & first.prog[,1] != all.alleles[,each.prog:ncol(all.alleles),2])
        a1.matches.b2.but.not.b1 = which(first.prog[,1] == all.alleles[,each.prog:ncol(all.alleles),2] & first.prog[,1] != all.alleles[,each.prog:ncol(all.alleles),1])

        a2.matches.b1.b2 = which(first.prog[,2] == all.alleles[,each.prog:ncol(all.alleles),1] & first.prog[,2] == all.alleles[,each.prog:ncol(all.alleles),2])
        a2.matches.b1.but.not.b2 = which(first.prog[,2] == all.alleles[,each.prog:ncol(all.alleles),1] & first.prog[,2] != all.alleles[,each.prog:ncol(all.alleles),2])
        a2.matches.b2.but.not.b1 = which(first.prog[,2] == all.alleles[,each.prog:ncol(all.alleles),2] & first.prog[,2] != all.alleles[,each.prog:ncol(all.alleles),1])

        #These individuals are AA/AA compared to AA/AA
        homozygous.identical = a1.matches.b1.b2[which(a1.matches.b1.b2 %in% a2.matches.b1.b2)]
        if(length(homozygous.identical) > 0 ){answer[homozygous.identical] = 2}

        #These individuals are AA/BB compared to AA/BB
        heterozygous.identical = a1.matches.b1.but.not.b2[which(a1.matches.b1.but.not.b2 %in% a2.matches.b2.but.not.b1)]
        if(length(heterozygous.identical) > 0){answer[heterozygous.identical] = 2}

        #These indivdiuals are AA/BB compared to BB/AA
        other.het.identical = a1.matches.b2.but.not.b1[which(a1.matches.b2.but.not.b1 %in% a2.matches.b1.but.not.b2)]
        if(length(other.het.identical) > 0){answer[other.het.identical] = 2}

        #These individuals are AA/AA compared to AA/BB
        both.a.match.first.other = a1.matches.b1.but.not.b2[which(a1.matches.b1.but.not.b2 %in% a2.matches.b1.but.not.b2)]
        if(length(both.a.match.first.other) > 0){answer[both.a.match.first.other] = 1}

        #These indivdiuals are AA/AA compared to BB/AA
        both.a.match.second.other = a1.matches.b2.but.not.b1[which(a1.matches.b2.but.not.b1 %in% a2.matches.b2.but.not.b1)]
        if(length(both.a.match.second.other) > 0){answer[both.a.match.second.other] = 1}



        #These individuals are AA/BB compared to AA/AA
        if(length(which(a1.matches.b1.b2 %in% a2.matches.b1.b2)) > 0) {
          first.a.matches.both = a1.matches.b1.b2[-c(which(a1.matches.b1.b2 %in% a2.matches.b1.b2))]} else {first.a.matches.both = a1.matches.b1.b2}
        if(length(first.a.matches.both) > 0){answer[first.a.matches.both] = 1}

        #These indivdiuals are AA/BB compared to BB/BB
        if(length(which(a2.matches.b1.b2 %in% a1.matches.b1.b2)) > 0) {
          second.a.matches.both = a2.matches.b1.b2[-c(which(a2.matches.b1.b2 %in% a1.matches.b1.b2))]} else {second.a.matches.both = a2.matches.b1.b2}
        if(length(second.a.matches.both) > 0){answer[second.a.matches.both] = 1}



        #These indivdiuals are AA/BB compared to AA/CC
        only.first.allele.match = a1.matches.b1.but.not.b2[-c(which(a1.matches.b1.but.not.b2 %in% a2.matches.b2.but.not.b1))]
        if(length(only.first.allele.match) > 0){answer[only.first.allele.match] = 1}

        #These individuals are AA/BB compared to CC/AA
        only.first.allele.match2 = a1.matches.b2.but.not.b1[-c(which(a1.matches.b2.but.not.b1 %in% a2.matches.b1.but.not.b2))]
        if(length(only.first.allele.match2) > 0){answer[only.first.allele.match2] = 1}

        #These indivdiuals are AA/BB compared to BB/CC
        only.second.allele.match = a2.matches.b1.but.not.b2[-c(which(a2.matches.b1.but.not.b2 %in% a1.matches.b2.but.not.b1))]
        if(length(only.second.allele.match) > 0){answer[only.second.allele.match] = 1}

        #These indivdiuals are AA/BB compared to CC/BB
        only.second.alle.match2 = a2.matches.b2.but.not.b1[-c(which(a2.matches.b2.but.not.b1 %in% a1.matches.b1.but.not.b2))]
        if(length(only.second.alle.match2) > 0){answer[only.second.alle.match2] = 1}

        the.ans = matrix(answer,nrow = 120,ncol=dim(all.alleles[,each.prog:ncol(all.alleles),1])[2])


        ll <- colSums(the.ans)
        ll/240
      }},mc.cores=60)

    g.mat <- matrix(0,nrow=length(prog.list),ncol=length(prog.list))
    g.mat[lower.tri(g.mat,diag=T)] <- unlist(genetic.relationships)
    t.gmat2 <- t(g.mat)
    g.mat[upper.tri(g.mat,diag=F)] <- t.gmat2[upper.tri(t.gmat2,diag=F)]
    rm(t.gmat2); gc()
    rownames(g.mat) <- colnames(allele1);
    colnames(g.mat) <- colnames(allele1)
    l <- which(rownames(g.mat) %in% names(selection.phenos))
    rel.mat <- g.mat[l,l]
  }

  extraction.info<-list(relmat=rel.mat, delt.alleles=result, selections=the.selections, bulmer.effect=bulmer.effect,
                        select.EBVs = selection.EBVs, all.genetic.vals=all.genetic.vals,
                        selection.phenos=selection.phenos,ped=ped,
                        prog.inbred.level=progeny.inbreeding, select.inbred.level=selections.inbreeding,
                        genos.3d=new.parent.genos, num.parents=numselections,select.genval=new.pars.genval, fullped=full.ped,
                        par.ids=select.ids,select.ped.ids=select.ped.ids,all.markers=all.markers,all.phenos=all.phenos, cumulative.total=cross.design$cumul.total)
  cat("The returned object is a list containing a matrix of phenotypic data with\n")
  cat("the specified heritability, a vector of unscaled true genetic values,\n")
  cat("to the same total variance as the phenotypic values, and a vector of\n" )
  cat("selected individuals, four per family, with the highest phenotype value.\n" )
  cat("Both phenotypes and true genetic values include values for parents first followed by progeny.\n")
  return(extraction.info)
  # # # # #   End of createPhenos() function # # # # # # #
}
