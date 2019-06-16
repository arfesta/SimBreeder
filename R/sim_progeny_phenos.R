#' Simulate progeny phenos
#'
#' This function simulates progeny phenotypes
#' @param progeny.TGV Object returned from calc_progeny_TGV()
#' @param h2 Heritability used to simulate phenotypes
#' @param cross.design Object returned from create_cross_design()
#' @param E.var Environmental variance that can be used instead of heritability
#' @param folder Folder to output results to. Only matters if save=T
#' @param save logical.  Set to TRUE to save output.  Default: FALSE
#' @keywords sim progeny phenos
#' @export
#' @examples
#' progeny1.PHENOS <- sim_progeny_phenos(progeny.TGV = progeny1.TGV,h2 = .3,cross.design = cross.file)

####Create Unscaled Phenos####
sim_progeny_phenos <- function(progeny.TGV, h2, cross.design, E.var = NULL, save = F,
                               folder = NULL) {
  genetic.values <-progeny.TGV$genetic.values
  total.indiv <-length(genetic.values)
  if(length(E.var) > 0){
    phenos <- genetic.values + rnorm(total.indiv,mean = 0,sd = sqrt(E.var))
  } else {
    phenos <- genetic.values + rnorm(total.indiv,mean = 0,sd = sqrt(var(genetic.values)/h2))}
  true.h2 <- round(var(genetic.values)/var(phenos),2)

  num.progeny <- as.numeric(cross.design$cross.design[1,3])
  crosses <- nrow(cross.design$cross.design)
  fam.mean.g <- vector(); fam.mean.p <- vector()
  first <- 1
  last <- num.progeny
  for(i in 1:crosses){
    fam.mean.g <- c(fam.mean.g,mean(genetic.values[first:last]))
    fam.mean.p <- c(fam.mean.p,mean(phenos[first:last]))
    first <- last +1
    last <- last + num.progeny}
  fam.mean.h2 <-  var(fam.mean.g)/var(fam.mean.p)

  pheno.info<-list(phenos=phenos, genetic.values=genetic.values, true.h2 = true.h2,fam.mean.h2 = fam.mean.h2)
  return(pheno.info)
}
