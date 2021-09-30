#' Simulate progeny phenos
#'
#' This function simulates progeny phenotypes
#' @param TGV.object Object returned from calc_TGV()
#' @param h2 Heritability used to simulate phenotypes
#' @param E.var Environmental variance that can be used instead of heritability
#' @keywords sim progeny phenos
#' @export
#' @examples
#' progeny1.PHENOS <- sim_phenos(TGV.object = progeny1.TGV,h2 = .3,cross.design = cross.file)

####Create Unscaled Phenos####
sim_phenos <- function(TGV.object, h2, E.var = NULL) {
  genetic.values <-TGV.object$genetic.values
  total.indiv <-length(genetic.values)
  if(length(E.var) > 0){
    phenos <- genetic.values + rnorm(total.indiv,mean = 0,sd = sqrt(E.var))
  } else {
    phenos <- genetic.values + rnorm(total.indiv,mean = 0,sd = sqrt((var(genetic.values)/h2) - var(genetic.values)))}
  
  pheno.info<-list(phenos=phenos, genetic.values=genetic.values)
  return(pheno.info)
}
