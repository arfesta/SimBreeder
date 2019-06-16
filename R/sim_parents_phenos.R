#' Simulate parental phenotypes
#'
#' This function generates parental phenotypes and is used after creating the founder population
#' @param parents.TGV Object returned from calc_parent_TGV()
#' @param h2 Heritability used to simulate phenotypes
#' @param E.var Environmental variance that can be used instead of heritability
#' @param folder Folder to output results to. Only matters if save=T
#' @param save logical.  Set to TRUE to save output.  Default: FALSE
#' @keywords sim progeny phenos
#' @export
#' @examples
#' parent.PHENOS <- sim_parents_phenos(parents.TGV = parent.TGV,h2 = .3)

####Create Unscaled Parent Phenos####
sim_parents_phenos <- function(parents.TGV, h2, E.var = NULL, save = F, folder = NULL) {
  genetic.values <-parents.TGV$genetic.values
  total.indiv <-length(genetic.values)
  if(length(E.var) > 0){
    phenos <- genetic.values + rnorm(total.indiv,mean = 0,sd = sqrt(E.var))
  } else {
    phenos <- genetic.values + rnorm(total.indiv,mean = 0,sd = sqrt(var(genetic.values)/h2))}
  true.h2 <- round(var(genetic.values)/var(phenos),2)
  pheno.info<-list(phenos=phenos, genetic.values=genetic.values, true.h2=true.h2 )
  return(pheno.info)
}
