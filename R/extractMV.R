#' Extracting Microbiota value
#'
#' @description
#' This function extracts the microbiota values of a population
#'
#' @param pop An object of class \code{\link{MapPop-class}}.
#' @param sym A value of 1 to consider the symbiosis effect on the microbiota.
#' @param globalSP An object of class \code{\link{GlobalSP}}. If \code{NULL}, it will be retrieved from the global environment.
#'
#' @return A matrix representing the microbiota composition.
#' @export
#' @usage extractM(pop,
#'                 sym,
#'                 globalSP = NULL)
#'
#' @examples
#' \dontrun{
#' microbiota <- extractM(Pop, sym = 0)
#' }
#'
#'
#'
extractM <- function(pop, sym, globalSP = NULL){

  if(is.null(globalSP)){
    globalSP = get("gSP", envir = .GlobalEnv)
  }

  if(sym == 1){
    microbiome <- matrix(0,nrow = nInd(pop), ncol = length(globalSP$nSpecies))
    for (ind in 1:nInd(pop)) {
      microbiome[ind,] <- unlist(pop@misc[[ind]]$mv_sym)
    }
  }else{
    ##Saving microbiome
    microbiome <- matrix(0,nrow = nInd(pop), ncol = length(globalSP$nSpecies))
    for (ind in 1:nInd(pop)) {
      microbiome[ind,] <- unlist(pop@misc[[ind]]$M)
    }

  }

  if(is.null(name)){
    name <- paste("microbiota_",sym,sep="")
  }
  write.table(microbiome,paste(path,name,sep=""),
              col.names = TRUE, row.names = FALSE, quote = FALSE)

  return(microbiome)
}
