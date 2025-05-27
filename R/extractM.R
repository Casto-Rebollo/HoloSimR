#' Extracting Microbiota Composition
#'
#' @description
#' This function extracts the microbiota composition of a population and generates a \emph{txt} file with the data.
#'
#' @param pop An object of class \code{\link{MapPop-class}}.
#' @param sym A value of 1 to consider the symbiosis effect on the microbiota.
#' @param globalSP An object of class \code{\link{GlobalSP}}. If \code{NULL}, it will be retrieved from the global environment.
#' @param name The name of the \emph{txt} file to be generated.
#' @param path The path where the \emph{txt} file will be saved. \emph{Default} is the current directory ("./").
#'
#' @return A matrix representing the microbiota composition.
#' @export
#' @usage extractM(pop,
#'                 sym,
#'                 globalSP = NULL,
#'                 name = NULL,
#'                 path = "./")
#'
#' @examples
#' \dontrun{
#' microbiota <- extractM(Pop, sym = 0, name = "M_gen3_sym")
#' }
#'
#'
#'
extractM <- function(pop, sym, globalSP = NULL, name = NULL, path = "./"){

  if(is.null(globalSP)){
    globalSP = get("gSP", envir = .GlobalEnv)
  }

  if(sym == 1){
    microbiome <- matrix(0,nrow = nInd(pop), ncol = globalSP$nSpecies)
    for (ind in 1:nInd(pop)) {
      microbiome[ind,] <- unlist(pop@misc[[ind]]$M_sym)
    }
  }else{
    ##Saving microbiome
    microbiome <- matrix(0,nrow = nInd(pop), ncol = globalSP$nSpecies)
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
