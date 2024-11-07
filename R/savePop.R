#' Save Population
#'
#' @description
#' This function allows populations to be stored in a common list.
#'
#' @param pop An object of class \code{\link{MapPop-class}}.
#' @param model The model used for the simulation. Valid values are \dQuote{G}, \dQuote{M}, \dQuote{NMH}, \dQuote{LMH}, \dQuote{MMH}, or \dQuote{HMH}.
#' @param g The current generation of selection. A value is required.
#' @param save A list of previously saved populations. If a variable called \emph{sortedPop} exists, it will be retrieved from the environment and used to continue storing the populations.
#'
#' @return A list of \code{\link{MapPop-class}} objects, ordered by generation.
#' @export
#' @usage savePop(pop, model, g, save = NULL)
#'
#' @examples
#' \dontrun{
#' # Saving the population
#' storedPop <- savePop(Pop, model = "NMH")
#' }
#'
savePop <- function(pop, model,g,save = NULL){

  if(is.null(pop)){
    stop("Population not provided")
  }

  if(is.null(model)){
    stop("Model not provided")
  }

  save_tmp <- setNames(list(pop), paste(model, g, sep = "_"))
  if(exists("sortedPop", envir = parent.frame())){
    save <- get("sortedPop", envir = parent.frame())
    save <- c(save, save_tmp)
  }else{
    if(exists("sortedPop", envir = .GlobalEnv)){
      save <- get("sortedPop", envir = .GlobalEnv)
      save <- save <- c(save, save_tmp)
    }else{
      save <- save_tmp
    }
  }

  return(save)

}
