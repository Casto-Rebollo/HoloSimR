#' Fill Species Information
#'
#' @description
#' This function incorporates microbial species abundances and the microbiota value into an object of class \code{\link{Pop-class}}.
#'
#' @param pop An object of class \code{\link{Pop-class}} generated using the \pkg{\link{AlphaSimR}} package, representing the current population.
#' @param mbiome A microbiota abundance matrix generated using the function \code{\link{makeM}}. If \code{NULL}, \code{mbiome} will be extracted from the global environment.
#' @param w A vector with the \eqn{\omega} effect scaled according to the base population.
#' @param sym A value of 1 to consider the symbiosis effect on the microbiota. \emph{Default} is 0.
#'
#' @return An object of class \code{\link{Pop-class}} with updated microbiota information.
#' @export
#' @usage fillSp(pop, w = NULL, mbiome = NULL, sym = 0)
#'
fillSp <- function(pop, w = NULL,mbiome = NULL, sym = 0) {

  if(exists("basePop")){
    basePop <- get("basePop", envir = .GlobalEnv)
    if(is.null(mbiome)){
      mbiome = get("mbiome", envir = .GlobalEnv)
    }
    if(is.null(w)){
      if(sym == 1){
        w = basePop[["w_scale1"]]
        mean_base = basePop[["mean_base_sym"]]

      }else{
        w = basePop[["w_scale0"]]
        mean_base = basePop[["mean_base"]]

      }
    }

  }else{
    if(is.null(mbiome)){
      mbiome = get("mbiome", envir = parent.frame())
      
    }

    if(sym == 1){
        mean_base = get(mean_base_sym, envir = parent.frame())

      }else{
        mean_base = get(mean_base, envir = parent.frame())

      }
  }

  scale_mbiome <- sweep(mbiome, 2, mean_base, "-")

  if(identical(pop@misc,list())){
    pop@misc = vector(length = nInd(pop), mode = "list")
    miscNULL = list(M = NA, M_sym = NA,mv = NA, mv_sym = NA)

    pop@misc <- lapply(1:nInd(pop), function(ind) {
      pop@misc[[ind]] = miscNULL
    })

  }


  mbiome <- data.frame(mbiome)
  
  if(sym==0){

    pop@misc <- lapply(1:nInd(pop), function(ind) {
      pop@misc[[ind]]$M <- mbiome[ind, ]
      pop@misc[[ind]]$mv <- (as.matrix(scale_mbiome[ind, ]) %*% w)
     # pop@misc[[ind]]$M_sym <- NA
    #  pop@misc[[ind]]$mv_sym <- NA
      pop@misc[[ind]]
    })
  }

  if(sym == 1){

    pop@misc <- lapply(1:nInd(pop), function(ind) {
     # pop@misc[[ind]]$M <- pop@misc[[ind]]$M
    #  pop@misc[[ind]]$mv <- pop@misc[[ind]]$mv
      pop@misc[[ind]]$M_sym <- mbiome[ind, ]
      pop@misc[[ind]]$mv_sym <- (as.matrix(scale_mbiome[ind, ]) %*% w)
      pop@misc[[ind]]
    })
  }

  return(pop)
}
