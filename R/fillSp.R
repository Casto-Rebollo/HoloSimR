#' Fill Species Information
#'
#' @description
#' This function incorporates microbial species abundances and the microbiota value into an object of class \code{\link{Pop-class}}.
#'
#' @param pop An object of class \code{\link{Pop-class}} generated using the \pkg{\link{AlphaSimR}} package, representing the current population.
#' @param mbiome A microbiota abundance matrix generated using the function \code{\link{makeM}}. If \code{NULL}, \code{mbiome} will be extracted from the global environment.
#' @param w A vector with the \eqn{\omega} effect scaled according to the base population.
#' @param sym A value of 1 to consider the symbiosis effect on the microbiota. \emph{Default} is 0.
#' @param mean Vector of microbiota mean abundance from base population
#'
#' @return An object of class \code{\link{Pop-class}} with updated microbiota information.
#' @export
#' @usage fillSp(pop, w = NULL, mbiome = NULL, sym = 0)
#'
fillSp <- function(pop, w = NULL,mbiome = NULL, sym = 0, mean = NULL) {

if (exists("basePop", envir = .GlobalEnv)) {
  basePop_env <- .GlobalEnv
} else if (exists("basePop", envir = parent.frame())) {
  basePop_env <- parent.frame()
} else {
  basePop_env <- NULL
}

if (!is.null(basePop_env)) {
  basePop <- get("basePop", envir = basePop_env)

  if (is.null(mbiome)) {
    if (exists("mbiome", envir = basePop_env)) {
      mbiome <- get("mbiome", envir = basePop_env)
    } else {
      stop("mbiome not found in basePop environment.")
    }
  }

  if (is.null(w)) {
    if (sym == 1) {
      w <- matrix(basePop[["w_scale1"]], ncol=1)
    } else {
      w <- matrix(basePop[["w_scale0"]], ncol=1)
    }
  }
  if(is.null(mean)){
    if (sym == 1) {
      mean <- basePop[["mean_base_sym"]]
    } else {
      mean <- basePop[["mean_base"]]
    }
  }

  } else {
  # If basePop is not found, get from parent.frame()
  if (is.null(mbiome)) {
    mbiome <- get("mbiome", envir = parent.frame())
  }
}

  w <- matrix(w, ncol=1, )
  scale_mbiome <- sweep(mbiome, 2, mean, "-")

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
      pop@misc[[ind]]$mv <- as.numeric(scale_mbiome[ind, , drop = FALSE] %*% w)
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
      pop@misc[[ind]]$mv_sym <- as.numeric(scale_mbiome[ind, , drop = FALSE] %*% w)
      pop@misc[[ind]]
    })
  }

  return(pop)
}
