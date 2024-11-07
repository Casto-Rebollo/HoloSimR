#' Compute Phenotypes
#'
#' @description
#' This function computes the phenotype of a population according to a specified scenario and generation.
#'
#' @param pop An object of \code{\link{MapPop-class}} for computing the phenotype.
#' @param model A string indicating the model to simulate. Only the values \dQuote{G}, \dQuote{M}, \dQuote{NMH}, \dQuote{LMH}, \dQuote{MMH}, or \dQuote{HMH} are available.
#' @param sym A value of 1 to consider the symbiosis effect on microbiota. \emph{Default} is 0.
#' @param limTrait A vector with the minimum and maximum values for the trait. If \code{NULL}, it will be obtained from the global environment.
#' @param sex Required if the phenotype is attributed to only one sex (e.g., litter size). Use \dQuote{M} for male and \dQuote{F} for female. Default is \dQuote{both} for both sexes.
#' @param globalSP An object of class \code{\link{GlobalSP}}. If \code{NULL}, it will be obtained from the global environment.
#' @param rndSeed A seed for random sampling. If it is \code{NULL}, a random seed will be set.
#'
#' @return An object of \code{\link{MapPop-class}} that includes the computed phenotypes.
#' @export
#' @importFrom stats rnorm
#' @usage
#' makeP(pop,
#'       model,
#'       sym = 0,
#'       limTrait = NULL,
#'       sex = "both",
#'       globalSP = NULL,
#'       rndSeed = NULL)
#'
#' @examples
#' # Set simulation global parameters
#' gSP <- GlobalSP$new(nPop = 1000, nyear = 10, nQTLchr = 100)
#'
#' # Define the microbiota
#' gSP$setSpecies(nSpecies = 1000, nSp0 = 600, nSpEff = 100, symbiosis = c(0, 1))
#'
#' # Define the trait
#' gSP$setTrait(meanP = 10, varP = 3, h2 = 0.10, m2 = 0.10)
#'
#' # Create founder genetic
#' founderG <- setFounderG(globalSP = gSP)
#'
#' # Set simulation genetic parameters
#' SP <- essentialSP(founder = founderG, minSnpFreq = 0.05, nSnpChr = 10000)
#'
#' # Create founder microbiota
#' founderM <- setFounderM(globalSP = gSP)
#'
#' # Create the founder population
#' founderPop <- newPop(founderG)
#'
#' # Create the base population
#' basePop <- simBasePop(model = "NMH")
#'
#' # Compute phenotype of the base population
#' Pop <- makeP(pop = basePop$Pop, model = "NMH", sym = 1, limTrait = c(0, 18), sex = "F")

makeP <- function(pop, model, sym = 0,
                  limTrait = NULL,  sex = "both",
                  globalSP = NULL, rndSeed = NULL){


  if(is.null(globalSP)){
    globalSP = get("gSP", envir = .GlobalEnv)
  }

  if(is.null(rndSeed)){
    rndSeed <- sample(c(1000:2000),1)
  }

  if(is.null(limTrait)){
    limTrait = globalSP$getPrivate()$limit.Trait
  }

  SP = get("SP", envir = .GlobalEnv)

  vE <- switch(
    model,
    G = globalSP$varP*(1-globalSP$h2),
    M = globalSP$varP*(1-globalSP$m2),
    NMH = globalSP$varP*(1-globalSP$m2 - globalSP$h2),
    LMH = globalSP$varP*(1-globalSP$m2 - globalSP$h2),
    MMH = globalSP$varP*(1-globalSP$m2 - globalSP$h2),
    HMH = globalSP$varP*(1-globalSP$m2 - globalSP$h2)
  )

  if(model != "G"){
    if(sym == 1){
      pop_mv <- unlist(lapply(pop@misc, `[[`, "mv_sym"))
    }else{
      pop_mv <- unlist(lapply(pop@misc, `[[`, "mv"))
    }
  }else{
    pop_mv <- rep(0,nInd(pop))
  }

  indx <- switch(
    sex,
    F = which(pop@sex=="F"),
    M = which(pop@sex=="M"),
    both = which(pop@sex == "F" | pop@sex == "M"),
  )


  if(model == "M"){
    set.seed(rndSeed)
    pop@pheno[indx] <- round(globalSP$meanP + pop_mv[indx] + rnorm(length(indx),sd=sqrt(vE)),2)
  }else{
    set.seed(rndSeed)
    pop@pheno[indx] <- round(globalSP$meanP + pop@gv[indx] + pop_mv[indx] + rnorm(length(indx),sd=sqrt(vE)),2)
  }


  if(!is.null(limTrait)){
    pop@pheno[pop@pheno > limTrait[2]] <- limTrait[2]
    pop@pheno[pop@pheno < limTrait[1]] <- limTrait[1]
  }

  return(pop)
}
