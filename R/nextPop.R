#' Generate the Next Population
#' @description
#' This function simulates animals for the next generation. A mating plan of class \code{\link{matingPLAN}} is required.
#'
#' @param pop An object of \code{\link{MapPop-class}}. If \code{NULL}, it will be obtained from the global environment (\code{Pop}).
#' @param selType Selection type (e.g., \dQuote{Divergent}, \dQuote{Low}, \dQuote{High}). If \code{NULL}, it will be obtained from the global environment.
#' @param nSon Number of offspring per cross. If \code{NULL}, it will be obtained from the global environment.
#' @param crossPlan A \code{list} or \code{data.frame} of type \code{\link{matingPLAN}} with the crosses to perform. If \code{NULL}, it will be obtained from the global environment.
#' @param LS \code{TRUE} if the phenotype being evaluated is litter size (LS). Default is \code{FALSE}.
#' @param globalSP An object of class \code{\link{GlobalSP}}. If \code{NULL}, it will be obtained from the global environment.
#'
#' @return A list including two objects of \code{\link{MapPop-class}} for the population with high (\dQuote{PopHigh}) and low (\dQuote{PopLow}) phenotypes.
#'
#' @export
#' @import AlphaSimR
#'
#' @usage
#' nextPop(pop,
#'         crossPlan,
#'         selType = NULL,
#'         nSon = NULL,
#'         LS = FALSE,
#'         globalSP = NULL)
#'
#' @examples
#'
#' # Set simulation global parameters
#' gSP = GlobalSP$new(nPop = 1000, nyear = 10, nQTLchr = 100)
#'
#' # Define the microbiota
#' gSP$setSpecies(nSpecies = 1000, nSp0 = 600, nSpEff = 100, symbiosis = c(0,1))
#'
#' # Define trait
#' gSP$setTrait(meanP = 10, varP = 3, h2 = 0.10, m2 = 0.10)
#'
#' # Create founder genetic
#' founderG = setFounderG(globalSP = gSP)
#'
#' # Set simulation genetic parameters
#' SP <- essentialSP(founder = founderG,
#'                   minSnpFreq = 0.05, nSnpChr = 10000)
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
#' # Compute phenotype in a base population
#' Pop <- makeP(pop = basePop$Pop, model = "NMH", sym = 1, limTrait = c(0,18), sex = "F")
#'
#' # Select breeding animals from basePop
#' parent <- selectBreeding(pop = Pop, nDam = 125, nSire = 25, sex = "F", g0 = TRUE)
#'
#' # Create mating plan for next generation
#' crossPlan <- matingPLAN(parent = parent)
#'
#' # Generate new generation
#' Pop <- nextPop(pop = Pop, crossPlan = crossPlan, nSon = 8, LS = TRUE)
nextPop <- function(pop, crossPlan,
                    selType = NULL, nSon = NULL,
                    LS = FALSE ,globalSP = NULL){

  if(is.null(globalSP)){
    globalSP = get("gSP", envir = .GlobalEnv)
  }

  if(is.null(nSon)){
    nSon = globalSP$nSon
  }

  newPop <- list(PopLow = NULL, PopHigh = NULL)

  if(is.null(selType)){
    selType = globalSP$selType
  }

  selOrder <- switch(
    selType,
    Divergent = c(FALSE,TRUE),
    Low = FALSE,
    High = TRUE
  )


  for(value in selOrder){

    Poptmp <- if (inherits(pop,"list")) {
      if (!value) pop$PopLow else pop$PopHigh
    } else {
      pop
    }

    matePlan <- if (!value) crossPlan$PopLow else crossPlan$PopHigh

    if(LS == TRUE){
      nSon = round(max(matePlan$Sel_value,na.rm = T))
    }

    if(is.null(globalSP$nCross)){
      gSP$nCross <<- globalSP$nDam/globalSP$nSire
      globalSP$nCross <- globalSP$nDam/globalSP$nSire
    }

    Pop <- makeCross(Poptmp, as.matrix(matePlan[,1:2]), nProgeny = nSon)

    ##Adapt population for litter size
    index_all <- NULL

    if(LS == TRUE){
      for (mom in matePlan$Dam) {
        nlit <- matePlan$Sel_value[matePlan$Dam == mom]
        index_all <- c(index_all, sample(which(Pop@mother == mom)))
      }

      Pop <- Pop[index_all, ]

    }

    if (!value){
      newPop$PopLow <- Pop
    }
    if(value){
      newPop$PopHigh <- Pop
    }


  }

  return(newPop)
}
