#' Generate Founder Population Demography
#'
#' @description
#' This function calls \code{\link{runMacs}} or \code{\link{runMacs2}} from the \pkg{\link{AlphaSimR}} package to generate the founder haplotypes. For more details, refer to the documentation for those functions.
#' This function allows the implementation of the rabbit's demographic history using \code{\link{runMacs2}}, based on \cite{Carneiro et al. (2014)}. Saving this object as \var{founderGenomes} will enable access by default in other functions.
#'
#' @param globalSP An object of class \code{\link{GlobalSP}} containing the global parameters of the simulation. If assigned as \code{gSP}, the function will access it directly.
#' @param rndSeed A seed for random sampling. If it is \code{NULL}, a random seed will be set.
#' @usage setFounderG(globalSP = NULL, rndSeed = NULL)
#' @return An object of class \code{\link{MapPop-class}}.
#' @export
#'
#' @examples
#' # Set global simulation parameters
#' gSP <- GlobalSP$new(nPop = 1000, nyear = 10, nQTLchr = 100)
#'
#' # Define the microbiota
#' gSP$setSpecies(nSpecies = 1000, nSp0 = 600, nSpEff = 100, symbiosis = c(0,1))
#'
#' # Define the trait
#' gSP$setTrait(meanP = 10, varP = 3.8, h2 = 0.15, m2 = 0.15)
#'
#' # Create founder genetics
#' founderG <- setFounderG(globalSP = gSP)
#'
#' @references Carneiro M, Albert FW, Afonso S, Pereira RJ, Burbano H, Campos R, et al. The genomic architecture of population divergence between subspecies of the European rabbit. PLoS Genet. 2014;10:e1003519
#'
setFounderG <- function(globalSP = NULL,
                        rndSeed=NULL){

  if(is.null(globalSP)){
    globalSP = get("gSP", envir = .GlobalEnv)
  }

  if(!globalSP$animal%in%c("RABBIT","GENERIC","CATTLE")){
    stop("The animal species indicated is wrong. Please, check that it match with the incorporated in this package: 'RABBIT', 'CATTLE', 'GENERIC'")
  }


  if(is.null(rndSeed)){
    rndSeed <- sample(c(1000:2000),1)
  }

  if(globalSP$animal == "RABBIT"){

    ## Genome
    BPchr <- 127543145
    CM <- 1.34
    mutR <- 1.74e-08

    ## Demography history of Oryctolagus cuniculus
    Ne <- 100
    histNE <- c(100, 125, 150, 175, 200, 225, 250, 275,
                300, 325, 400, 500, 600, 700, 800, 900,
                1000, 1100, 1200, 100000)

    histG <- c(0, 20, 30, 40, 50, 60, 70, 80, 90, 100,
               200, 400, 600, 800, 1000, 1200, 1400,
               1600, 1800, 2000)

    set.seed(rndSeed)
    founderGenomes = runMacs2(nInd = globalSP$nPop,
                              nChr = globalSP$nChr,
                              Ne = Ne,
                              segSites = globalSP$segSITESchr,
                              bp = BPchr,
                              genLen = CM,
                              mutRate = mutR,
                              histNe = histNE,
                              histGen = histG,
                              nThreads = globalSP$nt)

  }else{
    set.seed(rndSeed)
    founderGenomes = runMacs(nInd = globalSP$nPop,
                             nChr = globalSP$nChr,
                             segSites = globalSP$segSITESchr,
                             species = globalSP$animal,
                             nThreads = globalSP$nt)
  }

  return(founderGenomes)

}
