#' @title Global Simulation Parameters
#'
#' @description
#' The \sQuote{GlobalSP} class represents global simulation parameters for a specific
#' simulation scenario. It includes general parameters, trait parameters,
#' and microbiota parameters.
#'
#' @field nPop (numeric)
#' Number of individuals in the founder population.
#' 
#' @field dataM (TXT file)
#' Real microbiota data (count data). TXT file separated by tabulations. Indicate the complete path if the file is not in the working directory
#'
#' @field nChr (numeric)
#' Total number of chromosomes to simulate using \code{runcMacs} from \code{\link{AlphaSimR}}.
#'
#' @field nQTLchr (numeric)
#' Number of QTL per chromosome affecting the simulated trait.
#'
#' @field segSITESchr (numeric)
#' Total number of segregating sites to simulate using \code{runcMacs} from \code{\link{AlphaSimR}}.
#' If it is \code{NULL}, all are retained.
#'
#' @field nt (numeric)
#' If OpenMP is available, this allows for simulating chromosomes in parallel. \emph{Default} is 1.
#'
#' @field nSim (numeric)
#' Number of repetitions of the simulation. \emph{Default} is 1.
#'
#' @field animal (character)
#' Type of animal. Currently, three species histories are available: \dQuote{GENERIC}, \dQuote{CATTLE},
#' and \dQuote{RABBIT}. All except the rabbit demographic histories are implemented in \pkg{\link{AlphaSimR}}.
#' The rabbit demographic is implemented in this package following \cite{Carneiro et al. (2014)}.
#'
#' @field selType (character)
#' Selection type (e.g., \dQuote{Divergent}, \dQuote{Low}, \dQuote{High}).
#'
#' @field model (character)
#' Simulation model. There are six models implemented:
#' \itemize{
#'  \item \dQuote{G}: Effect of the genome only
#'  \item \dQuote{M}: Effect of microbiota only
#'  \item \dQuote{NMH}: Hologenome without microbial heritability (Null Microbial heritability)
#'  \item \dQuote{LMH}: Hologenome considering Low Microbial Heritability
#'  \item \dQuote{MHM}: Hologenome considering Medium Microbial Heritability
#'  \item \dQuote{HMH}: Hologenome considering High Microbial Heritability
#'  \item \dQuote{H}: Hologenome considering variability of Microbial heritability
#' }
#' All models can run simultaneously starting from the same base population if \dQuote{All} is specified.
#'
#' @field nyear (numeric)
#' Number of years.
#'
#' @field nSire (numeric)
#' Number of sires.
#'
#' @field nDam (numeric)
#' Number of dams.
#'
#' @field nCross (numeric)
#' Number of crosses.
#'
#' @field nSon (numeric)
#' Number of progeny per cross.
#'
#' @field scale (character)
#' Type of transformation to perform in the data. \emph{Default} is log ()
#' Centered-log ratio transformation is also allowed ("clr")
#'
#' @field meanP (numeric)
#' Mean value for the simulated trait.
#'
#' @field varP (numeric)
#' Variance of the simulated trait.
#'
#' @field h2 (numeric)
#' Heritability of the simulated trait.
#'
#' @field m2 (numeric)
#' Microbiability for the simulated trait.
#'
#' @field nSpecies (numeric)
#' Total number of microbial species simulated in a community.
#'
#' @field nSp0 (numeric)
#' Total number of simulated microbial species in the base population.
#'
#' @field nSpEff (numeric)
#' Number of microbial species affecting the simulated trait.
#'
#' @field PM (numeric)
#' Proportion of the initial parental microbiota (dam microbiota) exposed to the offspring.
#'
#' @field EM (numeric)
#' Proportion of the initial environmental microbiota exposed to the offspring.
#'
#' @field propMH (numeric)
#' Proportion of microbial species affected by the host genome.
#'
#' @field propQTL (numeric)
#' Proportion of QTLs affecting the abundance of the microbial species.
#'
#' @field propMH.wSp (numeric)
#' Proportion of microbial species with effects on the simulated trait that are affected by the host genome.
#'
#' @field MH.G (numeric)
#' Microbial heritability in the G scenario.
#'
#' @field MH.M (numeric)
#' Microbial heritability in the M scenario.
#'
#' @field MH.low (numeric)
#' Microbial heritability in the LMH scenario.
#'
#' @field MH.medium (numeric)
#' Microbial heritability in the MHM scenario.
#'
#' @field MH.high (numeric)
#' Microbial heritability in the HMH scenario.
#'
#' @field MH.H (numeric vector)
#' Vector of variable microbial heritability
#'
#' @field meanMH (numeric)
#' Mean of microbial heritability distribution
#'
#' @field varMH (numeric)
#' Variance of microbial heritability
#'
#' @field symbiosis (numeric)
#' Value or vector to specify whether to simulate (1) or not (0) the microbial species interaction or symbiosis effect on microbial abundance.
#' There are five types of interactions simulated according to \cite{Coyte KZ and Rakoff-Nahoum S (2019)}: Commensalism, Mutualism, Ammensalism, Competition, and Exploitation/Pathogen.
#'
#' @field s2 (numeric)
#' Proportion of microbial species abundance due to the symbiosis or interaction between microbial species.
#'
#' @section Methods:
#'
#' \bold{\emph{Public methods}}
#'
#' \itemize{
#'  \item \code{GlobalSP$new()} Create a new GlobalSP object.
#'  \item \code{GlobalSP$setTrait()} Set trait parameters for the GlobalSP object.
#'  \item \code{GlobalSP$setSpecies()} Set species parameters for the GlobalSP object.
#'  \item \code{GlobalSP$getValues()} Retrieve values from the GlobalSP object.
#' }
#'
#' @aliases GlobalSP
#' @export
#' @examples
#' # Set simulation parameters
#' gSP <- GlobalSP$new(nPop = 1000, nQTLchr = 100, nyear = 10)
#'
#' # Defining the trait
#' gSP$setTrait(meanP = 8, varP = 2, h2 = 0.15, m2 = 0.15)
#'
#' # Defining the microbiota
#' gSP$setSpecies(nSpecies = 1000, nSp0 = 600, nSpEff = 35)
#'
#' @import R6
#' @export
#' @references Coyte KZ, Rakoff-Nahoum S. Understanding competition and cooperation within the mammalian gut microbiome. Curr Biol. 2019;29:R538-R544
GlobalSP <- R6Class(
  classname = "GlobalSP",
  public = list(

    # General parameters
    nPop = NULL,
    dataM = NULL,
    nChr = 1,
    nQTLchr = NULL,
    segSITESchr = NULL,
    animal = "GENERIC",
    selType = "Divergent",
    nyear = NULL,
    model = "NMH",
    nt = 1,
    nSim = 1,

    nSire = 25,
    nDam = 125,
    nCross = NULL,
    nSon = 8,
    scale = "log",

    # Trait parameters
    meanP = NULL,
    varP = NULL,
    h2 = NULL,
    m2 = NULL,

    # Microbiome parameters
    nSpecies = NULL,
    nSp0 = NULL,
    nSpEff = NULL,
    PM = 0.5,
    EM = 0.5,
    propMH = 0.1,
    propQTL = 0.1,
    propMH.wSp = 0.5,
    symbiosis = c(0, 1),
    s2 = 0.2,
    MH.G = 0.2,
    MH.M = 0.2,
    MH.low = 0.1,
    MH.medium = 0.3,
    MH.high = 0.6,
    MH.H = NULL,
    meanMH = NULL,
    varMH = NULL,

    #' @description
    #' Function to initialize the parameters needed for the simulation.
    #' Saving this object as \var{gSP} will facilitate the functionality of the other functions.
    #' @param nPop The total number of individuals to simulate in the founder population.
    #' @param dataM Path to TXT file format of microbiota composition (count data).
    #' @param nChr The total number of chromosomes to simulate using \code{runcMacs} from \code{\link{AlphaSimR}}. \emph{Default} is 1.
    #' @param nQTLchr A value or vector indicating the number of QTL per chromosome affecting the simulated trait.
    #' @param segSITESchr The total number of segregating sites per chromosome to simulate using \code{runcMacs} from \code{\link{AlphaSimR}}.
    #' If it is \var{NULL}, all are retained.
    #' @param animal The species to simulate. Can be \dQuote{GENERIC}, \dQuote{CATTLE}, or \dQuote{RABBIT}.
    #' All except the rabbit demographic histories are implemented in \pkg{\link{AlphaSimR}}.
    #' @param selType The selection type (e.g., \dQuote{Divergent}, \dQuote{Low}, \dQuote{High}).
    #' @param nyear The total number of years to simulate.
    #' @param model The simulation model.
    #' There are five models implemented: \dQuote{G}, \dQuote{M}, \dQuote{NMH}, \dQuote{LMH}, \dQuote{MHM}, and \dQuote{HMH}.
    #' @param nt If OpenMP is available, this allows for simulating chromosomes in parallel. \emph{Default} is 1.
    #' @param nSim The number of repetitions of the simulation. \emph{Default} is 1.
    #' @param nSire Number of breeding sires. \emph{Default} is 25.
    #' @param nDam Number of breeding dam. \emph{Default} is 125. 
    #' @param nCross Number of crosses per sire. \emph{Default} in \var{NULL}. If it is \var{NULL}, it will directly computed. 
    #' @param nSon Number of progeny per cross. \emph{Default} is 8.
    #' @param scale Transformation to perform in the data. \emph{Default} log() is taken.
    initialize = function(nPop, dataM = NULL, nChr = 1, nQTLchr, segSITESchr = NULL,
                          animal = "GENERIC", nSire = 25, nDam = 125,
                          nCross = NULL, selType = "Divergent", nyear, model = "NMH", nt = 1, nSim = 1,
                          scale = "log"){
      self$nPop <- nPop
      self$dataM <- dataM
      self$nChr <- nChr
      self$nQTLchr <- nQTLchr
      self$segSITESchr <- segSITESchr
      self$animal <- animal
      self$selType <- selType
      self$nyear <- nyear
      self$model <- model
      self$nt <- nt
      self$nSim <- nSim
      self$nSire <- nSire
      self$nDam <- nDam
      self$nCross <- nCross
      self$scale <- scale
    },

    #' @description
    #' Function to set trait parameters for the GlobalSP object.
    #' @param meanP The mean value for the trait.
    #' @param varP The variance of the trait.
    #' @param h2 Heritability of the trait.
    #' @param m2 Microbiability of the trait.
    setTrait = function(meanP, varP, h2, m2) {
      self$meanP <- meanP
      self$varP <- varP
      self$h2 <- h2
      self$m2 <- m2

      invisible(self)
    },

    #' @description
    #' Function to set species parameters for the GlobalSP object.
    #' @param nSpecies Total number of microbial species in the community. \emph{Default} is NULL
    #' @param nSp0 Number of initial microbial species. \emph{Default} is NULL
    #' @param nSpEff Number of effective microbial species affecting the trait.
    #' @param PM Proportion of parental microbiota exposed to offspring.
    #' @param EM Proportion of environmental microbiota exposed to offspring.
    #' @param propMH Proportion of species influenced by the host genome.
    #' @param propQTL Proportion of QTLs affecting the microbial species.
    #' @param propMH.wSp Proportion of species affecting the trait influenced by the genome.
    #' @param symbiosis Whether to simulate symbiotic effects (1) or not (0).
    #' @param s2 Proportion of microbial abundance explained by symbiosis.
    #' @param MH.G Microbial heritability in the G scenario.
    #' @param MH.M Microbial heritability in the M scenario.
    #' @param MH.low Microbial heritability in the LMH scenario.
    #' @param MH.medium Microbial heritability in the MHM scenario.
    #' @param MH.high Microbial heritability in the HMH scenario.
    #' @param MH.H Vector of microbial heritability in H scenario.
    #' @param meanMH Mean of microbial heritability in H scenario.
    #' @param varM Variance of microbial heritability in H scenario.
    setSpecies = function(nSpecies = NULL, nSp0 = NULL, nSpEff, PM = 0.5, EM = 0.5, propMH = 0.1,
                          propQTL = 0.1, propMH.wSp = 0.5, symbiosis = c(0, 1),
                          s2 = 0.2, MH.G = 0.2, MH.M = 0.2,
                          MH.low = 0.1, MH.medium = 0.3, MH.high = 0.6,
                          MH.H = NULL, meanMH = NULL, varMH = NULL) {
      self$nSpecies <- nSpecies
      self$nSp0 <- nSp0
      self$nSpEff <- nSpEff
      self$PM <- PM
      self$EM <- EM
      self$propMH <- propMH
      self$propQTL <- propQTL
      self$propMH.wSp <- propMH.wSp
      self$symbiosis <- symbiosis
      self$s2 <- s2
      self$MH.G <- MH.G
      self$MH.M <- MH.M
      self$MH.low <- MH.low
      self$MH.medium <- MH.medium
      self$MH.high <- MH.high
      self$MH.H <- MH.H
      self$meanMH <- meanMH
      self$varMH <- varMH

      invisible(self)
    },

    #' @description
    #' Function to retrieve the values of the GlobalSP object.
    getValues = function() {
      print(list(
        nPop = self$nPop,
        dataM = self$dataM,
        nChr = self$nChr,
        nQTLchr = self$nQTLchr,
        segSITESchr = self$segSITESchr,
        animal = self$animal,
        selType = self$selType,
        nyear = self$nyear,
        model = self$model,
        nt = self$nt,
        nSim = self$nSim,
        nSire = self$nSire,
        nDam = self$nDam,
        nCross = self$nCross,
        nSon = self$nSon,
        meanP = self$meanP,
        varP = self$varP,
        h2 = self$h2,
        m2 = self$m2,
        nSpecies = self$nSpecies,
        nSp0 = self$nSp0,
        nSpEff = self$nSpEff,
        PM = self$PM,
        EM = self$EM,
        propMH = self$propMH,
        propQTL = self$propQTL,
        propMH.wSp = self$propMH.wSp,
        symbiosis = self$symbiosis,
        s2 = self$s2,
        MH.G = self$MH.G,
        MH.M = self$MH.M,
        MH.low = self$MH.low,
        MH.medium = self$MH.medium,
        MH.high = self$MH.high,
        MH.H = self$MH.H,
        meanMH = self$meanMH,
        varMH = self$varMH
      ))
    },

    #' @description
    #' Function to retrieve private parameters.
    getPrivate = function() {
      return(list(
        limit.Trait = self$limit.Trait,
        nDam.Lit = 8,
        bioLimit = 10^8,
        propComm = 0.2,
        propMut = 0.2,
        propAmm = 0.2,
        propComp = 0.2,
        propExp = 0.2
      ))
    }
  )
)

