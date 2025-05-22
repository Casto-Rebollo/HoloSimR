#' Generate Microbiota Architecture
#'
#' @description
#' Generation of the architecture for the parental and environmental microbiota at founder population. The \pkg{\link{mobsim}} package is used to help simulate the species abundance distribution. Saving it as \code{founderM} will facilitate access by other functions.
#'
#' @param globalSP an object of \code{\link{GlobalSP}} with the global parameters of the simulation. If it is \code{NULL}, the globalSP will be extracted from the Global Environment
#' @usage setFounderM(globalSP = NULL)
#' @return a list including three dataframes:\enumerate{
#'   \item{\strong{architecture}: Dataframe with the architecture of the microbiota summarized in the following columns:\itemize{
#'        \item{\code{RA_EM}: Mean of the relative abundance of the environmental microbiota}
#'        \item{\code{Species}: Artificial name for microbial species}
#'        \item{\code{EM}: Mean of the real abundance of the environmental microbiota}
#'        \item{\code{PM0}: Vector indicating which microbial species will form the microbiota of the base population}
#'        \item{\code{PM}: Mean of the real abundance of the parental microbiota}
#'        \item{\code{GR_PM}: Mean of the growth rate of the microbial species when belonging to the parental microbiota}
#'        \item{\code{SD}: Standard deviation of the growth rate}
#'        \item{\code{RA_PM0}: Mean of the relative abundance of the microbial species in the parental microbiota}
#'        \item{\code{w}: Vector of the effect of species abundance on the trait}
#'        \item{\code{beta}: Vector of 0 or 1 indicating the species abundance influenced by the host genome}}}
#'    \item{\strong{beta}: Dataframe with the values of the genetic effect on the abundance of each species. Dimensions \eqn{nQTL} x \eqn{nSpecies}}
#'    \item{\strong{symbiosis}: Dataframe with the symbiosis effect among species. Dimensions \eqn{nSpecies} x \eqn{nSpecies}}}
#' @import mobsim
#' @importFrom stats runif rgamma
#' @export
#' @details
#' The simulation of \eqn{\beta} and \eqn{\omega} is based on \cite{Perez-Enciso et al. (2021)}.
#' \eqn{\beta} is the genetic effect of the QTLs on the species abundance, which follows a gamma distribution \eqn{\Gamma(shape = 0.2, scale=5)} as suggested in \cite{Caballero et al. (2015)}. \eqn{\omega} is the effect of the species abundance on the phenotype which follows a gamma distribution \eqn{\Gamma(shape = 1.4, scale=3.8)} (\cite{Pérez-Enciso et al., 2021}). The hyperparameters of this distribution consider that few microbial species have an effect on the phenotype (\cite{Duvallet et al., 2017}; \cite{Casto-Rebollo et al., 2023}}).
#'
#' @examples
#'
#' # Set simulation parameters
#' gSP <- GlobalSP$new(nPop = 1000, nQTLchr = 100, nyear = 10)
#'
#' # Defining the microbiome
#' gSP$setSpecies(nSpecies = 1000, nSp0 = 600, nSpEff = 35)
#'
#' # Making microbiome architecture
#' founderM <- setFounderM(globalSP = gSP)
#'
#' @references Pérez-Enciso M, Zingaretti LM, Ramayo-Caldas Y et al. Opportunities and limits of combining microbiome and genome data for complex trait prediction. Genet Sel Evol. 2021;53:65. \doi{10.1186/s12711-021-00658-7}
#' @references Caballero A, Tenesa A, Keightley PD. The nature of genetic variation for complex traits revealed by GWAS and regional heritability mapping analyses. Genetics. 2015;201:1601–13. \doi{10.1534/genetics.115.177220}
#' @references Duvallet C, Gibbons SM, Gurry T et al. Meta-analysis of gut microbiome studies identifies disease-specific and shared responses. Nat Commun 2017;8:1784. \doi{10.1038/s41467-017-01973-8}
setFounderM <- function(globalSP = NULL){

  if(is.null(globalSP)){
    globalSP = get("gSP", envir = .GlobalEnv)
  }

  simEM <- sim_sad(s_pool = globalSP$nSpecies,
                   n_sim = globalSP$getPrivate()$bioLimit,
                   sad_type = "nbinom",
                   sad_coef = list(size = 1, prob = 0.0001),
                   fix_s_sim = FALSE,
                   drop_zeros = T
  )

  globalSP$nSpecies <- length(simEM)

  # Stable parameters of simulate species at founder population

  founderM <- data.frame(RA_EM = simEM / sum(simEM),
                         Species = paste("Sp", 1:globalSP$nSpecies, sep = ""))
  founderM$EM <- sample(simEM)
  # Parental microbiome
  founderM$PM0 <- sample(c(rep(1, globalSP$nSp0),
                           rep(0, globalSP$nSpecies - globalSP$nSp0))) #Position of species which are present in the founder population

  ###################################################################
  founderM$PM <- sample(simEM) #Mean of species abundance in the parental microbiome

  founderM$GR_PM <- founderM$PM/log2(founderM$PM)

  founderM$SD <- round(sample(runif(globalSP$nSpecies, 0,2)) * founderM$GR_PM,2) #Standard deviation of species abundance

  founderM$RA_PM0[founderM$PM0 == 1] <- founderM$PM[founderM$PM0 == 1] / sum(founderM$PM[founderM$PM0 == 1])
  founderM$RA_PM0[is.na(founderM$RA_PM0)] <- 0

  #################################
  # Simulation of species effects #
  #################################

  founderM$w <- 0

  founderM$w[founderM$PM0==1] <- sample(c(rgamma(globalSP$nSpEff,
                                                 shape = 1.4,
                                                 scale = 3.8),
                                          rep(0, globalSP$nSp0 - globalSP$nSpEff)))*sample(c(-1,1),globalSP$nSp0, replace = T)

  #Force species with effect to be variable
  founderM$SD[founderM$w!=0] <- round(sample(runif(globalSP$nSpEff,0.3,0.6))*founderM$GR_PM[founderM$w!=0],2)

  writeLines("1. Microbiota architecture DONE")

  ######################
  # HOST GENOME EFFECT #
  ######################

  # QTL effect (host genotype) on microbiota (Microbial heritability)
  if("H" %in% globalSP$model || globalSP$model == "All"){
    if(globalSP$MH.H == NULL){

      common_factor <- (globalSP$meanMH * (1 - globalSP$meanMH) / globalSP$varMH) - 1
      alpha <- globalSP$meanMH * common_factor
      beta <- (1 - globalSP$meanMH) * common_factor

      globalSP$MH.H <- rbeta(globalSP$nSpecies, alpha, beta)

    }
    MH.nSp_index <- which(globalSP$MH.H != 0)

  }else{
    MH.w <- round(globalSP$nSpEff * globalSP$propMH.wSp)
    MH.all <- round(globalSP$nSpecies * globalSP$propMH)

    if(MH.all < MH.w){
     MH.w <- MH.all * globalSP$propMH.wSp
    }

    MH.wSp_index <- sample(which(founderM$w != 0), MH.w) #MH of species with w
    founderM$beta <- 0
    founderM$beta[MH.wSp_index] <- 1

    MH.nSp <- MH.all - MH.w #Rest of MH species

    MH_index <- sample(which(founderM$beta == 0 ), MH.nSp)
    founderM$beta[MH_index] <- 1

    # All species with microbial heritability
    MH.nSp_index <- which(founderM$beta == 1)
  }



  # Pre-allocate beta matrix
  beta <- matrix(0, nrow = globalSP$nChr * globalSP$nQTLchr, ncol = globalSP$nSpecies)

  # Sample QTL.MH_index directly without replacement
  QTL.MH_index <- sample.int(globalSP$nQTLchr * globalSP$nChr,
                             size = round(globalSP$nQTLchr * globalSP$nChr * globalSP$propQTL),
                             replace = FALSE)

  # Generate random numbers and samples
  random_numbers <- rgamma(length(QTL.MH_index), shape = 0.2, scale = 5)
  samples <- sample(c(1, -1), round(globalSP$nQTLchr * globalSP$nChr * globalSP$propQTL), replace = TRUE)

  # Assign values to beta matrix
  beta[QTL.MH_index, MH.nSp_index] <- matrix(random_numbers * samples, nrow = length(QTL.MH_index))

  beta <- as.data.frame(beta)
  writeLines("2. Host genetic effect on abundance DONE")
  #############
  # Symbiosis #
  #############

  # How species is affected by the others in columns
  # How species affect the others in rows
  # Initialize founderMxM matrix with zeros

  founderMxM <- matrix(0, nrow = globalSP$nSpecies, ncol = globalSP$nSpecies)

  # Create a data frame for upper triangular indices
  upper_tri_indices <- data.frame(which(upper.tri(founderMxM), arr.ind = TRUE))

  # Create a vector of interaction types based on percentages
  interaction_types <- c(rep("commensal", round(globalSP$nSpecies * globalSP$getPrivate()$propComm)),
                         rep("mutualism", round(globalSP$nSpecies * globalSP$getPrivate()$propMut)),
                         rep("ammensalism", round(globalSP$nSpecies * globalSP$getPrivate()$propAmm)),
                         rep("competition", round(globalSP$nSpecies * globalSP$getPrivate()$propComp)),
                         rep("exploitation", round(globalSP$nSpecies * globalSP$getPrivate()$propExp)))

  # Sample interaction types
  interaction_types <- sample(interaction_types, size = nrow(upper_tri_indices), replace = TRUE)

  # Precompute interaction vectors
  interaction_vectors <- sapply(interaction_types, function(type) {
    switch(
      type,
      commensal = c(1, 0),
      mutualism = c(1, 1),
      competition = c(-1, -1),
      ammensalism = c(0, -1),
      exploitation = c(1, -1)
    )
  })

  #Loop
  for (i in seq_along(interaction_types)) {
    indices <- upper_tri_indices[i, ]

    # Assign interaction types to founderMxM matrix for both rows and columns
    founderMxM[indices$row, indices$col] <- interaction_vectors[1, i]
    founderMxM[indices$col, indices$row] <- interaction_vectors[2, i]
  }

  # Set the diagonal elements to 0
  diag(founderMxM) <- 0

  # Effect of species interaction
  founderMxM <- founderMxM * matrix(rgamma(globalSP$nSpecies^2, shape = 0.2, scale = 5), nrow = globalSP$nSpecies)
  founderMxM <- as.data.frame(founderMxM)

  writeLines("3. Symbiosis DONE")

  return(list("architecture" = founderM, "beta" = beta, "symbiosis" = founderMxM))

}
