#' Simulation base population
#'
#' @description
#' Simulation of the base population according to the parameters specified in the \code{\link{GlobalSP}} and the genome and microbiota created for the founder population with \code{\link{setFounderG}} and \code{\link{setFounderM}}. This function sets the base population, scales the host genetic (\eqn{\beta}), symbiosis (\eqn{\gamma}) and microbiota (\eqn{\omega}) effects according to the variance specified by the user, and simulates the microbiota in the base population.
#'
#' @param model a value for the model to simulate. Only one value is available according to \dQuote{G}, \dQuote{M}, \dQuote{NMH}, \dQuote{LMH}, \dQuote{MMH}, \dQuote{HMH}.
#' @param founderPop Founder population to use for founding the base population. If it is \code{NULL}, it will be taken from the environment.
#' @param founderM a list of \code{\link{setFounderM}}. If \code{NULL}, it will be retrieved from the global environment.
#' @param globalSP an object of class \code{\link{GlobalSP}}. If \code{NULL}, it will be retrieved from the global environment.
#' @param rndSeed rndSeed seed to initialise the random sampling. If it is \code{NULL}, it will be set randomly.
#' @param progressBar progressBar If TRUE, progress bar and text will be displayed. \emph{Default} \code{FALSE}
#'
#' @return a list including three dataframes:\enumerate{
#' \item{\strong{Pop}}: an object of \code{\link{MapPop-class}} the base population
#' \item{\strong{Beta_scale}}: a matrix with the scaled \eqn{\beta} effect
#' \item{\strong{Symbiosis_scale}}: a matrix with the scaled \eqn{gamma} effect. It is \code{NULL} if there is no symbiosis.
#' \item{\strong{varG}}: a vector with the genetic variance of each microbial species.
#' \item{\strong{varS}}: a vector with the symbiosis variance of each microbial species. It is \code{NULL} if there is no symbiosis.
#' \item{\strong{varE}}: a vector with the environmental variance of each microbial species.
#' \item{\strong{w_scale0}}: a vector with the scaled \eqn{\omega} effect without taking the symbiosis effect into account.
#' \item{\strong{w_scale1}}: a vector with the scaled \eqn{\omega} effect including the symbiosis effect. It is \code{NULL} if there is no symbiosis.
#' }
#'
#' @export
#' @importFrom stats var
#' @import copula
#' @usage
#' simBasePop(model,
#'            founderPop = NULL,
#'            founderM = NULL,
#'            globalSP = NULL,
#'            rndSeed=NULL,
#'            progressBar = FALSE)
#' @details
#' The \eqn{\beta} (host genetic) and \eqn{\gamma} (symbiosis) effect on the microbiota were scaled according to the variance set by the user in the object \code{\link{GlobalSP}}.
#' The scale values were obtained by applying a scale factor follows \eqn{\sqrt{\frac{\text{var}_{\text{set}}}{\text{var}_{\text{sim}}}}} where \eqn{\text{var}_{\text{set}}} is the variance set according the parameters on \code{\link{GlobalSP}}, and \eqn{\text{var}_{\text{set}}} is the variance on the simulated matrix (\code{\link{setFounderM}}), to be adjusted.
#'
#'
#' @examples
#' #Set simulation global parameters
#' gSP = GlobalSP$new(nPop = 1000, nyear = 10, nQTLchr = 100)
#'
#' #Defining the microbiota
#' gSP$setSpecies(nSpecies = 1000, nSp0 = 600, nSpEff = 100, symbiosis = c(0,1))
#'
#' #Defining trait
#' gSP$setTrait(meanP = 10,varP=3,h2=0.10,m2=0.10)
#'
#' #Create founder genetic
#' founderG = setFounderG(globalSP = gSP)
#'
#' #Set simulation genetic parameters
#' SP <- essentialSP(founder = founderG,
#'                   minSnpFreq = 0.05,rndSeed=NULL, nSnpChr = 10000)
#'
#' #Create founder microbiota
#' founderM <- setFounderM(globalSP = gSP)
#'
#' #Create the founder population
#' founderPop <- newPop(founderG)
#'
#' #Create the base population
#' basePop <- simBasePop(model = "NMH")
#'
simBasePop <- function(model, founderPop = NULL,
                       founderM = NULL, globalSP = NULL, rndSeed=NULL, progressBar = FALSE){

  if(is.null(globalSP)){
    globalSP <- get("gSP", envir = .GlobalEnv)
  }

  if(is.null(globalSP$nCross)){
    gSP$nCross <<- floor(globalSP$nDam / globalSP$nSire)
    globalSP$nCross <- floor(globalSP$nDam / globalSP$nSire)
  }

  if(is.null(founderPop)){
    founderPop <- get("founderPop", envir = .GlobalEnv)
  }

  if(is.null(founderM)){
    founderM <- get("founderM", envir = .GlobalEnv)
    founderMxM <- founderM[["symbiosis"]]
    beta <- founderM[["beta"]]
    founderM <- founderM[["architecture"]]
  }

  if(is.null(rndSeed)){
    rndSeed <- sample(c(1000:2000),1)
  }

  if(model == "G"){
    founderM$w <- 0
  }

  SP = get("SP", envir = .GlobalEnv)

  set.seed(rndSeed)
  Pop <- randCross(founderPop,
                       nCrosses = globalSP$nDam,
                       nProgeny = globalSP$nSon)
  if(progressBar == TRUE){
    writeLines("--> Base Population DONE")
  }


  acquiredSp <- acquiredSpecies(pop = Pop, rndSeed = rndSeed)


  if(progressBar == TRUE){
    writeLines("--> Acquisition matrix DONE")
  }


  nPop <- nInd(Pop)

  if(is.null(model)){
    stop("Please provide a value for model. Vector is not supported")
  }

  h <- switch(
    model,
    G = globalSP$MH.G,
    M = globalSP$MH.M,
    NMH = 0,
    LMH = globalSP$MH.low,
    MMH = globalSP$MH.medium,
    HMH = globalSP$MH.high,
    H = globalSP$MH.H
  )

  varG.sp <- (founderM$SD^2) * h
  varE.sp <- (founderM$SD^2) - varG.sp

  set.seed(rndSeed)

  #Extract QTL to compute abundance due to host genome
  geno <- pullQtlGeno(Pop,
                      trait = 1,
                      chr = NULL,
                      simParam = SP)

  # Step 1: Calculate allele frequencies (p) for each SNP
  p <- colMeans(geno) / 2

  # Step 2: Center the matrix by subtracting 2 * p
  geno <- sweep(geno, 2, 2 * p)

  #geno <- scale(geno, center = TRUE, scale = FALSE)

  #Scale the genetic effect on the microbiome abundance
  scaleGxM <- geno %*% as.matrix(beta)
  varGxM <- apply(scaleGxM,2,var)

  baseGxM.scaled <- beta
  baseGxM.scaled <- sweep(beta, 2, sqrt(varG.sp / varGxM), `*`)

  baseGxM.scaled[is.na(baseGxM.scaled)] <- 0
  geno.biome <- geno %*% as.matrix(baseGxM.scaled)


  if(progressBar == TRUE){
    writeLines("--> Scaled beta matrix DONE")
  }


  #Compute the microbiome without symbiosis
  set.seed(rndSeed)
  
  #Environment

  mbiome_VE <- mvrnorm(nInd(Pop),mu = rep(0, globalSP$nSpecies) ,Sigma = diag(c(varE.sp)))
  
  mbiome <- matrix((acquiredSp + geno.biome + mbiome_VE),
                   nrow = nInd(Pop), ncol = length(founderM$w),
                   dimnames = list(NULL, founderM$Species))
  
  #Avoid extreme values
  quantile_99 <- apply(mbiome, 2, quantile, probs = 0.99)
  # Step 2: Truncate values above the 99th percentile
  mbiome <- sweep(mbiome, 2, quantile_99, FUN = pmin)
  
  #mbiome <- sweep(mbiome, 2, founderM$PM, `+`)
  
  #Scale microbiota without symbiosis
  scale_mbiome <- scale(mbiome, center = TRUE, scale = FALSE)
  mean_base <- attr(scale(mbiome, center = TRUE, scale = FALSE), "scaled:center")
  
  mv.raw <- scale_mbiome %*% founderM$w
  wScale0 <- as.vector(sqrt((globalSP$m2*globalSP$varP) / var(mv.raw))) * founderM$w

  mv.base <- mean(scale_mbiome %*% wScale0)

  symbiosis <- 0
  pop <- fillSp(pop=Pop, w = wScale0, mbiome = mbiome, sym = symbiosis, mean = mean_base)

  #Check for symbiosis effect
  if(1%in%globalSP$symbiosis){
    #Scale the interaction matrix according to the base population microbiota abundance
    
    if(is.null(globalSP$MH.H)){
      if((globalSP$s2+h) >1){
        stop("Please check the values for s2 and h2. The must sum 1 as maximun")
      }
    }
    
    if(!is.null(globalSP$MH.H)){
       varI.sp <- varE.sp * globalSP$s2
       varE_sym.sp <- varE.sp - varI.sp
       total <- varG.sp + varI.sp + varE_sym.sp
       if(total[1] > (founderM$SD[1]^2)[1]){
        print(total[1], (founderM$SD^2)[1])
        stop("Please check the values for s2. Remember that simulating H scenario, s is a proportion of VE")
       }
      
    }

    varI.sp <- varE.sp * globalSP$s2
    varE_sym.sp <- varE.sp - varI.sp

    set.seed(rndSeed)
    mbiome_VE <- mvrnorm(nInd(pop),mu = rep(0, globalSP$nSpecies), Sigma = diag(c(varE_sym.sp)))

    mbiome_sym <- matrix((geno.biome + mbiome_VE),
                     nrow = nInd(pop), ncol = length(founderM$w),
                     dimnames = list(NULL, founderM$Species))

    # Implement multivariate gamma with copula
    gauss_cop <- normalCopula(param = P2p(founderMxM), 
                              dim = globalSP$nSpecies, dispstr = "un")

    # Generate correlated uniforms via Gaussian copula
    set.seed(rndSeed)
    u <- copula::rCopula(globalSP$nSpecies, gauss_cop)
    
    gamma.sym <- NULL
    for(i in 1:globalSP$nSpecies){
      gamma.sym <- cbind(gamma.sym,
        qgamma(u[,i], shape=0.2,rate=5)
      )
    }

    #Generate matrix of abscence/presence of specie in the individual    
    mbiome_acq <- ifelse(mbiome_sym<0,0,1)
    
    scale_mbiome <- scale(mbiome_acq, center = TRUE, scale = FALSE)
    scaleMxM <- scale_mbiome %*% as.matrix(gamma.sym)
    varMxM <- apply(scaleMxM, 2, var)
    
    baseMxM.scaled <- gamma.sym
    baseMxM.scaled <- sweep(gamma.sym, 2, sqrt(varI.sp / varMxM), `*`)

    mbiome.sym <- scale_mbiome%*%baseMxM.scaled          

    #Scale microbiota with symbiosis
    mbiome_total <- matrix((acquiredSp + mbiome_sym + mbiome.sym),
                         nrow = nInd(pop), ncol = length(founderM$w),
                         dimnames = list(NULL, founderM$Species))

    #Avoid extreme values
    quantile_99 <- apply(mbiome_total, 2, quantile, probs = 0.99)
    # Step 2: Truncate values above the 99th percentile
    mbiome_total <- sweep(mbiome_total, 2, quantile_99, FUN = pmin)

    #mbiome_total <- sweep(mbiome_total, 2, founderM$PM, `+`)

    scale_mbiome <- scale(mbiome_total, center = TRUE, scale = FALSE)
    mean_base_sym <- attr(scale(mbiome_total, center = TRUE, scale = FALSE), "scaled:center")
    
    mv.raw <- scale_mbiome %*% founderM$w
    wScale <- as.vector(sqrt((globalSP$m2*c(globalSP$varP)) / var(mv.raw))) * founderM$w

    mv.base_sym <- mean(scale_mbiome %*% wScale)

    symbiosis <- 1
    pop <- fillSp(pop=pop, mbiome = mbiome_total, w = wScale, sym = symbiosis, mean = mean_base_sym)

    if(progressBar == TRUE){
      writeLines("--> Scaled symbiosis matrix DONE")
    }


  }else{
    baseMxM.scaled = NULL
    varE_sym.sp = NULL
    wScale = NULL
    mv.base_sym = NULL
    varI.sp = NULL
  }

  if(progressBar == TRUE){
    writeLines("--> Scale microbiota effect on the trait DONE")
  }


  ## Centered the mean of the microbiome values to 0 for base pop

  if(0%in%globalSP$symbiosis){
    pop_mv <- lapply(pop@misc, `[[`, "mv")
    pop_mv <- lapply(pop_mv, function(x) x)
    pop@misc <- Map(function(x, y) {x$mv <- y; x}, pop@misc, pop_mv)

  }

  if(1%in%globalSP$symbiosis){
    pop_mv <- lapply(pop@misc, `[[`, "mv_sym")
    pop_mv <- lapply(pop_mv, function(x) x)
    pop@misc <- Map(function(x, y) {x$mv_sym <- y; x}, pop@misc, pop_mv)
  }


  return(list("Pop" = pop,"Beta_scale" = baseGxM.scaled,"Symbiosis_scale" = baseMxM.scaled,
              "varG" = varG.sp , "varS" = varI.sp, "varE" = varE.sp,"varE_sym" = varE_sym.sp,
              "w_scale0"=wScale0,"w_scale1"=wScale,"mu_mv"=mv.base,"mu_mv_sym"=mv.base_sym,
              "mean_base_sym"= mean_base_sym, "mean_base" = mean_base))
}
