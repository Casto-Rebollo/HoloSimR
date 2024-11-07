#' Create Mating Plan
#' @description
#' Creates crosses based on the selected breeding animals using \code{\link{selectBreeding}}.
#'
#' Selection of breeding animals is according to the Phenotype (\dQuote{phe}) or the direct genetic value (\dQuote{gv}).
#' @param parent List of selected breeding animals obtained using the function \code{\link{selectBreeding}}.
#' @param nCross Number of females mating with each male. If \code{NULL}, it will be obtained from the global environment.
#' @param selType Selection type (e.g., \dQuote{Divergent}, \dQuote{Low}, \dQuote{High}). If \code{NULL}, it will be obtained from the global environment.
#' @param globalSP An object of class \code{\link{GlobalSP}}. If \code{NULL}, it will be obtained from the global environment.
#'
#' @return A list including two dataframes with the mating plans for increasing and decreasing the phenotype:\enumerate{
#'   \item{\strong{PopLow}}: Mating plan for decreasing the phenotype. ID of breeding females (Dam) and males (Sire), and value for the genetic or phenotypic value (Sel_value).
#'   \item{\strong{PopHigh}}: Mating plan for increasing the phenotype. ID of breeding females (Dam) and males (Sire), and value for the genetic or phenotypic value (Sel_value).
#'   }
#' @export
#'
#' @usage
#' matingPLAN(parent,
#'            nCross = NULL,
#'            selType = NULL,
#'            globalSP = NULL)
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
matingPLAN <- function(parent, nCross = NULL,
                       selType = NULL, globalSP = NULL){


  if(is.null(globalSP)){
    globalSP = get("gSP", envir = .GlobalEnv)
  }
  if(is.null(nCross)){
    nCross = globalSP$nCross
  }

  if(is.null(selType)){
    selType <- globalSP$selType
  }

  crossTotal <- list(PopLow = NULL, PopHigh = NULL)

  selOrder <- switch(
    selType,
    Divergent = c(FALSE,TRUE),
    Low = FALSE,
    High = TRUE
  )

  for(value in selOrder){
    if(value == FALSE){
      nDam <- nrow(parent$Dam_low)
      crossPlan <- data.frame(Dam = character(nDam), Sire = character(nDam))
      crossPlan$Dam <- parent$Dam_low$Dam
      crossPlan$Sel_value <- as.numeric(parent$Dam_low$Sel_value)
      male <- parent$Male_low
      female <- parent$Dam_low
    }else{
      nDam <- nrow(parent$Dam_high)
      crossPlan <- data.frame(Dam = character(nDam), Sire = character(nDam))
      crossPlan$Dam <- parent$Dam_high$Dam
      crossPlan$Sel_value <- as.numeric(parent$Dam_high$Sel_value)
      male <- parent$Male_high
      female <- parent$Dam_high
    }
    while(length(which(crossPlan$Sire==""))>0){

      cross.made <- 0
      for(m in 1:nrow(male)){

        sire <- male$father[m]
        dam <- male$mother[m]

        indx <- which(female$father != sire)
        indx <- indx[!indx %in% cross.made]

        if(length(indx) > nCross){
          indx <- sample(indx, nCross)
        }

        crossPlan$Sire[indx] <- male$Male[m]
        cross.made <- c(cross.made,indx)
      }
    }

    if(!value){
      crossTotal$PopLow <- crossPlan
    }
    if(value){
      crossTotal$PopHigh <- crossPlan
    }
  }
  return(crossTotal)
}

