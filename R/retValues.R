#' Retrieve Simulation Values
#'
#' @description
#' Retrieve information from the simulation process, including:
#' \enumerate{
#'   \item{\strong{Generation}}: Computed generation
#'   \item{\strong{Model}}: Computed model
#'   \item{\strong{iter}}: Iteration performed; relevant if the selection is replicated
#'   \item{\strong{P}}: Average phenotype of the population
#'   \item{\strong{varP}}: Phenotypic variance
#'   \item{\strong{gv}}: Average additive genetic value
#'   \item{\strong{varG}}: Genetic variance
#'   \item{\strong{S_gv}}: Mean genetic value of the breeding animals
#'   \item{\strong{mv}}: Average microbiota value
#'   \item{\strong{varM}}: Variance of microbiota value
#'   \item{\strong{S_mv}}: Mean microbiota value of the breeding animals
#'   }
#'
#'
#' @param model The model value to simulate. Only values corresponding to \dQuote{G}, \dQuote{M}, \dQuote{NMH}, \dQuote{LMH}, \dQuote{MMH}, or \dQuote{HMH} are available.
#' @param g Current generation of selection. A value is required.
#' @param sym A value of 1 to consider the symbiosis effect on microbiota. \emph{Default} is 0.
#' @param selType Selection type (\dQuote{Low}, \dQuote{High}). A value is required.
#' @param iter Current iteration. \emph{Default} is 1.
#' @param pop An object of class \code{\link{MapPop-class}}. If \code{NULL}, it will be obtained from the environment (\code{Pop}).
#' @param breedInd A list of selected breeding animals obtained using the function \code{\link{selectBreeding}}. If \code{NULL}, it will be obtained from the environment (\code{parent}). A list of type \code{selectBreeding} can be used.
#'
#' @return A list containing two sub-lists with data frames that include the above information for the populations with high (\dQuote{PopHigh}) and low (\dQuote{PopLow}) phenotypes. Each population has a data frame for values with (\dQuote{Symbiosis}) and without (\dQuote{No_symbiosis}) the symbiosis effect.
#' @export
#' @import AlphaSimR
#' @import dplyr
#' @importFrom stats var
#'
#' @usage
#' retValues(model,
#'           g,
#'           sym = 0,
#'           selType,
#'           iter = 1,
#'           pop = NULL,
#'           breedInd = NULL)
#'
#' @examples
#'
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
#' #Computing phenotype in a base population
#' Pop <- makeP(pop = basePop$Pop, model = "NMH", sym = 0, limTrait = c(0,18), sex = "F")
#'
#' #Selection breeding animals from basePop
#' parent <- selectBreeding(pop = Pop, nDam = 125, nSire = 25, sex = "F", g0 = TRUE)
#'
#' #Retrieving values
#' report <- retValues(model = "NMH", g = 1, sym = 0, selType = "Low", pop=Pop, breedInd= parent)

retValues <- function(model, g, sym = 0, selType,iter = 1,
                      pop = NULL, breedInd = NULL){

  if(is.null(pop)){
    if(exists("Pop", envir = parent.frame())){
      pop <- get("Pop", envir = parent.frame())
    }else{
      pop <- get("Pop", envir = .GlobalEnv)
    }
  }

  if(is.null(breedInd)){
    if(exists("parent", envir = parent.frame())){
      breedInd <- get("parent", envir = parent.frame())
    }else{
      breedInd <- get("parent", envir = .GlobalEnv)
    }
  }


  #Check if there is a data.frame report generated
  if(exists("report", envir = parent.frame())){
    trend <- get("report", envir = parent.frame())
  }else{
    if(exists("report", envir = .GlobalEnv)){
      trend <- get("report", envir = .GlobalEnv)
    }else{

      trendtmp <- data.frame(Generation = numeric(),
                             Model = character(),
                             iter = numeric(),
                             P = numeric(),
                             varP = numeric(),
                             gv  = numeric(),
                             varG = numeric(),
                             S_gv = numeric(),
                             mv = numeric(),
                             varM = numeric(),
                             S_mv = numeric()
      )


      trend <- list(PopLow = list(No_symbiosis = trendtmp, Symbiosis = trendtmp),
                    PopHigh = list(No_symbiosis = trendtmp, Symbiosis = trendtmp))
    }
  }

  trendtmp <- data.frame(Generation = numeric(1),
                         Model = character(1),
                         iter = numeric(1),
                         P = numeric(1),
                         varP = numeric(1),
                         gv  = numeric(1),
                         varG = numeric(1),
                         S_gv = numeric(1),
                         mv = numeric(1),
                         varM = numeric(1),
                         S_mv = numeric(1)
  )



  trendtmp$Generation <- g - 1
  trendtmp$Model <- model
  trendtmp$iter <- iter

  #Phenotype
  trendtmp$P <- round(mean(pheno(pop),na.rm = T),2)
  trendtmp$varP <- round(var(pheno(pop), na.rm = T),2)[1]

  #Genotype
  trendtmp$gv <- round(meanG(pop),2)[1]
  trendtmp$varG <- round(varG(pop),2)[1]

  #To compute the differential of selection
  if(selType == "Low"){
    dam <- breedInd$Dam_low$Dam
    male <- breedInd$Male_low$Male
  }

  if(selType == "High"){
    dam <- breedInd$Dam_high$Dam
    male <- breedInd$Male_high$Male
  }

  trendtmp$S_gv <- round(mean(pop@gv[pop@id %in% c(dam, male)]),2)

  #Microbiome

  if(sym == 1){
    pop_mv <- unlist(lapply(pop@misc, `[[`, "mv_sym"))
    trendtmp$mv <- round(mean(pop_mv),2)
    trendtmp$varM <- round(var(pop_mv),2)

    #Differential of selection for the mv
    pop_mv <- NULL
    pop_mv <- unlist(lapply(pop@misc[pop@id %in% c(dam, male)], `[[`, "mv_sym"))
    trendtmp$S_mv <- round(mean(pop_mv), 2)

    if(selType == "Low"){
      trend$PopLow$Symbiosis <- rbind(trend$PopLow$Symbiosis, trendtmp)
    }

    if(selType == "High"){
      trend$PopHigh$Symbiosis <- rbind(trend$PopHigh$Symbiosis, trendtmp)
    }

  }else{
    pop_mv <- unlist(lapply(pop@misc, `[[`, "mv"))
    trendtmp$mv <- round(mean(pop_mv),2)
    trendtmp$varM <- round(var(pop_mv),2)

    #Differential of selection for the mv
    pop_mv <- NULL
    pop_mv <- unlist(lapply(pop@misc[pop@id %in% c(dam, male)], `[[`, "mv"))
    trendtmp$S_mv <- round(mean(pop_mv), 2)

    if(selType == "Low"){
      trend$PopLow$No_symbiosis <- rbind(trend$PopLow$No_symbiosis, trendtmp)
    }

    if(selType == "High"){
      trend$PopHigh$No_symbiosis <- rbind(trend$PopHigh$No_symbiosis, trendtmp)
    }

  }
  return(trend)
}
