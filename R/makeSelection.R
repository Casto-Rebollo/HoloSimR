#' Perform Selection
#' @description
#' This function automates the selection process across different generations and returns a report of class \code{retValues}.
#'
#' @param g Number of generations for selection. If \code{NULL}, it will be obtained from the global environment.
#' @param model A value for the model to simulate. Only values according to \dQuote{G}, \dQuote{M}, \dQuote{NMH}, \dQuote{LMH}, \dQuote{MMH}, \dQuote{HMH}, \dQuote{All} are available. If \code{NULL}, it will be obtained from the global environment.
#' @param selType Selection type (e.g., \dQuote{Divergent}, \dQuote{Low}, \dQuote{High}). If \code{NULL}, it will be obtained from the global environment.
#' @param sym A value of 1 to consider the symbiosis effect on microbiota. If \code{NULL}, it will be obtained from the global environment.
#' @param limTrait Vector with the minimum and maximum values for the trait. If \code{NULL}, it will be obtained from the global environment.
#' @param sex Indicates if the phenotype belongs exclusively to females ("F"). Default is \dQuote{both}.
#' @param LS \code{TRUE} if the phenotype being evaluated is litter size (LS). Default is \code{FALSE}.
#' @param globalSP An object of class \code{\link{GlobalSP}}. If \code{NULL}, it will be obtained from the global environment.
#' @param save If not null, the base and last population of the simulation will be saved.
#' @param iter Number of iterations to perform. If \code{NULL}, just one iteration will be done.
#' @param progressBar If TRUE, a progress bar and text will be displayed. \emph{Default} is \code{FALSE}.
#'
#' @return A list of class \code{retValues} with the following information:\enumerate{
#'   \item{\strong{Generation}}: Generation computed
#'   \item{\strong{Model}}: Model computed
#'   \item{\strong{iter}}: Iteration performed. In case the selection will be replicated
#'   \item{\strong{P}}: Average of the population phenotype
#'   \item{\strong{varP}}: Variance of the phenotype
#'   \item{\strong{gv}}: Average of the additive genetic value
#'   \item{\strong{varG}}: Genetic variance
#'   \item{\strong{S_gv}}: Selection differential for additive genetic values
#'   \item{\strong{mv}}: Average of the microbiome values
#'   \item{\strong{varM}}: Variance of the microbiome values
#'   \item{\strong{S_mv}}: Selection differential for the microbiome values
#'   }
#'
#' @export
#' @import AlphaSimR
#' @import dplyr
#' @importFrom utils txtProgressBar setTxtProgressBar write.table
#'
#' @usage
#' makeSelection(g = NULL,
#'               model = NULL,
#'               selType = NULL,
#'               sym = NULL,
#'               limTrait = NULL,
#'               sex = "both",
#'               LS = FALSE,
#'               globalSP = NULL,
#'               save = c(),
#'               iter = 1,
#'               progressBar = FALSE)
#'
#' @examples
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
#' # Initialize the selection process
#' selection <- makeSelection(sex = "F", LS = TRUE, limTrait = c(0,18))

makeSelection <- function(g = NULL, model = NULL, selType = NULL, sym = NULL,
                          limTrait = NULL, sex = "both", LS = FALSE, globalSP = NULL,
                          save = c(), iter = 1, progressBar = FALSE){

  #Vector in model is supported
  if(is.null(globalSP)){
    globalSP = get("gSP", envir = .GlobalEnv)
  }

  if(is.null(model)){
    model = globalSP$model
  }

  if(is.null(g)){
    g = globalSP$nyear
  }

  if(is.null(limTrait)){
    limTrait = globalSP$limit.Trait
  }

  if(is.null(sym)){
    sym = globalSP$symbiosis
  }

  if(is.null(selType)){
    selType = globalSP$selType
  }

  if(is.null(save)){
    save <- c(0,g)
  }

  selOrder <- switch(
    selType,
    Divergent = c(FALSE,TRUE),
    Low = FALSE,
    High = TRUE
  )

  nson <- globalSP$nSon

  if(length(model) == 1){
    if(model == "All" && is.null(globalSP$MH.H)){
      model = c("G",
                "M",
                "NMH",
                "LMH",
                "MMH",
                "HMH")
    }
    if(model == "All" && !is.null(globalSP$MH.H)){
      model = c("G",
                "M",
                "NMH",
                "H")
    }
  }


  # Initialize time tracking
  start_time <- Sys.time()
  for(iter_n in 1:iter){
    rndSeed_set <- sample(c(1000:2000),1)

    count = 0
    for (iter_model in model) {

      if(progressBar == TRUE){
        #To optimize text
        n <- (20 - nchar(iter_model))%/%2
        if((20 - nchar(iter_model))/2!=0){
          n_init <- n
        }else{
          n_init <- n-2
        }


        writeLines(paste("##########################\n","#",
                         paste(rep(" ",n_init),collapse = ""),"Model ",iter_model,paste(rep(" ",n-1),collapse=""),"#\n",
                         "##########################",sep=""))

        writeLines(">>>> Base population <<<<<")
      }


      #Starting the model, use the first seed
      #Create the base population
      basePop <- simBasePop(model = iter_model, rndSeed = rndSeed_set[1],
                            progressBar = progressBar) #Vector for model is not supported

      for(value in selOrder){
        selType_iter <- if(!value) "Low" else "High"

        if(progressBar == TRUE){
          writeLines(">>> Starting Selection <<<")
          writeLines(paste(">>>> ",toupper(selType_iter)
                           , " population <<<<<<", sep=""))
        }


        for(iter_sym in sym){

          if(progressBar == TRUE){
            writeLines(paste("**** Symbiosis  ", iter_sym, "****"))
          }


          #Count the number of model executed (needed for seed)
          count <- count + 1


          #Compute phenotype
          Pop <- makeP(pop = basePop$Pop, model = iter_model,
                       sym = iter_sym, limTrait = limTrait,
                       sex = sex, rndSeed = rndSeed_set[1])
          parentPop <- Pop

          #Select breeding animals
          parent <- selectBreeding(pop = parentPop, sex = sex, g0 = TRUE, LS = LS)

          #Creation of mating plan for next generation
          crossPlan <- matingPLAN(parent = parent)

          report <- retValues(model = iter_model, g = 1,
                              sym = iter_sym, selType = selType_iter,
                              pop = parentPop,
                              breedInd = parent, iter = iter_n)

          if((any(save == "All") | any(save == 0)) & iter_n == iter){
            sortedPop <- savePop(Pop, model = iter_model, g = 0, save = sortedPop)
          }

          if(progressBar == TRUE){
            pb <- txtProgressBar(min = 1, max = g, style = 3)
          }


          for (iter_g in 1:g) {

            if(progressBar == TRUE){
              setTxtProgressBar(pb, iter_g)
            }

            if(iter_g == 1){
              Pop <- nextPop(pop = parentPop, crossPlan = crossPlan, nSon = nson, LS = LS)
            }else{
              #Create a new Population
              Pop <- nextPop(pop = Pop, crossPlan = crossPlan,
                             nSon = nson, selType = selType_iter, LS = LS)
            }


            #Save seed to sample the same values to homogenieze effect throught models
            if(count == 1){
              rndSeed_set[iter_g + 1] <- sample(c(1000:2000),1)
            }

            #New population based of breeding animals selected.

            Poptmp <- if (inherits(Pop,"list")) {
              if (!value) Pop$PopLow else Pop$PopHigh
            } else {
              Pop
            }

            acquiredSp <- acquiredSpecies(pop = Poptmp, sym = iter_sym,
                                          parentPop = parentPop,
                                          rndSeed = rndSeed_set[iter_g + 1])

            #Creating microbiota matrix
            Poptmp <- makeM(pop = Poptmp,
                            sym = iter_sym , rndSeed = rndSeed_set[iter_g + 1])

            #Simulation phenotype
            Pop <- makeP(pop = Poptmp,
                         model = iter_model,
                         sym = iter_sym, limTrait = limTrait ,
                         sex = sex, rndSeed = rndSeed_set[iter_g + 1])


            #Selection breeding animals and creation of mating plan
            parent <- selectBreeding(pop = Pop,
                                     sex = sex,
                                     selType = selType_iter, LS = LS)

            crossPlan <- matingPLAN(parent = parent, selType = selType_iter)

            #Save values of population generated
            report <- retValues(model = iter_model, g = iter_g + 1,
                                sym = iter_sym, selType = selType_iter,
                                pop = Pop,
                                breedInd = parent, iter = iter_n)
            #Save parental population
            parentPop <- Pop

          }
          if(progressBar == TRUE){
            close(pb)
          }


          ## Save populations
          if((any(save == iter_g) | any(save == "All")) & iter_n == iter){
            sortedPop <- savePop(Pop, model = iter_model,g = iter_g,save = sortedPop)
          }
        }
      }
    }
  }

  end_time <- Sys.time()

  execution_time <- difftime(end_time, start_time, units = "mins")
  print(paste("Execution time:", round(as.numeric(execution_time), 2), "minutes"))

  if(inherits(save, "numeric") | any(save == "All")){
    return(list("Report" = report,"Pop" = sortedPop))
  }else{
    return(report)
  }


}
