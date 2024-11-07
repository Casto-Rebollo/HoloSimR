#' Selection of Breeding Animals
#' @description
#' Selection of breeding animals according to the Phenotype (\dQuote{phe}) or the direct genetic value (\dQuote{gv}).
#'
#' @param pop an object of \code{\link{MapPop-class}}
#' @param nDam number of females to select. If \code{NULL}, it will be set by \code{nDam} from \code{\link{GlobalSP}}. Default 125.
#' @param nSire number of males to select. If \code{NULL}, it will be set by \code{nSire} from \code{\link{GlobalSP}}. Default 25.
#' @param selType Selection type (e.g., \dQuote{Divergent}, \dQuote{Low}, \dQuote{High}). If \code{NULL}, it will be obtained from the global environment
#' @param selby type of selection; Phenotypic selection (\dQuote{phe}; Default) or Genomic selection (\dQuote{gv})
#' @param sex indicate if the phenotype belongs exclusively to females ("F"). Default \dQuote{both}
#' @param maxFS maximum number of full siblings allowed when selecting females.
#' @param g0 \code{TRUE} for indicating that the animals will be selected from the base population. Default \code{FALSE}
#' @param LS \code{TRUE} if the phenotype being evaluated is litter size (LS). Default \code{FALSE}
#' @param globalSP an object of class \code{\link{GlobalSP}}. If \code{NULL}, it will be obtained from the global environment.
#'
#' @return a list including four dataframes:\enumerate{
#'   \item{\strong{Dam_low}}: ID of breeding females for decreasing the phenotype, mother ID, father ID and value for the genetic value or phenotype, depending on the selby
#'   \item{\strong{Dam_high}}: ID of breeding females for increasing the phenotype, mother ID, father ID and value for the genetic value or phenotype, depending on the selby
#'   \item{\strong{Male_low}}: ID of breeding males for decreasing the phenotype, mother ID, and father ID
#'   \item{\strong{Male_high}}: ID of breeding males for increasing the phenotype, mother ID, and father ID
#'   }
#'
#' @export
#'
#' @usage
#' selectBreeding(pop,
#'                nDam = NULL,
#'                nSire = NULL,
#'                selType = NULL,
#'                selby = "phe",
#'                sex = "both",
#'                maxFS = 2,
#'                g0 = FALSE,
#'                LS = FALSE,
#'                globalSP = NULL)
#'
#' @examples
#' # Set simulation global parameters
#' gSP = GlobalSP$new(nPop = 1000, nyear = 10, nQTLchr = 100)
#'
#' # Defining the microbiota
#' gSP$setSpecies(nSpecies = 1000, nSp0 = 600, nSpEff = 100, symbiosis = c(0,1))
#'
#' # Defining trait
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
#' # Computing phenotype in a base population
#' Pop <- makeP(pop = basePop$Pop, model = "NMH", sym = 1, limTrait = c(0,18), sex = "F")
#'
#' # Selection of breeding animals from basePop
#' parent <- selectBreeding(pop = Pop, nDam = 125, nSire = 25, sex = "F", g0 = TRUE)
#'
selectBreeding <- function(pop, nDam = NULL, nSire = NULL,
                           selType = NULL,selby = "phe", sex = "both",
                           maxFS = 2, g0 = FALSE, LS = FALSE, globalSP = NULL){
  if(is.null(globalSP)){
    globalSP = get("gSP", envir = .GlobalEnv)
  }
  if(is.null(nDam)){
    nDam = globalSP$nDam
  }
  if(is.null(nSire)){
    nSire = globalSP$nSire
  }

  if(is.null(selType)){
    selType <- globalSP$selType
  }


  ##Select breeding females

  selOrder <- switch(
    selType,
    Divergent = c(FALSE,TRUE),
    Low = FALSE,
    High = TRUE
  )

  parents <- list(Dam_low = NULL, Dam_high = NULL,
                  Male_low = NULL, Male_high = NULL)

  #Loop for select each breeding animals depending on the selection process
    for(value in selOrder){
        Poptmp <- pop[pop@sex=="F"]
        select.Dam <- data.frame(Dam = character(nDam),
                                 father = character(nDam),
                                 mother = character(nDam),
                                 Sel_value = numeric(nDam))

        #Order by gv or phenotype & considering the biological limit of the trait
        if(selby == "gv"){
          Poptmp <- Poptmp[order(Poptmp@gv, decreasing = value)]
          select.Dam$Sel_value <- as.numeric(Poptmp@gv[1:nDam])
        }else{
          Poptmp <- Poptmp[order(Poptmp@pheno, decreasing = value)]
          select.Dam$Sel_value <- as.numeric(Poptmp@pheno[1:nDam])
        }

        #Remove from the population to select the dam with no offspring

        if(LS == TRUE){
          Poptmp@pheno <- round(Poptmp@pheno)
          Poptmp <- Poptmp[Poptmp@pheno!=0,]
        }


      select.Dam$Dam <- Poptmp@id[1:nDam]
      select.Dam$mother <- Poptmp@mother[1:nDam]
      select.Dam$father <- Poptmp@father[1:nDam]


      #Select just a limit of full siblings
      nLitter <- data.frame(table(paste(select.Dam$mother,select.Dam$father,sep="")))
      nLitter$Var1 <- as.character(nLitter$Var1)
      nLitter$FS <- maxFS - nLitter$Freq

      extraFS <- nLitter$Var1[nLitter$FS<0]
      rm.index <- NULL
      if(length(extraFS) !=0){
        for(extra in extraFS){
          lit <- paste(select.Dam$mother,select.Dam$father,sep="")
          indx <- which(lit == extra)
          rm.index <- c(rm.index,indx[(maxFS+1):length(indx)])
        }

        select.Dam <- select.Dam[-rm.index,]
        nLitter$FS[nLitter$FS<0] <- 0

        ###Select rest of females we need
        nD <- nDam + 1
        while(length(select.Dam$Dam)!=nDam | nD == length(Poptmp@id)){

          lit <- paste(Poptmp@mother[nD],Poptmp@father[nD],sep="")
          cond <- FALSE
          if(length(nLitter$FS[nLitter$Var1==lit]) == 0){
            ##Split on two to update nLitter dataframe
            nLitter <- rbind(nLitter,c(lit,1,maxFS-1))
            cond <- TRUE

          }else if(nLitter$FS[nLitter$Var1==lit]>0){
            nLitter$FS[nLitter$Var1==lit] <- nLitter$FS[nLitter$Var1==lit]-1
            cond <- TRUE
          }

          if(cond == TRUE & selby == "gv"){
            select.Dam <- rbind(select.Dam,
                                c(Poptmp@id[nD],
                                  Poptmp@father[nD],
                                  Poptmp@mother[nD],
                                  Poptmp@gv[nD]))
          }
          if(cond == TRUE & selby == "phe"){
            select.Dam <- rbind(select.Dam,
                                c(Poptmp@id[nD],
                                  Poptmp@father[nD],
                                  Poptmp@mother[nD],
                                  Poptmp@pheno[nD]))
          }

          nD <- nD + 1
          nLitter$FS <- as.numeric(nLitter$FS)
        }

      }


      ## Selection breeding males
      select.m <- NULL
      tmp.male <- pop[pop@sex=="M",]

      select.male <- data.frame(Male = character(nSire),
                 father = character(nSire),
                 mother = character(nSire))

      #Selection of males if the phenotype come from just from the females
      if(sex == "F" & selby == "phe"){
        if(g0 == TRUE){
          lit.dam <- unique(paste(select.Dam$mother, select.Dam$father, sep=""))[1:nSire]
          lit <- paste(tmp.male@mother,tmp.male@father,sep="")
          for (lit.d in lit.dam) {

            select.m <- c(select.m,sample(tmp.male@id[which(lit==lit.d)],1))

          }
        }else{

          if (!exists("crossPlan", envir = parent.frame())) {
            if(!exists("crossPlan", envir = .GlobalEnv)){
              stop("Error in get(\"crossPlan\", envir = .GlobalEnv) : objeto 'crossPlan' no encontrado. Are you in the generation 0?")
            }
          }

          if(exists("crossPlan", envir = parent.frame())){
            crossplan = get("crossPlan", envir = parent.frame())
            old_parents = get("parent", envir = parent.frame())
          }else{
            crossplan = get("crossPlan", envir = .GlobalEnv)
            old_parents = get("parent", envir = .GlobalEnv)
          }


          ###Include information
          if(value == FALSE){
            male <- old_parents$Male_low$Male
            crossplan <- crossplan$PopLow
          }else{
            male <- old_parents$Male_high$Male
            crossplan <- crossplan$PopHigh
          }


          #Male from the best cross with the female

          for(sire in male){

            dam.sel <- "X"

            subset <- crossplan[crossplan$Sire == sire,]

            dam <- subset$Dam[order(subset$Sel_value, decreasing = value)]
            dam <- dam[dam!=dam.sel]

            n <- 1
            offspring <- pop@id[(pop@mother == dam[n] &
                                     pop@father == sire & pop@sex == "M")]



            #Select a female with males in their offspring
            while(length(offspring)==0  & n!=5){
              n = n + 1
              offspring <- pop@id[pop@mother == dam[n] &
                                       pop@father == sire &
                                       pop@sex == "M"]

            }
            select.m <- c(select.m,sample(offspring, 1))
            dam.sel <- dam[n]

          }

        }

        select.male$Male <- select.m
        select.male$mother <- pop@mother[match(select.m,pop@id)]
        select.male$father <- pop@father[match(select.m,pop@id)]

      }else{

        if(selby == "gv"){
          Poptmp <- tmp.male[order(tmp.male@gv, decreasing = value)]
          select.male$Sel_value <- as.numeric(Poptmp@gv[1:nSire])

        }else{
          Poptmp <- tmp.male[order(tmp.male@pheno, decreasing = value)]
          select.male$Sel_value <- as.numeric(Poptmp@pheno[1:nSire])
        }


        select.male$Male <- tmp.male@id[1:nSire]
        select.male$mother <- tmp.male@mother[1:nSire]
        select.male$father <- tmp.male@father[1:nSire]

        #Select just a limit of full siblings
        nLitter <- data.frame(table(paste(select.male$mother,select.male$father,sep="")))
        nLitter$Var1 <- as.character(nLitter$Var1)
        nLitter$FS <- maxFS - nLitter$Freq

        extraFS <- nLitter$Var1[nLitter$FS<0]
        rm.index <- NULL
        if(length(extraFS) !=0){
          for(extra in extraFS){
            lit <- paste(select.male$mother,select.male$father,sep="")
            indx <- which(lit == extra)
            rm.index <- c(rm.index,indx[(maxFS+1):length(indx)])
          }

          select.male <- select.male[-rm.index,]
          nLitter$FS[nLitter$FS<0] <- 0

          ###Select rest of females we need
          nD <- nSire + 1
          while(length(select.male$Male)!=nSire | nD == length(Poptmp@id)){

            lit <- paste(Poptmp@mother[nD],Poptmp@father[nD],sep="")
            cond <- FALSE
            if(length(nLitter$FS[nLitter$Var1==lit]) == 0){
              ##Split on two to update nLitter dataframe
              nLitter <- rbind(nLitter,c(lit,1,maxFS-1))
              cond <- TRUE

            }else if(nLitter$FS[nLitter$Var1==lit]>0){
              nLitter$FS[nLitter$Var1==lit] <- nLitter$FS[nLitter$Var1==lit]-1
              cond <- TRUE
            }

            if(cond == TRUE & selby == "gv"){
              select.male <- rbind(select.male,
                                   c(Poptmp@id[nD],
                                     Poptmp@father[nD],
                                     Poptmp@mother[nD],
                                     Poptmp@gv[nD]))
            }
            if(cond == TRUE & selby == "phe"){
              select.male <- rbind(select.male,
                                   c(Poptmp@id[nD],
                                     Poptmp@father[nD],
                                     Poptmp@mother[nD],
                                     Poptmp@pheno[nD]))
            }

            nD <- nD + 1
            nLitter$FS <- as.numeric(nLitter$FS)
          }

        }

      }

      if(value == FALSE){
        parents$Dam_low <- select.Dam
        parents$Male_low <- select.male
        if(sex!="F"){
          parents$Male_low$Sel_value <- as.numeric(parents$Male_low$Sel_value)
        }
        parents$Dam_low$Sel_value <- as.numeric(parents$Dam_low$Sel_value)
      }else{
        parents$Dam_high <- select.Dam
        parents$Male_high <- select.male
        if(sex!="F"){
          parents$Male_high$Sel_value <- as.numeric(parents$Male_high$Sel_value)
        }
        parents$Dam_high$Sel_value <- as.numeric(parents$Dam_high$Sel_value)
      }
    }
  return(parents)
}
