#-------------------------------------------------------------------------
## Mandatory execution
#-------------------------------------------------------------------------

##Setting simulation parameters
gSP <- GlobalSP$new(nPop = 100, nChr = 3, nQTLchr = 100, 
                    nyear = 10,nSire = 20, nDam = 500)

gSP$setTrait(meanP = 3.66, varP = 0.46, h2 = 0.16, m2 = 0.13)
gSP$setSpecies(nSpEff = 25, varMH = 0.028, meanMH = 0.47)

gSP$model <- "H"
gSP$dataM <- "data/example.txt"

##Founder population 
founderG = setFounderG(globalSP=gSP)
founderM = setFounderM(globalSP=gSP)

SP = essentialSP(founder = founderG, nSnpChr = 10000, minSnpFreq = 0.05)
founderPop = newPop(founderG)
#founderPop@pheno

gSP$s2 <- 0.5 #50% of VE. For real values of h2 (non-constant), the s2 is a proportion of the VE to avoid VG + VS + VE > 1


#-------------------------------------------------------------------------
## Create Base population manually
#-------------------------------------------------------------------------
basePop <- simBasePop(model = "H")

# Compute phenotype in a base population
Pop <- makeP(pop = basePop$Pop, model = "H", sym = 1,
             limTrait = c(0,NA))
#hist(Pop@pheno); range(Pop@pheno)

# Select breeding animals from basePop
parent <- selectBreeding(pop = Pop, g0 = TRUE, sym = 1, selby="gv_mv")
#parent$Male_low ; parent$Dam_low

# Create mating plan for next generation
crossPlan <- matingPLAN(parent = parent)

# Generate first generation of selection
Pop <- nextPop(pop = Pop, crossPlan = crossPlan)
#hist(Pop$PopLow@gv)

acquiredSp <- acquiredSpecies(pop = Pop$PopLow, sym = 1)

#Creating microbiota matrix
Poptmp <- makeM(pop = Pop$PopLow, sym = 1)

#Simulation phenotype
Pop <- makeP(pop = Poptmp, model = "H", sym = 1)

#-------------------------------------------------------------------------
## Automatize all selection (Paper)
#-------------------------------------------------------------------------
# rm(list = ls())
# After founderPop creation
selection <- makeSelection(iter = 10, model = "H")
