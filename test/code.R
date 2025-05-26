##Setting simulation parameters
gSP <- GlobalSP$new(nPop = 100, nChr = 3, nQTLchr = 100, 
                    nyear = 10,nSire = 20, nDam = 500,
                    selType = "Low")

gSP$setTrait(meanP = 3.66, varP = 0.46, h2 = 0.16, m2 = 0.13)
gSP$setSpecies(nSpEff = 25, varMH = 0.028, meanMH = 0.47)

gSP$model <- "H"
gSP$dataM <- "data/example.txt"

##Founder population 
founderG = setFounderG(globalSP=gSP)
founderM = setFounderM(globalSP=gSP)

SP = essentialSP(founder = founderG, nSnpChr = 10000, minSnpFreq = 0.05)
founderPop = newPop(founderG)
founderPop@pheno

##Create Base population manually
gSP$s2 <- 0.5 #50% of VE. For real values of h2 (non-constant), the s2 is a proportion of the VE to avoid VG + VS + VE > 1
basePop <- simBasePop(model = "H")

# Compute phenotype in a base population
Pop <- makeP(pop = basePop$Pop, model = "H", sym = 1,
            limitTrait = c(0, NA))
hist(Pop@pheno); range(Pop@pheno)
