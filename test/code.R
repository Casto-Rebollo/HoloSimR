##Setting simulation parameters
gSP <- GlobalSP$new(nPop = 100, nChr = 3, nQTLchr = 100, nyear = 100,nSire = 20, nDam = 500)
gSP$setTrait(meanP = 3.66, varP = 0.46, h2 = 0.16, m2 = 0.13)
gSP$setSpecies(nSpEff = 25, varMH = 0.03, meanMH = 0.47)

gSP$model <- "H"
gSP$dataM <- "data/example.txt"

##Founder population 
founderG = setFounderG(globalSP=gSP)
founderM = setFounderM(globalSP=gSP)
