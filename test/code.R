##Setting simulation parameters
gSP <- GlobalSP$new(nPop = 100, nChr = 3, nQTLchr = 100, nyear = 100,nSire = 20, nDam = 500)
gSP$setTrait(meanP = 3.66, varP = 0.46, h2 = 0.16, m2 = 0.13)
gSP$setSpecies(nSpEff = 25, varMH = 0.016, meanMH = 0.45)

gSP$model <- "H"
gSP$dataM <- "data/example.txt"

##Founder population 
founderG = setFounderG(globalSP=gSP)
founderM = setFounderM(globalSP=gSP)

##Check simulation of the distribution of microbial heritability fits real data
#h2_microbiota <- read.table("data/h2_microbioa.txt",header=T,sep=",")
#
#df <- data.frame(
#  value = c(h2_microbiota$h2asv, gSP$MH.H),
#  type = factor(c(rep("Real", length(h2_microbiota$h2asv)), 
#                  rep("Simulated", length(gSP$MH.H))))
#)

# Plot
#ggplot(df, aes(x = value, fill = type, color = type)) +
#  geom_density(alpha = 0.3, size = 0.75) +
#  labs(
#    title = "Density Plot: Real vs Simulated Data",
#    x = "Value",
#    y = "Density",
#    fill = "Data Type",
#    color = "Data Type"
#  ) +
#  scale_color_manual(values = c("#8d1c1a","#4c83c8"))+
#  scale_fill_manual(values = c("#8d1c1a","#4c83c8"))+
#  scale_x_continuous(labels = function(x) round(x, 1)) +
#  xlim(c(0,1))+

#  theme_classic()
#ggsave("test/microbial_heritability.tiff",width = 10, height = 10, units = "cm",dpi=600)

SP = essentialSP(founder = founderG, nSnpChr = 10000, minSnpFreq = 0.05)
founderPop = newPop(founderG)
founderPop@pheno

##Create Base population manually
gSP$s2 <- 0.5 #50% of VE. For variable h2, the s2 is a proportion of the VE to avoid VG + VS + VE > 1
basePop <- simBasePop(model = "H")


