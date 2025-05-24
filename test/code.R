##Setting simulation parameters
gSP <- GlobalSP$new(nPop = 100, nChr = 3, nQTLchr = 100, nyear = 100,nSire = 20, nDam = 500)
gSP$setTrait(meanP = 3.66, varP = 0.46, h2 = 0.16, m2 = 0.13)
gSP$setSpecies(nSpEff = 25, varMH = 0.028, meanMH = 0.47)

gSP$model <- "H"
gSP$dataM <- "data/example.txt"

##Founder population 
founderG = setFounderG(globalSP=gSP)
founderM = setFounderM(globalSP=gSP)

##Check simulation of the distribution of microbial heritability fits real data
h2_microbiota <- read.table("data/h2_microbioa.txt",header=T,sep=",")
#
df <- data.frame(
  value = c(h2_microbiota$h2asv, gSP$MH.H),
  type = factor(c(rep("Real", length(h2_microbiota$h2asv)), 
                  rep("Simulated", length(gSP$MH.H))))
)

# Plot
ggplot(df, aes(x = value, fill = type, color = type)) +
  geom_density(alpha = 0.3, size = 0.75) +
  labs(
    title = "Density Plot: Real vs Simulated Data",
    x = "Value",
    y = "Density",
    fill = "Data Type",
    color = "Data Type"
  ) +
  scale_color_manual(values = c("#8d1c1a","#00a7c7"))+
  scale_fill_manual(values = c("#8d1c1a","#00a7c7"))+
  scale_x_continuous(labels = function(x) round(x, 1)) +
  xlim(c(0,1))+

  theme_classic()
ggsave("test/microbial_heritability.tiff",width = 10, height = 10, units = "cm",dpi=600)

SP = essentialSP(founder = founderG, nSnpChr = 10000, minSnpFreq = 0.05)
founderPop = newPop(founderG)
founderPop@pheno

##Create Base population manually
gSP$s2 <- 0.5 #50% of VE. For variable h2, the s2 is a proportion of the VE to avoid VG + VS + VE > 1
basePop <- simBasePop(model = "H")

##Plotting simulated distribution of microbiota
df <- data.frame(
  value = c(data[,100],
            mbiome_sym[,100]),
  type = factor(c(rep("Real", nrow(data)),
                  rep("Simulated", nrow(mbiome)))
  )
)

##Plot
ggplot(df, aes(x = value, fill = type, color = type)) +
  geom_density(alpha = 0.3, size = 0.75) +
  labs(
    title = names(data)[100],
    x = "Value",
    y = "Density",
    fill = "Data Type",
    color = "Data Type"
  ) +
  scale_color_manual(values = c("#00a7c7","#8d1c1a"))+
  scale_fill_manual(values =c("#00a7c7","#8d1c1a"))+
  scale_x_continuous(labels = function(x) round(x, 1))+

  theme_classic()
ggsave("test/Prevotella_sym.tiff",width = 10, height = 10, units = "cm",dpi=600)


cor1 <- as.matrix(founderMxM)
colnames(gamma.sym) <- colnames(data)
cor2 <- cor(mbiome.sym)

corrplot(cor1, method = "color", order = "hclust", tl.pos="n")
corrplot(cor2, method = "color", order = "hclust", tl.pos="n")

# Assuming your matrices are named cor_mat1 and cor_mat2

# Save both plots to variables
corr_plot1 <- recordPlot({
  corrplot(cor1, method = "color", type = "lower",order="hclust", tl.cex = 0.6,
           title = "Original Correlation Matrix", mar = c(0,0,2,0))
})

corr_plot2 <- recordPlot({
  corrplot(cor_mat2, method = "color", type = "lower", order="hclust", tl.cex = 0.6,
           title = "Simulated Correlation Matrix", mar = c(0,0,2,0))
})

library(reshape2)
library(ggplot2)

# Step 1: Get clustering order (e.g., from original matrix)
hc <- hclust(as.dist(1 - cor1))
ord <- hc$labels[hc$order]

# Step 2: Reorder matrices
cor1_ord <- cor1[ord, ord]
cor2_ord <- cor2[ord, ord]

x <- cor1_ord[lower.tri(cor1_ord)]
y <- cor2_ord[lower.tri(cor2_ord)]

cor(x, y) 

# Step 3: Melt to long format
cor1_long <- melt(cor1_ord)
cor2_long <- melt(cor2_ord)

# Step 4: Create combined matrix showing:
# - cor1 in lower triangle
# - cor2 in upper triangle
combined <- cor1_long
combined$value <- ifelse(as.numeric(combined$Var1) > as.numeric(combined$Var2),
                         cor1_long$value,
                         ifelse(as.numeric(combined$Var1) < as.numeric(combined$Var2),
                                cor2_long$value, NA))

# Step 5: Plot
ggplot(combined, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_gradient2(low = "#00a7c7", mid = "white", high = "#8d1c1a",
                       midpoint = 0, limits = c(-1, 1), name = NULL) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(4, "cm"),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(fill = guide_colorbar(
    barwidth = 10, barheight = 0.5,
    ticks = FALSE, frame.colour = NA,
    title.position = "top"
  )) +
  labs(title = "Lower: Original | Upper: Simulated")
ggsave("test/cor_matrix.png",width = 15, height = 15, units = "cm",dpi=600)


pdf("All_species_density_plots_G+E+I.pdf", width = 8, height = 6)  

# Loop para cada especie (columna)
for (i in 1:ncol(data)) {
  
  df <- data.frame(
    value = c(data[, i], exp(mbiome_total[, i])),
    type = factor(c(rep("Real", nrow(data)), rep("Simulated", nrow(mbiome_total))))
  )
  
  p <- ggplot(df, aes(x = value, fill = type, color = type)) +
    geom_density(alpha = 0.3, size = 0.75) +
    labs(
      title = names(data)[i],
      x = "Value",
      y = "Density",
      fill = "Data Type",
      color = "Data Type"
    ) +
    scale_color_manual(values = c("#00a7c7","#8d1c1a")) +
    scale_fill_manual(values = c("#00a7c7","#8d1c1a")) +
    scale_x_continuous(labels = function(x) round(x, 1)) +
    theme_classic()
  
  print(p)  # Imprime el plot en el PDF
}

dev.off()  
