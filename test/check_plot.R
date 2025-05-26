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

##Checking population structure added
cor1 <- as.matrix(founderMxM)
colnames(mbiome_total) <- colnames(data)
cor2 <- cor(mbiome_total)

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


## Species distribution
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