#' Plot Trend
#'
#' @description
#' \code{plotREND} generates a linear plot of genetic, phenotypic, and microbiome trends across generations for one scenario or for all, focusing on a single effect.
#'
#' @param data A data.frame from the \code{retValues} object.
#' @param plot The type of plot to generate (\dQuote{boxplot} or \dQuote{line}). The default is \dQuote{line}, which creates trends for all effects within a single model when \code{draw} is set to \dQuote{Model} (the \code{model} parameter is required). To plot trends of all models for a specific effect, use the \dQuote{Effect} option in \code{draw}. The default is \dQuote{P}, which plots the phenotype for all evaluated models. To change it, use the \code{effect} argument with \dQuote{gv} or \dQuote{mv}. For boxplot, just data is required.
#' @param draw To graph all the effect for a specific model (\dQuote{Model}; requires \code{model} argument) or all models for a specific effect (\dQuote{Effect}; requires \code{effect} argument). Default \dQuote{Model}
#' @param model The model for which to generate the plot line for \code{draw} option equal to \dQuote{Model}.
#' @param effect The effect to plot when \code{plot} is set to \dQuote{Effect}. The default is \dQuote{P} for Phenotype. Other options are \dQuote{gv} for genetic values and \dQuote{mv} for microbiota values.
#' @param type_se Graph style for plotting the standard error (\dQuote{ribbon} or \dQuote{classic}). Default is \dQuote{classic}.
#' @param data2 A data.frame from the \code{retValues} object. Required if \code{selType} is set to \dQuote{Divergent}.
#' @param selType The population to plot. The default is \dQuote{High}.
#' @param name The name of the graph file.
#' @param path The path where the graph file will be saved.
#' @param show If TRUE, the plot will be displayed in the current R session.
#' @param dpi The resolution of the image. The default is 600 dpi.
#' @param width The width of the image. The default is 17 cm.
#' @param height The height of the image. The default is 12 cm.
#' @param units The units for the width and height.
#'
#' @return A graph file and a \code{\link{ggplot2}} object.
#' @import ggplot2
#' @importFrom stats sd
#' @export
#'
#' @usage plotREND(data,
#'                 plot = "line",
#'                 draw = "Model",
#'                 model = NULL,
#'                 effect = "P",
#'                 type_se = "classic",
#'                 data2 = NULL,
#'                 selType = "High",
#'                 name = NULL,
#'                 path = "./",
#'                 show = FALSE,
#'                 dpi = 600,
#'                 width = 17,
#'                 height = 12,
#'                 units = "cm")
#' @examples
#' \dontrun{
#' plotREND(data=selection[["Report"]][["PopHigh"]]$Symbiosis,
#'          plot = "boxplot",
#'          selType = "High",
#'          show = TRUE,
#'          name = "lastG_boxplot_symbiosis_high.tiff")
#' }
plotREND <- function(data,
                     plot = "line",
                     draw = "Model",
                     model = NULL,
                     effect = "P",
                     type_se = "classic",
                     data2 = NULL,
                     selType = "High",
                     name = NULL, path = "./",
                     show = FALSE,
                     dpi = 600, width = 17, height = 12, units = "cm"){

  if(selType == "Divergent" & is.null(data2)){
    stop("Two datasets are required when selType is set to 'Divergent'")
  }

  iter <- max(data$iter)
  nyear <- max(data$Generation)

  if(plot == "boxplot"){
    data.tmp <- data
    for(m in unique(data.tmp$Model)){
      for(i in 1:iter){
        data.tmp$P[data.tmp$Model == m & data.tmp$iter == i] <- data.tmp$P[data.tmp$Model == m & data.tmp$iter == i] -
          data.tmp$P[data.tmp$Model == m & data.tmp$iter == i][1]
      }
    }
    data.tmp <- data.tmp[data.tmp$Generation == nyear,]


    color_scale <- c("G" = "#00a7c7", "M" = "#8d1c1a", "NMH" = "#708090",
                     "LMH" = "#4c83c8", "MMH" = "#835faa", "HMH" = "#9e3773")
    #effect_order <- c("G", "M", "NMH", "LMH", "MMH", "HMH")

    color_effect <- c("Phenotypic" = "#AD70AD", "Genetic" = "#00a7c7", "Microbiota" = "#8d1c1a")
    effect_order <- c("Phenotypic", "Genetic", "Microbiota")

    data.plot <- data.frame(Value= c(data.tmp$P,data.tmp$gv,data.tmp$mv))
    data.plot$Effect <- c(rep("Phenotypic", length(data.tmp$Generation)),rep("Genetic",length(data.tmp$Generation)),rep("Microbiota",length(data.tmp$Generation)))
    data.plot$Scenario <- rep(data.tmp$Model,3)

    y.min <- min(data.plot$Value) + min(data.plot$Value)/10
    y.max <- max(data.plot$Value) + max(data.plot$Value)/10

    # Ensure Effect is a factor with the correct order
    data.plot$Effect <- factor(data.plot$Effect, levels = c("Phenotypic","Genetic", "Microbiota"))

    # Create the plot
    p <- ggplot(data.plot, aes(x = Scenario, y = Value)) +
      # Add background shading
      annotate(geom = "rect",
               xmin = c(0.55, 1.55, 2.55, 3.55, 4.55, 5.55),
               xmax = c(1.45, 2.45, 3.45, 4.45, 5.45, 6.45),
               ymin = -Inf, ymax = Inf, alpha = 0.1,
               fill = color_scale) +
      # Add boxplots
      geom_boxplot(aes(fill = Effect), alpha = 1, outlier.size = 0.5) +
      # Set scales
      scale_fill_manual(values = color_effect) +
      scale_x_discrete(limits = c("G", "M", "NMH", "LMH", "MMH", "HMH")) +
      # Set labels
      labs(
        x = "Scenario",
        y = "Value"
      ) +
      # Set theme and other aesthetics
      theme_classic() +
      theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.25, 'cm'),
        legend.position = "right",
        panel.grid.major.y = element_line(colour = "grey", size = 0.2),
        panel.grid.minor.y = element_line(colour = "grey", size = 0.2),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
      ) +
      # Set y-axis limits
      ylim(y.min, y.max)

  }else{

    #Number of data to plotting
    if(draw == "Model"){
      if(is.null(model)){
        stop("A model is required")
      }

      n <- 3*(nyear+1)

    }else{
      n <- length(unique(data$Model))*(nyear+1)
    }

    #Set databases
    data.tmp <- data.frame(Generation = numeric(n),
                           Population = rep("data",n),
                           Value = numeric(n),
                           CI_min = numeric(n),
                           CI_max = numeric(n),
                           Effect = character(n))

    if(iter > 1){
      tmp <- summarySel(data,c(0:nyear))
      if(draw == "Model"){
        tmp <- tmp[tmp$Model == model,]
      }

    }else{
      if(draw == "Model"){
        tmp <- data[data$Model == model,]
      }
      names(tmp)[which(names(tmp)=="P")] <- "P_mean"
      names(tmp)[which(names(tmp)=="gv")] <- "gv_mean"
      names(tmp)[which(names(tmp)=="mv")] <- "mv_mean"
      tmp$gv_sd <- 0
      tmp$mv_sd <- 0
      tmp$P_sd <- 0


    }

    if(selType == "Divergent"){
      data.tmp2 <- data.frame(Generation = numeric(n),
                              Population = rep("data2",n),
                              Value = numeric(n),
                              CI_min = numeric(n),
                              CI_max = numeric(n),
                              Effect = character(n))
      if(iter > 1){
        tmp2 <- summarySel(data2,c(0:nyear))
        if(draw == "Model"){
          tmp2 <- tmp2[tmp2$Model == model,]
        }
      }else{
        if(draw == "Model"){
          tmp2 <- data2[data2$Model == model,]
        }
        names(tmp2)[which(names(tmp2)=="P")] <- "P_mean"
        names(tmp2)[which(names(tmp2)=="gv")] <- "gv_mean"
        names(tmp2)[which(names(tmp2)=="mv")] <- "mv_mean"
        tmp2$gv_sd <- 0
        tmp2$mv_sd <- 0
        tmp2$P_sd <- 0
      }
    }

    #Generate database for plotting depend of the type of plot
    if(draw == "Model"){
      ng <- nyear + 1
      data.tmp$Effect <-   c(rep("P",nyear+1),rep("gv",nyear+1),rep("mv", nyear+1))
      for(i in 1:ng){
        tmp$P_mean <- tmp$P_mean - tmp$P_mean[1]
        data.tmp$Generation[i] <- i-1
        data.tmp$Generation[ng+i] <- i-1
        data.tmp$Generation[2*ng+i] <- i-1
        data.tmp$Value[i] <- tmp$P_mean[i]
        data.tmp$CI_min[i] <- tmp$P_mean[i] - 1.96*(tmp$P_sd[i]/sqrt(iter))
        data.tmp$CI_max[i] <- tmp$P_mean[i] + 1.96*(tmp$P_sd[i]/sqrt(iter))
        data.tmp$Value[ng+i] <- tmp$gv_mean[i]
        data.tmp$CI_min[ng+i] <- tmp$gv_mean[i] - 1.96*(tmp$gv_sd[i]/sqrt(iter))
        data.tmp$CI_max[ng+i] <- tmp$gv_mean[i]+ 1.96*(tmp$gv_sd[i]/sqrt(iter))
        data.tmp$Value[2*ng+i] <- tmp$mv_mean[i]
        data.tmp$CI_min[2*ng+i] <- tmp$mv_mean[i] - 1.96*(tmp$mv_sd[i]/sqrt(iter))
        data.tmp$CI_max[2*ng+i] <- tmp$mv_mean[i] + 1.96*(tmp$mv_sd[i]/sqrt(iter))

        if(selType == "Divergent"){
          ## Other dataset
          data.tmp2$Effect <-   c(rep("P",nyear+1),rep("gv",nyear+1),rep("mv", nyear+1))
          tmp2$P_mean <- tmp2$P_mean - tmp2$P_mean[1]
          data.tmp2$Generation[i] <- i-1
          data.tmp2$Generation[ng+i] <- i-1
          data.tmp2$Generation[2*ng+i] <- i-1
          data.tmp2$Value[i] <- tmp2$P_mean[i]
          data.tmp2$CI_min[i] <- tmp2$P_mean[i] - 1.96*(tmp2$P_sd[i]/sqrt(iter))
          data.tmp2$CI_max[i] <- tmp2$P_mean[i] + 1.96*(tmp2$P_sd[i]/sqrt(iter))
          data.tmp2$Value[ng+i] <- tmp2$gv_mean[i]
          data.tmp2$CI_min[ng+i] <- tmp2$gv_mean[i] - 1.96*(tmp2$mv_sd[i]/sqrt(iter))
          data.tmp2$CI_max[ng+i] <- tmp2$gv_mean[i]+ 1.96*(tmp2$mv_sd[i]/sqrt(iter))
          data.tmp2$Value[2*ng+i] <- tmp2$mv_mean[i]
          data.tmp2$CI_min[2*ng+i] <- tmp2$mv_mean[i] - 1.96*(tmp2$mv_sd[i]/sqrt(iter))
          data.tmp2$CI_max[2*ng+i] <- tmp2$mv_mean[i] + 1.96*(tmp2$mv_sd[i]/sqrt(iter))
        }
      }
    }else{
      for(m in unique(tmp$Model)){
        tmp$P_mean[tmp$Model == m] <- tmp$P_mean[tmp$Model == m] - tmp$P_mean[tmp$Model == m][1]
        if(selType == "Divergent"){
          tmp2$P_mean[tmp$Model == m] <- tmp2$P_mean[tmp2$Model == m] - tmp2$P_mean[tmp2$Model == m][1]
        }
      }

      data.tmp$Generation <- tmp$Generation
      data.tmp$Effect <- tmp$Model
      data.tmp$Value <- switch (effect,
                                P = tmp$P_mean,
                                gv = tmp$gv_mean,
                                mv = tmp$mv_mean
      )

      data.tmp$CI_min <- switch(effect,
                                P = (tmp$P_mean - 1.96*tmp$P_sd/sqrt(iter)),
                                gv = (tmp$gv_mean - 1.96*tmp$gv_sd/sqrt(iter)),
                                mv = (tmp$mv_mean - 1.96*tmp$mv_sd/sqrt(iter))
      )
      data.tmp$CI_max <- switch(effect,
                                P = (tmp$P_mean + 1.96*tmp$P_sd/sqrt(iter)),
                                gv = (tmp$gv_mean + 1.96*tmp$gv_sd/sqrt(iter)),
                                mv = (tmp$mv_mean + 1.96*tmp$mv_sd/sqrt(iter))
      )

      if(selType == "Divergent"){
        data.tmp2$Generation <- tmp2$Generation
        data.tmp2$Effect <- tmp2$Model
        data.tmp2$Value <- switch (effect,
                                   P = tmp2$P_mean,
                                   gv = tmp2$gv_mean,
                                   mv = tmp2$mv_mean
        )

        data.tmp2$CI_min <- switch(effect,
                                   P = (tmp2$P_mean - 1.96*tmp2$P_sd/sqrt(iter)),
                                   gv = (tmp2$gv_mean - 1.96*tmp2$gv_sd/sqrt(iter)),
                                   mv = (tmp2$mv_mean - 1.96*tmp2$mv_sd/sqrt(iter))
        )
        data.tmp2$CI_max <- switch(effect,
                                   P = (tmp2$P_mean + 1.96*tmp2$P_sd/sqrt(iter)),
                                   gv = (tmp2$gv_mean + 1.96*tmp2$gv_sd/sqrt(iter)),
                                   mv = (tmp2$mv_mean + 1.96*tmp2$mv_sd/sqrt(iter))
        )
      }
    }

    #Unify datasets of divergent populations
    if(selType == "Divergent"){
      data.plot <- rbind(data.tmp,data.tmp2)
      x <- sum(data.tmp2$Value[data.tmp2$Generation == nyear],na.rm = TRUE)
      y <- sum(data.tmp$Value[data.tmp$Generation == nyear],na.rm = TRUE)

      if(x > y){
        data.plot$Population[data.plot$Population == "data2"] <- "High"
        data.plot$Population[data.plot$Population == "data"] <- "Low"
      }else{
        data.plot$Population[data.plot$Population == "data2"] <- "Low"
        data.plot$Population[data.plot$Population == "data"] <- "High"
      }

    }else{
      data.plot <- data.tmp
      data.plot$Population <- rep(selType,nrow(data.plot))
    }


    y.min <- min(data.plot$CI_min,na.rm = T) + min(data.plot$CI_min,na.rm = T)/10
    y.max <- max(data.plot$CI_max,na.rm = T) + max(data.plot$CI_max,na.rm = T)/10

    #### Asignar con switch dependiendo de si modelo o efecto
    data.plot$col<-paste(data.plot$Population,"_",data.plot$Effect,sep="")
    if(draw == "Model"){
      color_scale <- c("P" = "#AD70AD", "gv" = "#00a7c7", "mv" = "#8d1c1a")
      effect_order <- c("P", "gv", "mv")
    } else {
      color_scale <- c("G" = "#00a7c7", "M" = "#8d1c1a", "NMH" = "#708090",
                       "LMH" = "#4c83c8", "MMH" = "#835faa", "HMH" = "#9e3773")
      effect_order <- c("G", "M", "NMH", "LMH", "MMH", "HMH")
    }

    # Set the factor levels for Effect
    data.plot$Effect <- factor(data.plot$Effect, levels = effect_order)

    p <- ggplot(data = data.plot, aes(x = factor(Generation), y = Value, group = col, color = Effect, fill = Effect))

    if(type_se == "classic"){
      p <- p + geom_errorbar(aes(ymin = CI_min, ymax = CI_max), width = 0.2)
    }else{
      p <- p + geom_ribbon(aes(ymin = CI_min, ymax = CI_max), alpha = 0.3)
    }

    p <- p +  geom_line(aes(linetype = Population), linewidth = 0.5) +
      geom_point(size = 2, shape = 21) +
      scale_color_manual(values = color_scale, name = plot) +
      scale_fill_manual(values = color_scale, name = plot) +
      ylim(y.min, y.max) +
      scale_x_discrete(expand = c(0, 0)) +
      geom_hline(yintercept = 0, colour = "#2F2F28", linetype = "dashed", size = 0.1) +
      labs(x = "Generation") +
      theme_classic() +
      theme(
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.position = "right",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
      )
  }




  if(is.null(name)){
    if(plot == "boxplot"){
      name <- paste(plot,selType,"last_genetation",".tiff",sep="_")
    }else{
      if(draw=="Model"){
        name <- paste(plot,draw,selType,model,".tiff",sep="_")
      }else{
        name <- paste(plot,draw,selType,effect,".tiff",sep="_")
      }

    }

  }
  ggsave(paste(path,name, sep = ""),p,dpi = dpi, width = width, height = height, units = units)

  if(show == TRUE){
    p
  }
  return(p)
}


