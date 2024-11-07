#' Summary of a generation of selection
#'
#' @description
#' Extract the mean, standard deviation and range of genetic, microbiota and phenotypic valuesof the scenarios, taking into account all iterations performed.
#'
#' @param data data.frame from \code{retValues} object.
#' @param g generation for calculating the statistics
#'
#' @return a data.frame with the following columns: \enumerate{
#'   \item{\strong{Generation}}: Generation computed
#'   \item{\strong{Model}}: Scenario computed
#'   \item{\strong{Gv_mean}}: Model computed
#'   \item{\strong{Gv_sd}}: Iteration performed. In case the selection will be replicated
#'   \item{\strong{Gv_range}}: Average of the population phenotype
#'   \item{\strong{Gvba_mean}}: Mean of the microbiota value for the breeding animals
#'   \item{\strong{Gvba_sd}}: Standard deviation of the microbiota value for the breeding animals
#'   \item{\strong{Gvba_range}}: Range of the microbiota value for the breeding animals
#'   \item{\strong{Mv_mean}}: Variance of the phenotype
#'   \item{\strong{Mv_sd}}: Average of the additive genetic value
#'   \item{\strong{Mv_range}}: Genetic variance
#'   \item{\strong{Mvba_mean}}: Mean of the microbiota value for the breeding animals
#'   \item{\strong{Mvba_sd}}: Standard deviation of the microbiota value for the breeding animals
#'   \item{\strong{Mvba_range}}: Range of the microbiota value for the breeding animals
#'   \item{\strong{P_mean}}: phenotypic mean
#'   \item{\strong{P_sd}}: standard deviation of the phenotype
#'   \item{\strong{P_range}}: range of the phenotype
#'   \item{\strong{Pba_mean}}: Mean of the phenotypic value for the breeding animals
#'   \item{\strong{Pba_sd}}: Standard deviation of the phenotypic value for the breeding animals
#'   \item{\strong{Pba_range}}: Range of the phenotypic value for the breeding animals
#'   }
#' @export
#' @usage summarySel(data, g)
#'
#' @examples
#' \dontrun{
#' summary_lastg <- summarySel(selection$Report$PopHigh$Symbiosis, globalSP$nyear)
#' }
summarySel <- function(data, g){

  summary <- data.frame(Model = character(),
                        gv_mean = numeric(),
                        gv_sd = numeric(),
                        gv_range = numeric(),
                        mv_mean = numeric(),
                        mv_sd = numeric(),
                        mv_range = numeric(),
                        P_mean = numeric(),
                        P_sd = numeric(),
                        P_range = numeric())

  for (g_iter in g) {
    count = 0
    data_tmp <- data[data$Generation==g_iter,]
    n <- length(unique(data_tmp$Model))
    tmp <- data.frame(Model = character(n),
                          gv_mean = numeric(n),
                          gv_sd = numeric(n),
                          gv_range = numeric(n),
                          mv_mean = numeric(n),
                          mv_sd = numeric(n),
                          mv_range = numeric(n),
                          P_mean = numeric(n),
                          P_sd = numeric(n),
                          P_range = numeric(n))

    for (model in unique(data_tmp$Model)) {
      count <- count + 1
      tmp$Generation[count] <- g_iter
      tmp$Model[count] <- model
      tmp$gv_mean[count] <-  mean(data_tmp$gv[data_tmp$Model==model])
      tmp$gv_sd[count] <-  sd(data_tmp$gv[data_tmp$Model==model])
      tmp$gv_range[count] <-  paste(range(data_tmp$gv[data_tmp$Model==model]),collapse = " ")

      tmp$mv_mean[count] <-  mean(data_tmp$mv[data_tmp$Model==model])
      tmp$mv_sd[count] <-  sd(data_tmp$mv[data_tmp$Model==model])
      tmp$mv_range[count] <-  paste(range(data_tmp$mv[data_tmp$Model==model]),collapse = " ")

      tmp$P_mean[count] <-  mean(data_tmp$P[data_tmp$Model==model])
      tmp$P_sd[count] <-  sd(data_tmp$P[data_tmp$Model==model])
      tmp$P_range[count] <-  paste(range(data_tmp$P[data_tmp$Model==model]),collapse = " ")
    }
    summary <- rbind(summary,tmp)

  }


  return(summary)
}


