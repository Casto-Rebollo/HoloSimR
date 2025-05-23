% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeM.R
\name{makeM}
\alias{makeM}
\title{Compute Microbiota Abundances}
\usage{
makeM(pop,
      sym = 0,
      founderM = NULL,
      acquiredSp = NULL,
      rndSeed = NULL,
      globalSP = NULL)
}
\arguments{
\item{pop}{An object of \code{\link{MapPop-class}} or a list from \code{\link{simBasePop}}.}

\item{sym}{A value of 1 to consider the symbiosis effect on microbiota. \emph{Default} is 0.}

\item{founderM}{A list from \code{\link{setFounderM}}. If \code{NULL}, it will be obtained from the global environment.}

\item{acquiredSp}{A matrix of the starting microbial species in the current population from \code{\link{acquiredSpecies}}. If \code{NULL}, it will be obtained from the global environment.}

\item{rndSeed}{A seed for random sampling. If it is \code{NULL}, a random seed will be set.}

\item{globalSP}{An object of class \code{\link{GlobalSP}}. If \code{NULL}, it will be obtained from the global environment.}
}
\value{
A \code{\link{MapPop-class}} object that includes microbiota abundance information (M) and the microbiota value (mv) for each individual, according to \code{\link{fillSp}}.
}
\description{
The microbiota abundance is computed according to the variance set by the user in the object \code{\link{GlobalSP}}.
}
\details{
The microbiota abundance is computed considering that the variability of each microbial species is due to:
\eqn{\text{var}_{\text{m}_{\text{k}}} = \text{var}_{\text{hg}_{\text{k}}} + \text{var}_{\text{s}_{\text{k}}} + \text{var}_{\text{e}_{\text{k}}}},
where \eqn{\text{var}_{\text{m}_{\text{k}}}} is the total abundance variance of species \eqn{k},
\eqn{\text{var}_{\text{hg}_{\text{k}}}} is the variance due to the host genetic effect on \eqn{k},
\eqn{\text{var}_{\text{s}_{\text{k}}}} is the variance of the symbiosis effect, and
\eqn{\text{var}_{\text{e}_{\text{k}}}} is the residual variance of species \eqn{k}.

The microbiota abundance is considered to be mediated by two components: a stable abundance generated by the host genetics and the environment (\eqn{\text{x}_{\text{ik}_{\text{0}}}}), and the abundance due to the symbiosis in this stable abundance:
\eqn{\text{x}_{\text{ik}} = (\text{x}_{\text{ik}})_{\text{0}} + (\text{x}_{\text{ik}})_{\text{s}}},

\eqn{\text{x}_{\text{ik}_{\text{0}}} = \mu_{\text{k}} + \sum_{j=1}^n \text{z}_{\text{ij}} \beta_{\text{jk}} + \text{e}_{\text{ik}}},
where \eqn{\text{x}_{\text{ik}_{\text{0}}}} is the stable abundance of individual \eqn{i} and species \eqn{k},
\eqn{\mu_{\text{k}}} is the mean abundance for species \eqn{k},
\eqn{\text{z}_{\text{ij}}} is the allele dosage of individual \eqn{i} for SNP \eqn{j},
\eqn{\beta_{\text{jk}}} is the host genetic effect of SNP \eqn{j} on the species abundance \eqn{k},
and \eqn{\text{e}_{\text{ik}}} is the residual. This is based on \cite{Pérez-Enciso et al. (2021)}.

\eqn{(\text{x}_{\text{ik}})_{\text{s}} = (\text{x}_{\text{i}})_{\text{0}}' \gamma_{\text{k}}},
where \eqn{(\text{x}_{\text{ik}})_{\text{s}}} is the abundance due to the symbiosis effect for individual \eqn{i} and species \eqn{k},
and \eqn{(\text{x}_{\text{i}})_{\text{0}}'} is the vector of the stable microbial abundance of individual \eqn{i}.
}
\examples{
# Set simulation global parameters
gSP <- GlobalSP$new(nPop = 1000, nyear = 10, nQTLchr = 100)

# Define the microbiota
gSP$setSpecies(nSpecies = 1000, nSp0 = 600, nSpEff = 100, symbiosis = c(0, 1))

# Define the trait
gSP$setTrait(meanP = 10, varP = 3, h2 = 0.10, m2 = 0.10)

# Create founder genetic
founderG <- setFounderG(globalSP = gSP)

# Set simulation genetic parameters
SP <- essentialSP(founder = founderG, minSnpFreq = 0.05, nSnpChr = 10000)

# Create founder microbiota
founderM <- setFounderM(globalSP = gSP)

# Create the founder population
founderPop <- newPop(founderG)

# Create the base population
basePop <- simBasePop(model = "NMH")

# Create a new population from the base population
Pop <- randCross(basePop$Pop, nCrosses = gSP$nDam, nProgeny = gSP$nSon)

# Make acquisition matrix of the new population
acquiredSp <- acquiredSpecies(pop = Pop)

# Compute the microbiota
Pop <- makeM(pop = Pop, globalSP = gSP)

}
\references{
Pérez-Enciso M, Zingaretti LM, Ramayo-Caldas Y, et al. Opportunities and limits of combining microbiome and genome data for complex trait prediction. Genet Sel Evol. 2021;53:65. \url{https://doi.org/10.1186/s12711-021-00658-7}
}
