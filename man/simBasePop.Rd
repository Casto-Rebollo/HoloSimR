% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simBasePop.R
\name{simBasePop}
\alias{simBasePop}
\title{Simulation base population}
\usage{
simBasePop(model,
           founderPop = NULL,
           founderM = NULL,
           globalSP = NULL,
           rndSeed=NULL,
           progressBar = FALSE)
}
\arguments{
\item{model}{a value for the model to simulate. Only one value is available according to \dQuote{G}, \dQuote{M}, \dQuote{NMH}, \dQuote{LMH}, \dQuote{MMH}, \dQuote{HMH}.}

\item{founderPop}{Founder population to use for founding the base population. If it is \code{NULL}, it will be taken from the environment.}

\item{founderM}{a list of \code{\link{setFounderM}}. If \code{NULL}, it will be retrieved from the global environment.}

\item{globalSP}{an object of class \code{\link{GlobalSP}}. If \code{NULL}, it will be retrieved from the global environment.}

\item{rndSeed}{rndSeed seed to initialise the random sampling. If it is \code{NULL}, it will be set randomly.}

\item{progressBar}{progressBar If TRUE, progress bar and text will be displayed. \emph{Default} \code{FALSE}}
}
\value{
a list including three dataframes:\enumerate{
\item{\strong{Pop}}: an object of \code{\link{MapPop-class}} the base population
\item{\strong{Beta_scale}}: a matrix with the scaled \eqn{\beta} effect
\item{\strong{Symbiosis_scale}}: a matrix with the scaled \eqn{gamma} effect. It is \code{NULL} if there is no symbiosis.
\item{\strong{varG}}: a vector with the genetic variance of each microbial species.
\item{\strong{varS}}: a vector with the symbiosis variance of each microbial species. It is \code{NULL} if there is no symbiosis.
\item{\strong{varE}}: a vector with the environmental variance of each microbial species.
\item{\strong{w_scale0}}: a vector with the scaled \eqn{\omega} effect without taking the symbiosis effect into account.
\item{\strong{w_scale1}}: a vector with the scaled \eqn{\omega} effect including the symbiosis effect. It is \code{NULL} if there is no symbiosis.
}
}
\description{
Simulation of the base population according to the parameters specified in the \code{\link{GlobalSP}} and the genome and microbiota created for the founder population with \code{\link{setFounderG}} and \code{\link{setFounderM}}. This function sets the base population, scales the host genetic (\eqn{\beta}), symbiosis (\eqn{\gamma}) and microbiota (\eqn{\omega}) effects according to the variance specified by the user, and simulates the microbiota in the base population.
}
\details{
The \eqn{\beta} (host genetic) and \eqn{\gamma} (symbiosis) effect on the microbiota were scaled according to the variance set by the user in the object \code{\link{GlobalSP}}.
The scale values were obtained by applying a scale factor follows \eqn{\sqrt{\frac{\text{var}_{\text{set}}}{\text{var}_{\text{sim}}}}} where \eqn{\text{var}_{\text{set}}} is the variance set according the parameters on \code{\link{GlobalSP}}, and \eqn{\text{var}_{\text{set}}} is the variance on the simulated matrix (\code{\link{setFounderM}}), to be adjusted.
}
\examples{
#Set simulation global parameters
gSP = GlobalSP$new(nPop = 1000, nyear = 10, nQTLchr = 100)

#Defining the microbiota
gSP$setSpecies(nSpecies = 1000, nSp0 = 600, nSpEff = 100, symbiosis = c(0,1))

#Defining trait
gSP$setTrait(meanP = 10,varP=3,h2=0.10,m2=0.10)

#Create founder genetic
founderG = setFounderG(globalSP = gSP)

#Set simulation genetic parameters
SP <- essentialSP(founder = founderG,
                  minSnpFreq = 0.05,rndSeed=NULL, nSnpChr = 10000)

#Create founder microbiota
founderM <- setFounderM(globalSP = gSP)

#Create the founder population
founderPop <- newPop(founderG)

#Create the base population
basePop <- simBasePop(model = "NMH")

}
