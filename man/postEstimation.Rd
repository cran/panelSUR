\name{postEstimation}
\alias{postEstimation}
\title{Compute post-estimation indicators}

\description{
This function aims to compute some post estimation indicators.
}

\usage{
postEstimation(modelFrame=modelFrame,
               firstEstimate=firstEstimate,
               system=system)
}

\arguments{
\item{modelFrame}{an object of the class \code{prepareData},}
\item{firstEstimate}{an object of the class \code{preliminaryEstimate},}
\item{system}{an object of the class \code{system}.}
}

\value{
An object of class \code{postEstimation}, which is the list of R-squared obtained for each equation of the estimated system.
}

