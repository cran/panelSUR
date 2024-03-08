\name{obtainSigmas}
\alias{obtainSigmas}
\title{Compute errors' variance-covariance matrices}

\description{
This function aims to obtain the errors' variance-covariance matrices.
}

\usage{
obtainSigmas(modelFrame=modelFrame,
             firstEstimate=firstEstimate,
             method=method)
}

\arguments{
\item{modelFrame}{an object of the class \code{prepareData},}
\item{firstEstimate}{an object of the class \code{preliminaryEstimate},}
\item{method}{the estimation method to be used, one of "1wayWB", "2wayWB", or "2wayQUE".}
}

\value{
An object of class \code{obtainSigmas}, which is a list of the following elements:
\item{Sigma_u}{remainder error variance-covariance matrix,}
\item{Sigma_mu}{individual error variance-covariance matrix,}
\item{Sigma_nu}{time error variance-covariance matrix.}
}

