\name{preliminaryEstimate}
\alias{preliminaryEstimate}
\title{Obtain preliminary system equation estimates}

\description{
This function aims to obtain the preliminary (single within one or two way) estimate of the system equations.
}

\usage{
preliminaryEstimate(modelFrame=modelFrame,
                    method=method)
}

\arguments{
\item{modelFrame}{an object of the class \code{prepareData},}
\item{method}{the estimation method to be used, one of "1wayWB", "2wayWB", or "2wayQUE".}
}

\value{
An object of class \code{preliminaryEstimate}, which is a list of the following elements:
\item{f1w}{centered residuals of the \emph{oneway within} estimation obtained for each single equation of the system,}
\item{f2w}{centered residuals of the \emph{twoways within} estimation obtained for each single equation of the system,}
\item{mi_f1w}{individual means of the centered \emph{oneway within} residuals,}
\item{mi_f2w}{individual means of the centered \emph{twoways within} residuals,}
\item{mt_f2w}{time means of the centered \emph{twoways within} residuals,}
\item{m_f1w}{mean of the centered \emph{oneway within} residuals,}
\item{reglist}{list of the regressor matrix of each equation,}
\item{reglist2}{list of the regressors data frame of each equation,}
\item{regnames}{a vector whose elements are the names of all the variables included in each equation of the system,}
\item{final_regnames}{a vector whose elements are the names of the variables considered in the system equations, taking into account only the first appearance of those affected by restrictions on the coefficients.}
}

