\name{SURest}
\alias{SURest}

\title{EC SUR System Models Estimation on (Unbalanced) Panel Data}

\description{
\code{SURest} is used to estimate one-way and two-way SUR systems on unbalanced panel data by GLS
    estimator also allowing cross-equation restrictions.
}

\usage{
SURest(data = data,
       eqlist = eqlist,
       restrictions = NULL,
       method="1wayWB")
}

\arguments{
\item{eqlist}{a \code{list} containing the equations making up the SUR system. They should be object of the class "formula" and necessarly include the intercept,}
\item{restrictions}{a vector containing constraints on the equation coefficients, which should be expressed in the form "equation_name$variable_name". Any spaces should be excluded from the restrictions definition. If one of the constraints includes an intercept term, the variable_name will be simply 'const'. Only simple restrictions involving equality between two parameters are considered, and not linear combinations involving more than two parameters,}
\item{method}{the estimation method to be used, one of "1wayWB", "2wayWB", or "2wayQUE" (see details),}
\item{data}{a data frame of the class "pdata.frame" (mandatory).}
}

\details{
\code{SURest} is a function for the GLS estimation of SUR system
models on (unbalanced) panel data.  It supports the following estimation methods: one-way error component procedure based on the Biorn (2004)’s procedure (\code{1wayWB}), two-way error component procedure based on the Biorn (2004)’s procedure (\code{2wayWB}), and the two-way QUE estimation procedure by Platoni \emph{et al.} (2012) (\code{2wayQUE}).
}

\value{
An object of class \code{SURest}, which is a list of the following elements:
\item{Sigma_u}{remainder error variance-covariance matrix,}
\item{Sigma_mu}{individual error variance-covariance matrix,}
\item{Sigma_nu}{time error variance-covariance matrix,}
\item{varnames}{a vector whose elements are the names of the variables considered in the system equations, taking into account only the first appearance of those affected by restrictions on the coefficients,}
\item{Estimate}{vector of the coefficient estimates of the system equations, taking into account only the first appearance of those affected by restrictions,}
\item{std_error}{vector of the standard errors of the coefficient estimates, taking into account only the first appearance of those affected by restrictions,}
\item{tstat}{vector of the t-statistics associated to the coefficient estimates, taking into account only the first appearance of those affected by restrictions,}
\item{pvalue}{vector of the p-values associated to the t-statistics,}
\item{infoSample}{information on the considered dataset obtained trougth the \code{pdim} command of the \code{plm} package,}
\item{neq}{number of the system equations,}
\item{Rsquared}{list of R-squared obtained for each equation of the estimated system,}
\item{method}{method choosen for the system estimation.}
}

\references{
Biorn E, (2004),
\emph{Regression Systems for Unbalanced Panel Data: a Stepwise Maximum Likelihood Procedure},
Journal of Econometrics, 122(2), 181--291.

Platoni S, Sckokai P, Moro D (2012),
\emph{A Note on Two-way ECM Estimation of SUR Systems on Unbalanced Panel Data},
Econometric Reviews, 31(2), 119--141.
}


\examples{
data("SURdata", package="panelSUR")

## Data preparation
library(plm)
datap <- pdata.frame(data, index=c("IND", "TIME"))

## Equations specification
eq1<-Y1~X1+X2
eq2<-Y2~X1+X2+X3
eq3<-Y3~X2+X3
eqlist<-c(eq1,eq2,eq3)

## Constraints specification
constraints<-c("eq1$X2=eq2$X1","eq2$X3=eq3$X2")

## System estimation
mod1<-SURest(eqlist=eqlist,restrictions=constraints,method="2wayQUE",data=datap)
}
