\name{printSUR}
\alias{printSUR}
\title{Print summary of estimated equation system}

\description{
This function prints a summary of the estimated equation system.
}

\usage{
printSUR(object)
}

\arguments{
\item{object}{an object of class \code{SURest}.}
}

\value{
No values are returned from the \code{printSUR} function. However, when called, it generates a visual output in the console, consisting of a formatted table containing the results of the SUR estimation and other relevant information.
}

\examples{
data("SURdata", package="panelSUR")

## Data preparation
library(plm)
datap <- pdata.frame(data, index=c("IND", "TIME"))

## Equations specification
eq1<-Y1~X1+X2
eq2<-Y2~X1+X2+X3
eqlist<-c(eq1,eq2)

## System estimation
mod1<-SURest(eqlist=eqlist,method="1wayWB",data=datap)

## Summary of estimation results
printSUR(mod1)
}


\keyword{models}
