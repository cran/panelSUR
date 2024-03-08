\name{prepareData}
\alias{prepareData}
\title{Prepare data for use}

\description{
This function prepares data that have to be used.
}

\usage{
prepareData(data=data,
            restrictions=NULL,
            eqlist=eqlist)
}

\arguments{
\item{eqlist}{a \code{list} containing the equations making up the SUR system. They should be object of the class "formula" and necessarly include the intercept,}
\item{restrictions}{a vector containing constraints on the equation coefficients, which should be expressed in the form "equation_name$variable_name". Any spaces should be excluded from the restrictions definition. If one of the constraints includes an intercept term, the variable_name will be simply 'const'. Only simple restrictions involving equality between two parameters are considered, and not linear combinations involving more than two parameters,}
\item{data}{a data frame of the class "pdata.frame" (mandatory).}
}

\value{
An object of class \code{prepareData}, which is a list of the following elements:
\item{eqlist}{list of the equations of the system,}
\item{neq}{number of the system equations,}
\item{varlist}{list of the system variables,}
\item{ncoeff}{number of the system coefficients,}
\item{sumreg}{position of the first variable of each equation, including the constant, in the ordered list of the variables of the system,}
\item{nconstr}{number of contraints,}
\item{constr}{a matrix with as many rows as constraints, and whose row elements indicate the position, in the sorted list of model variables, of the variables affected by each constraint,}
\item{nind}{total number of individuals,}
\item{nt}{total number of individuals observed in each period,}
\item{psur}{table reporting the number of times each individual is observed,}
\item{psurmax}{maximum number of times the individuals are observed in the panel,}
\item{tmax}{number of period included in the panel,}
\item{sumTi}{sum of squares of the numbers of times each individual is observed,}
\item{sumNt}{sum of squares of the numbers of individuals observed in each time period,}
\item{vectorTi}{vector containing the number of times each individual is observed,}
\item{sysdata}{subset of the original data frame containing only the variables used in the estimated system,}
\item{infoSample}{information on the \code{sysdata} dataframe obtained trougth the \code{pdim} command of the \code{plm} package.}
}
