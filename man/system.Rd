\name{system}
\alias{system}
\title{Build and solve system for beta coefficient estimates}

\description{
This function aims to built and solve the system in order to obtain beta coefficient
estimates.
}

\usage{
system(modelFrame=modelFrame,
       firstEstimate=firstEstimate,
       matrices=matrices)
}

\arguments{
\item{modelFrame}{an object of the class \code{prepareData},}
\item{firstEstimate}{an object of the class \code{preliminaryEstimate},}
\item{matrices}{an object of the class \code{obtainSigmas}.}
}

\value{
An object of class \code{system}, which is a list of the following elements:
\item{BsurQ}{vector of the coefficient estimates of the system equations, taking into account only the first appearance of those affected by restrictions,}
\item{std_error}{vector of the standard errors of the coefficient estimates, taking into account only the first appearance of those affected by restrictions,}
\item{t_stat}{vector of the t-statistics associated to the coefficient estimates, taking into account only the first appearance of those affected by restrictions,}
\item{p_value}{vector of the p-values associated to the t-statistics.}
}

