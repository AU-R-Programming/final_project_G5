% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lmg5.R
\name{lmg5}
\alias{lmg5}
\title{Linear regression function for final}
\usage{
lmg5(response, covariates, alpha = 0.05)
}
\arguments{
\item{response}{A \code{vector}  of length n to denote response variable.}

\item{alpha}{A \code{numeric} used to denote the significance level.}

\item{covariate}{A \code{matrix} of dimension n*p denoting covariate}
}
\value{
A \code{list} containing the following attributes:
\describe{
     \item{beta}{Estimated coefficients}
     \item{sigma2}{Estimated variance}
     \item{variance_beta}{variance of beta}
     \item{ci}{Confidence interval}
     \item{r2}{Goodness of fit}
     \item{Cp}{Mallows' Cp}
     \item{Fteststat}{F-test statistic}
     \item{pvalue}{p-value of F-test statistic}
     }
}
\description{
This package contains the basic functions to perform linear
regression without using any of the linear regression functions already
available
}
\examples{
data(cars)
lmg5(cars$speed, cars$dist, alpha = 0.5)
}
\author{
Group 5
}
