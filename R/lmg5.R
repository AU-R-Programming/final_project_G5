#' @title Linear regression function for final
#'
#' @description This package contains the basic functions to perform linear
#' regression without using any of the linear regression functions already
#' available
#' @param response A \code{vector}  of length n to denote response variable.
#' @param covariate A \code{matrix} of dimension n*p denoting covariate
#' @param alpha A \code{numeric} used to denote the significance level.
#' @return A \code{list} containing the following attributes:
#' \describe{
#'      \item{beta}{Estimated coefficients}
#'      \item{sigma2}{Estimated variance}
#'      \item{variance_beta}{variance of beta}
#'      \item{ci}{Confidence interval}
#'      \item{r2}{Goodness of fit}
#'      \item{Cp}{Mallows' Cp}
#'      \item{Fteststat}{F-test statistic}
#'      \item{pvalue}{p-value of F-test statistic}
#'      }
#' @author Group 5
#' @importFrom ggplot2 ggplot
#' @export
#' @examples
#' data(cars)
#' lmg5(cars$speed, cars$dist, alpha = 0.5)
lmg5 = function(response, covariates, alpha = 0.05) {

  library(ggplot2)

  #Function to minimize to achieve beta estimate

  beta_hat_func = function(y, x, beta){

    return((t(y-as.matrix(x)%*%beta))%*%(y-as.matrix(x)%*%beta))
  }

  # Make sure data formats are appropriate
  response <- as.vector(response)
  covariates <- as.matrix(covariates)

  # Define parameters
  n <- length(response)
  p <- dim(covariates)[2]
  df <- n - p

  # Estimate beta through optimization()
  x= covariates
  x<- cbind(1,x)
  covariates = x
  y= response
  beta_zero<- NA
  p<- dim(x)[2]
  beta_zero[1] <- mean(y)

  #continuing to fill in beta_zero

  for (i in 2:p) {
    beta_zero[i]<- cov(y,x[,i])/var(x[,i])
  }


  #make sure beta hat function works for y x and betazero

  beta_hat_func(y,x,beta_zero)

  betahat_optimize<- optim(beta_zero, function (z) beta_hat_func(y,x,z))
  beta.hat<- betahat_optimize $par

  # Estimate of the residual variance (sigma2) from Eq. (6.3)
  # Compute residuals
  resid <- response - covariates%*%as.matrix(beta.hat)
  sigma2.hat <- (1/df)*t(resid)%*%resid

  # Estimate of the variance of the estimated beta from Eq. (6.2)

  var.beta <- as.numeric(sigma2.hat)*solve(t(covariates)%*%covariates)

  # Estimate of the confidence interval based on alpha
  quant <- 1 - alpha/2
  ci.beta <- c(beta.hat - qnorm(p = quant)*sqrt(var.beta), beta.hat +
                 qnorm(p = quant)*sqrt(var.beta))
  # Estimate yhat
  yhat<- x %*% as.matrix(beta.hat)

  #calulate SSE
  SSE<- sum((y-yhat)^2)


  #calculate SST
  SST<- sum((y-mean(y))^2)

  #Calulate R^2

  Rsquared<- 1-(SSE/SST)

  #calculate Cp

  cp<- SSE + 2*p*sigma2.hat

  #calculate DFM
  DFM<- p-1

  #calculate DFE

  DFE<- n-p

  #calculate SSM

  SSM<- sum((yhat-mean(y))^2)

  #calculate MSM

  MSM<- SSM/DFM

  #calculate MSE

  MSE<- SSE/ DFE

  #calculate Fstar
  fstar<- MSM/MSE

  #calculate fpvalue
  fpvalue<- pf(fstar,DFM, DFE)

  #Plot residuals vs yhat
  resid_vs_fitted <- ggplot() + geom_point(aes(x= yhat, y= resid))

  plot(resid_vs_fitted)

  # Plot qq-plot of residuals
  qqnorm(resid)

  # Histogram of residuals
  hist(resid)

  # Return all estimated values
  return(list(beta = beta.hat, sigma2 = sigma2.hat,
              variance_beta = var.beta, ci = ci.beta, r2= Rsquared, Cp= cp, Fteststat= fstar, pvalue= fpvalue ))
}








