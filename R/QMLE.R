#################################################
## Code QMLE : Quantile Maximization Likelihood Estimation
## Yousri SLAOUI and Abir EL HAJ January 2021
## Clara SOLIER and Cyril PERRET January 2021
## abir_hajj_1993@hotmail.com
## Yousri.slaoui@math.univ-poitiers.fr
## caroline.bordes@univ-poitiers.fr
## cyril.perret@univ-poitiers.fr
## jean.dumoncel@univ-poitiers.fr
## The maintenance of the package is done by Jean Dumoncel
#################################################

## We need the following packages:

require(pracma)
require(stats)
require(nloptr)
require(fitdistrplus)
require(invgamma)
require(dlm)
require(gamlss.dist)


#exemple of use:

#for QMLE:
#3 algorithms are possible for the QMLE function (neldermead, fminsearch, nlminb)
#neldermead = "NEMD"
#fminsearch = "FMIN"
#nlminb = "NLMI"

#when you have chosen the algorithm you can use the QMLEEstim function
#y<-rexGAUS(n=100, mu = 500, sigma = 150, nu = 100)
#QMLEEstim(y, 'NEMD')
#then you have mu, simga and tau in the result variable
#you can select only on of then : result$mu, result$sigma, result$tau


#for BayesianExgaussian:
#you will have to give him the number of sample n and the dataset

#y<-rexGAUS(n=100, mu = 500, sigma = 150, nu = 100)
#BayesianExgaussian(100, y)
#As the QMLEEstim function, you can select the variable that you need
#result$mu, result$sigma, result$tau




## Function determining the quantile estimates (q1) and the number of elements in each quantile (N)
theQuantile<-function(y){
#this function can estimate q1 quantile and the number of elements in each quantile (N)
#its parameter is y the data vector
  
  #packages require
  # library(pracma)
  # library(stats)
  # library(nloptr)
  # library(fitdistrplus)

  
  
  qq=quantile(y, probs = seq(0, 1, 0.05), na.rm = FALSE, names = TRUE, type = 5)

  probs = seq(0, 1, 0.05) # increasing set of proportions corresponding to the cumulative probabilities for each quantile
  data1=sort(y); # data sorted in ascending order
  #WARNING Y MUST BE A VECTOR
  n=length(y); #length of data
  N=matrix(0,(length(qq)),1)
  N[1]=NaN #NAN = missing
  # Number of elements in each quantile
  for (j in 2:length(qq)){
    N[j]=(probs[j]-(probs[j-1]))*n} 

  I=matrix(0,(length(qq)-1),1)
  mI=matrix(0,(length(qq)-1),1)
  pI=matrix(0,(length(qq)-1),1)
  I[1]=NaN
  mI[1]=NaN
  pI[1]=NaN
  q1=matrix(0,length(qq),1)

  # Compute the quantile estimates
  for (j in 2:length(qq)-1){
    I[j]= qq[j]*n+1/2
    mI[j]= floor(I(j)) # The largest integer that is smaller than or equal to I(j)
    pI[j]=ceiling(I(j)) # the smallest integer equal to or above I(j)
    

    q1[j]=y[mI[j]]+(y[pI[j]]-y[mI[j]])*(I[j]-mI[j])
  
  }

  q1[1]=qq[1] #  extreme quantile
  q1[length(qq)]=qq[length(qq)] # extreme quantile

  out=list(q1=q1,N=N)# returns q1 and N values

}




#############################################################################################################################################################################

## Define the Quantile maximum likelihood function

## The QMLE function


QMLE<-function(y, lower=NULL, upper=NULL, start1=NULL, q1, N, func){
#determines the quantile maximum likelihood of the data set

   
  #packages require
  # library(pracma)
  # library(stats)
  # library(nloptr)
  # library(fitdistrplus)



  Lvalue=matrix(0,length(q1),1)

  # function calculating (-1)*quantile likelihood
  
  
  dexgauss<-function(y){
    d=dnorm(y)+dexp(y)  
    return(d)
  }
  
  
  objFun<-function(start1=start1){
    zero <- .Machine$double.eps^20
    for (j in 2:length(q1)){
      integ=integrate(dexgauss ,q1[j-1], q1[j])$value # Compute the bounded integral of the ex-gaussian pdf over the interval [q1[j-1], q1[j]]

      if(integ<zero){ Lvalue[j-1]=0} else
      {Lvalue[j-1]=(N[j]*log(integ))}
    }
    Lvalue[is.nan(Lvalue)]=0;
    return(-sum(Lvalue))
  }

  
  if(func == "NEMD"){
    neldermead(start1, objFun, lower, upper, nl.info = FALSE)$par # Compute the minimum of the (-1)*quantile likelihood
  }else if(func == "FMIN"){
    fminsearch(objFun, start1, minimize = TRUE, maxiter = 1000, tol = .Machine$double.eps^(2/3));
  }else if(func == "NLMI"){
    nlminb(start1, objFun, gradient = NULL, hessian = NULL, scale = 1, control = list(), lower = lower, upper = upper)
  }
}


################################################################################################################################################################################
## Application on simulated data

#'Ex-Gaussian Quantile Maximum Likelihood Estimate
#' @import pracma stats nloptr fitdistrplus
#' 
#' @description Estimates the mu, sigma, and tau parameters of an ex-Gaussian distribution. 3 different algorithms can be used : neldermead ('NEMD'), fminsearch ('FMIN') and nlminb ('NLMI').

#' 
#' @param y the data. Must be a vector with no missing values
#' @param func the function selected for the estimation. 'NEMD' for neldermead, 'FMIN' for fminsearch, and 'NLMI' for nlminb.
#'
#' @return  QMLEEstim() returns an object "valEstim" which is a list with components: estimates of mu, sigma, and tau
#' 
#' @references 
#' Brown, S., & Heathcote, A. (2003). QMLE: Fast, robust, and efficient estimation of distribution functions based on quantiles. \emph{Behavior Research Methods, Instruments, & Computers}, \strong{35}, 485-492. 
#' 
#' McCormack, P. D., & Wright, N. M. (1964). The positive skew observed in reaction time distributions. \emph{Canadian Journal of Psychology/Revue canadienne de psychologie}, \strong{18}, 43-51. 
#' 
#' Van Zandt, T. (2000). How to fit a response time distribution. \emph{Psychonomic Bulletin & Review}, \strong{7}, 424-465.
#' 
#' El Haj, A., Slaoui, Y., Solier, C., & Perret, C. (2021). Bayesian Estimation of The Ex-Gaussian distribution. \emph{Statistics, Optimization & Information Computing}, \strong{9(4)}, 809-819.
#' 
#' Gilks, W. R., Best, N. G., & Tan, K. K. C. (1995). Adaptive rejection Metropolis sampling within Gibbs sampling. \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, \strong{44}, 455-472.
#' 
#' @examples
#' library(gamlss.dist)
#' set.seed(2703)
#' data<-rexGAUS(n=100, mu = 500, sigma = 150, nu = 100)
#' QMLEEstim(data, 'NEMD')
#' 
#' 
#' @export
QMLEEstim<-function(y,func){
  if (is.vector(y) == F) {
    stop ("'y' must be a vector")
  } else if (is.numeric(y) == F | is.integer(y) == T) {
    stop ("'y' must be numeric")
  } else if (NA %in% y) {
    stop ("'y' must not contain NAs")
  } else if (length(y) < 20) {
    stop("'y' contains less than 20 values")
  } else if ((func %in% c("NEMD", "FMIN", "NLMI")) == F) {
    stop("this is not a valid argument, must be 'NEMD', 'FMIN', or 'NLMI'")
  }

  #packages require for the called function
   # library(pracma)
   # library(stats)
   # library(nloptr)
   # library(fitdistrplus)
   # library(gamlss.dist)


  ## Application

  ## Determing the quantiles q1 and the number of elements in each quantile N

  q1=theQuantile(y)$q1
  N=theQuantile(y)$N

  ##Lower and Upper bound of the parameters

  lower <- c(min(y), 0, 0)
  upper <- c(max(y), (max(y)-min(y)), (max(y)-min(y)))

  ##Initial values of the parameters

  start1 <-c(
    mu      = mean(y)-(sd(y)*0.8),
    sigma   = sqrt(var(y)-(sd(y)*0.8)^2),
    tau     = sd(y)*0.8
  )
  ## minimize the (-1)*quantile likelihood to estimate the parameters
  Estim<-QMLE(y, lower=lower, upper=upper, start1=start1, q1, N, func) ## the last line "Optimal value of controls" is the estimated values of the parameters mu, sigma, and tau respectively


  #fetching of the estimate values of mu, sigma and tau for each function estimate
  if(func == "FMIN"){
    #fminsearch
    valEstim=list(
     mu= Estim$xmin[1],
     sigma=Estim$xmin[2],
     tau=Estim$xmin[3]
     );
  }else if(func == "NEMD"){
    #nedelmead
    valEstim=list(
     mu= Estim[1],
     sigma=Estim[2],
     tau=Estim[3]
    );
  }else if(func == "NLMI"){
    #nlminb
    valEstim=list(
      mu=Estim$par["mu"],
      sigma=Estim$par["sigma"],
      tau=Estim$par["tau"]
    );
  }

  return(valEstim) #return the three estimate values for mu, sigma and tau
}




















###################################################################################################################################################################################################################################
## Bayesian Function (input: data size, data, number of Samples, burn-in | output: the three estimated parameters: mu, sigma and tau)


#' Bayesian Ex-gaussian Estimate
#' @import stats invgamma dlm fitdistrplus
#'
#' @description Estimates the mu, sigma, and tau parameters of an ex-Gaussian distribution using a bayesian method.
#' @param n the data size
#' @param x the data. Must be a vector, with no missing values
#' @param nSamples number of Samples
#' @param Ti burn-in
#'
#' @return BayesianExgaussian() returns an object "theta" which is a list with components: estimates of mu, sigma, and tau
#' 
#' @references 
#' Brown, S., & Heathcote, A. (2003). QMLE: Fast, robust, and efficient estimation of distribution functions based on quantiles. \emph{Behavior Research Methods, Instruments, & Computers}, \strong{35}, 485-492. 
#' 
#' McCormack, P. D., & Wright, N. M. (1964). The positive skew observed in reaction time distributions. \emph{Canadian Journal of Psychology/Revue canadienne de psychologie}, \strong{18}, 43-51. 
#' 
#' Van Zandt, T. (2000). How to fit a response time distribution. \emph{Psychonomic Bulletin & Review}, \strong{7}, 424-465.
#' 
#' El Haj, A., Slaoui, Y., Solier, C., & Perret, C. (2021). Bayesian Estimation of The Ex-Gaussian distribution. \emph{Statistics, Optimization & Information Computing}, \strong{9(4)}, 809-819.
#' 
#' Gilks, W. R., Best, N. G., & Tan, K. K. C. (1995). Adaptive rejection Metropolis sampling within Gibbs sampling. \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, \strong{44}, 455-472.
#' 
#' @examples
#' \donttest{library(gamlss.dist)
#' set.seed(2703)
#' data<-rexGAUS(n=100, mu = 500, sigma = 150, nu = 100)
#' BayesianExgaussian(n = 100, x = data)}
#' 
#' @export
#' 
#' 
BayesianExgaussian <-function(n,x,nSamples=5000,Ti=2500){
  if (is.vector(x) == F) {
    stop ("'x' must be a vector")
  } else if (is.numeric(x) == F) {
    stop ("'x' must be numeric")
  } else if (length(levels(as.factor(x))) < 4 ) {
    stop("'x' has less than 4 distinct values")
  } else if (length(x) < 19) {
    stop("'x' contains less than 19 values")
  } else if (NA %in% x) {
    stop("'x' contains NAs")
  } else if (n != length(x)) {
    warning ("'n' is not equal to 'x' length !")
  }
  #packages require
  # library(fitdistrplus)
  # library(invgamma)
  # library(dlm)
  # library(gamlss.dist)

  
  ## Initialization of hyper-parameters

  delta1 = 1e-04
  delta2 = 0.1
  tau1 = 0.1
  tau2 = 0.1

  ## Defining the vectors of the parameters

  mu <- rep(NA, nSamples)
  zeta <- rep(NA, nSamples) # zeta=(sigma^2)
  lambda <- rep(NA, nSamples) # lambda=1/tau

  ## Initialization of the parameters

  mu[1] = mean(x) - (sd(x) * 0.5)
  zeta[1] = (var(x) - (sd(x) * 0.5)^2)
  lambda[1] = 1/(sd(x) * 0.5)

  y1 = mean(mu) - (sd(mu) * 0.5)
  y2 = var(mu) - (sd(mu) * 0.5)^2

  ## Lower and upper bounds

  minnmu = min(x)-1000
  maxxmu = max(x) + 1000
  minnzeta = 0.01
  maxxzeta = ((max(x) - min(x))/4)^2
  minnlambda = 0.001/mean(x)
  maxxlambda = 100/((max(x) - min(x))/2)

  # The posterior functions

  mufunc <- function(mu1) {return(sum((mu1-x))*lambda[i-1]+sum(log(pnorm(((x-mu1)/(sqrt(zeta[i-1])))-(lambda[i-1]*sqrt(zeta[i-1])))))-(((mu1-y1)^2)/(2*y2^2)))}
  zetafunc <- function(zeta1) {return(((n*lambda[i-1]^2*zeta1)/(2))+sum(log(pnorm(((x-mu[i])/(sqrt(zeta1)))-(lambda[i-1]*sqrt(zeta1)))))-(delta1+1)*log(zeta1)-(delta2/zeta1))}
  lambdafunc <- function(lambda1) {return((n*log(lambda1))+(sum((mu[i]-x))*lambda1)+((n*lambda1^2*zeta[i])/(2))+sum(log(pnorm(((x-mu[i])/(sqrt(zeta[i])))-(lambda1*sqrt(zeta[i])))))+(tau1-1)*log(lambda1)-(tau2*lambda1))}

  ## Defining the vectors of the parameters

  for(i in 2:nSamples) {

    mu[i] <- arms(mu[i-1], mufunc,function(mu1) (mu1>minnmu)*(mu1<maxxmu), 1)
    zeta[i] <- arms(zeta[i-1], zetafunc,function(zeta1) (zeta1>minnzeta)*(zeta1<maxxzeta), 1)
    lambda[i] <- arms(lambda[i-1], lambdafunc,function(lambda1) (lambda1>minnlambda)*(lambda1<maxxlambda), 1)

  }

  #The estimation of the parameters

  muest=sum(mu[Ti:nSamples])/(nSamples-Ti)
  zetaest=sum(zeta[Ti:nSamples])/(nSamples-Ti)
  sigmaest=sqrt(zetaest)
  lambdaest=sum(lambda[Ti:nSamples])/(nSamples-Ti)
  tauest=1/lambdaest
  theta=list(mu=muest,sigma=sigmaest,tau=tauest)
  return(theta)
}


#######################################################################################################################################################################################################