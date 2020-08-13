####################################################################################
# Nonparametric Bayesian latent factor model for multivariate functional data
####################################################################################
# by Daewon Yang
####################################################################################
# 2020.04.06
####################################################################################
# simulation study 1 and 2
####################################################################################


#####################
# load library
#####################
library(mvtnorm)
library(MCMCpack)
library(plsgenomics)
library(lpSolve)
library(corrplot)
library(fda)
library(Funclustering)
library(funHDDC)
library(maps)
library(expm)
library(R.matlab)
library(Rcpp)
library(RcppArmadillo)



#####################
# data
#####################
DATA <- 1
# 1: univariate
# 2: multivariate

if( DATA == 1 ){
  
  filename <- "data\\data1.R"
  
}else if(DATA ==2){
  
  filename <- "data\\data2.R"
}



#####################
# source
#####################
sourceCpp("functions\\fastNEW.cpp")
source("functions\\CCR.R")

#####################
# main
#####################

# iteration
niter <- 25000; 
burn <- 20000
thin <- 1
burnthin <- seq(burn+1, niter, thin)-1

n <- 40
NUM <- 100 # 100

MFLFM0result <- MFLFM2result <- MFLFM1result <- matrix(0, NUM, n)



##### MFLFM - only functional data
for(qwe in 1:NUM){
  
  set.seed(qwe)
  
  # data generation
  source(filename)
  
  # data preprocessing
  containCOV <- 0
  containCAT <- 0
  containBETA <- 1
  containMULT <- 0
  if(DATA == 2){ containMULT <- 1 }
  
  source("functions\\preprocess.R")
  
  # clustering
  result <- getCLUST( burnthin, niter, eps,
                      n, M, tempKpq, tempKp, tempK, ppi, llamb, cc, containBETA, containCOV, containCAT,
                      Hp, kkind, nobs, B, KKs, Ydat, Xdat, Zdatstar, qq, qdims, qind, qind1,
                      alpha, beta, ZETA, LL, ZZ, VV, Eta, ums, Sigmas, mm, wms, Zmu, PP, nu,
                      psi, psi0, tempTHETA, gis,
                      a_u, b_u, a_w, b_w, a_alpha, b_alpha,
                      a_beta, b_beta, a_nu, b_nu, a_sigma, b_sigma, a_psi, b_psi )
  
  MFLFM0result[qwe,] <- result
}


##### MFLFM - functional data + covariate
for(qwe in 1:NUM){
  
  set.seed(qwe)
  
  # data generation
  source(filename)
  
  # data preprocessing
  containCOV <- 1
  containBETA <- 1
  containMULT <- 0
  if(DATA == 2){ containMULT <- 1 }
  
  source("functions\\preprocess.R")
  
  # clustering
  result <- getCLUST( burnthin, niter, eps,
                      n, M, tempKpq, tempKp, tempK, ppi, llamb, cc, containBETA, containCOV, containCAT,
                      Hp, kkind, nobs, B, KKs, Ydat, Xdat, Zdatstar, qq, qdims, qind, qind1,
                      alpha, beta, ZETA, LL, ZZ, VV, Eta, ums, Sigmas, mm, wms, Zmu, PP, nu,
                      psi, psi0, tempTHETA, gis,
                      a_u, b_u, a_w, b_w, a_alpha, b_alpha,
                      a_beta, b_beta, a_nu, b_nu, a_sigma, b_sigma, a_psi, b_psi )
  
  MFLFM1result[qwe,] <- result
}




###################
# result
###################
true_cluster <- c( rep(1,15), rep(2, 15), rep(3,10) )

rate1 <- CCRcomplete(MFLFM1result, true_cluster)
rate0 <- CCRcomplete(MFLFM0result, true_cluster)


rate <- cbind(rate1, rate0)
colnames(rate) <- c("MF-LFM1", "MF-LFM0")

boxplot(rate, main = "", ylab="Correct Classification Rate (%)")
