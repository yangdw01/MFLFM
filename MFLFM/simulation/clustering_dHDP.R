####################################################################################
# Nonparametric Bayesian latent factor model for multivariate functional data
####################################################################################
# by Daewon Yang
####################################################################################
# 2020.04.06
####################################################################################
# simulation study 3
####################################################################################
# dynamic clustering model using dDHP
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
filename <- "data\\data3.R"


#####################
# source
#####################
sourceCpp("functions\\fast_dHDP_MFLFM.cpp")


#####################
# main
#####################

# iteration
niter <- 25000; 
burn <- 20000
thin <- 1
burnthin <- seq(burn+1, niter, thin)-1

NUM <- 100 # 100
n <- 120

MFLFM0result <- MFLFM1result <- matrix(0, NUM, n)


# main function
mainex <- function(){
  
  Zji_result <- matrix(0, length(burnthin), sum(nts))
  
  for(wert in 1:niter){
    
    cat(wert, "th iteration \n")
    
    ### part 1
    tempzv = getZandV( M, alpha, beta, tempKpq, ppi, llamb, ZETA, a_u, b_u, a_w, b_w, zzji, ZZ, VV, Eta, ums, Sigmas, mm, wms, Zmu )
    ZZ = tempzv[1:tempKpq,]
    VV = tempzv[1:tempKpq+tempKpq,]
    Eta = t( tempzv[1:sum(nts)+2*tempKpq,] )
    ums = c( tempzv[2*tempKpq+sum(nts)+1,] )
    mm = c( tempzv[2*tempKpq+sum(nts)+2,] )
    wms = c( tempzv[2*tempKpq+sum(nts)+3,] )
    Zmu = t( tempzv[1:M+2*tempKpq+sum(nts)+3,] )
    
    LFLM = ZZ * VV
    r = nrow(Eta)
    
    ums = getUMS(r, ZZ, VV, a_u, b_u);
    Eta = getEta(r, LFLM, Sigmas, ZETA, zzji, Zmu );
    
    alpha = getAlpha(ZZ, containBETA, tempKpq, beta, a_alpha, b_alpha, Hp);
    beta = getBeta(ZZ, containBETA, tempKpq, a_beta, b_beta, beta, alpha);
    
    ### part 2
    til_wl = UPDATE__til_wl( nts, rrji, aw, bw )
    wjl = UPDATE__wjl( length(nts), til_w0, til_wl )
    til_pilk = UPDATE__til_pilk( nts, M, betak, alpha0j, rrji, zzji, eps )
    pilk = UPDATE__pilk( length(nts), M, til_pilk )
    rrji = c(UPDATE__rji( nts, wjl, pilk, zzji, t(Eta), t(Zmu) ))
    zzji = c(UPDATE__zji( nts, pilk, rrji, t(Eta), t(Zmu) ))
    alpha0j = UPDATE__alpha0j( pilk, c0, d0 )
    betak = UPDATE__betak(M, nts, eps, zzji, gam)
    gam = UPDATE__gamma( betak, gam01, gam02 )
    Zmu = getZmu(zzji, Eta, r, wms, mm, M);
    mm = getmm(M, r, wms, Zmu);
    wms = getwms(M, Zmu, a_w, b_w, mm, r);
    
    ### part 3
    Sigmas = getSigmas( tempKpq, a_sigma, b_sigma, Eta, ZZ, VV, ZETA );
    
    if( containCAT == 1 ){ gis = getgis( tempKp, qind1, qind, qdims, Zdat, gis, LFLM, Eta, Sigmas, qq ); }
    if( cc>1 ){
      
      tempTHETA = getTHETA( cc, kkind, nobs, B, Sigmas, KKs, psi, LFLM, Eta, Ydat );
      psi = getpsi( a_psi, b_psi, cc, nobs, KKs, tempTHETA, Ydat, B );
      
      if(containCOV == 1){
        
        ZETA = rbind(tempTHETA, Xdat)
        
      }else{
        
        ZETA = tempTHETA
      }
      
    }else{
      
      tempTHETA = getTHETA0( nobs, B, Sigmas, tempK, psi0, LFLM, Eta, Ydat );
      psi0 = getpsi0( a_psi, b_psi, nobs, tempTHETA, Ydat, B );
      
      if(containCOV == 1){
        
        ZETA = rbind(tempTHETA, Xdat)
        
      }else{
        
        ZETA = tempTHETA
      }
    }
    
    if( containCAT == 1 ){ ZETA = rbind( ZETA, gis ); }
    
    
    if( wert %in% burnthin ){
      
      Zji_result[which(burnthin==wert),] = zzji
    }
  }
  
  clust_result <- Dahlclust(Zji_result)
  clusters <- clust_result[,1]
  
  return(clusters)
}




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
  containMULT <- 1
  
  source("functions\\preprocess_dHDP.R")
  
  # clustering
  result <- mainex()
  
  # result
  MFLFM0result[qwe,] <- result
}



##### MFLFM - functional data + covariate 
for(qwe in 1:NUM){

  set.seed(qwe)

  # data generation
  source(filename)

  # data preprocessing
  containCOV <- 1
  containCAT <- 0
  containBETA <- 1
  containMULT <- 0
  containMULT <- 1
  
  source("functions\\preprocess_dHDP.R")

  # clustering
  result <- mainex()

  # result
  MFLFM1result[qwe,] <- result
}


###################
# result
###################
true_cluster <- c( c(rep(1,10),rep(2,10),rep(3,10)), c(rep(1,5),rep(2,15),rep(3,10)), c(rep(2,15),rep(3,10),rep(4,5)), c(rep(2,10),rep(3,15),rep(5,5)) )

rate1 <- CCRcomplete(MFLFM1result, true_cluster)
rate0 <- CCRcomplete(MFLFM0result, true_cluster)


rate <- cbind(rate1, rate0)
colnames(rate) <- c("MF-LFM1", "MF-LFM0")

boxplot(rate, main = "", ylab="Correct Classification Rate (%)")
