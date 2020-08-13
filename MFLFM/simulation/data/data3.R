library(mvtnorm)
library(fda); library(KFAS); library(MCMCpack); 
#library(stochvol);
library(truncnorm) ; 
library(dlm)

getSplineInfo = function(tau, KNOTS=20, intKnots = NULL){
  
  allTaus = sort(unique(tau)) 	# all observation points
  a = min(allTaus)        	# lower endpoint
  b = max(allTaus)        	# upper endpoint
  
  numIntKnots = KNOTS 	# number of interior knots (= M in paper)
  if( is.null(intKnots) ){ intKnots = quantile(allTaus,seq(0,1,length= (numIntKnots+2))[-c(1,(numIntKnots+2))])}  	# interior knots 
  
  basis = create.bspline.basis(c(a,b),breaks=c(a,intKnots,b))						# spline basis
  blin = create.monomial.basis(c(a,b), nbasis=2) 								# linear basis
  
  Phi = bs(allTaus,knots=intKnots,degree=3,Boundary.knots=c(a,b),intercept=TRUE) 		# basis matrix
  Omega = eval.penalty(basis, Lfdobj=2, rng=c(a,b)) 					 		# derivative matrix
  Jkl = eval.penalty(basis, Lfdobj=0, rng=c(a,b))						 		# integral matrix
  
  # Now, transform these matrices to align with the linear/nonlinear decomposition of d_k:
  # The first two entries of d_k are the linear components, the remaining are nonlinear (see Wand and Ormerod, 2008)
  eigOmega = eigen(Omega)
  indsZ = 1:(numIntKnots+2)
  UZ = eigOmega$vectors[, indsZ] 			# Linear part
  LZ = t(t(UZ)/sqrt(eigOmega$values[indsZ])) 	# Nonlinear part
  
  # Basis matrices
  PhiL = cbind(1, allTaus)				# Linear part
  PhiN = Phi%*%LZ						# Nonlinear part
  Phi = cbind(PhiL, PhiN)
  
  return(Phi)
}


# Preda
TT <- 4
nts <- rep(30, TT)
te <- seq(0,1,length.out = 128) #seq(1,21,0.05)

# observation numbers
nobs1 <- nobs2 <- rep(length(te),sum(nts)) #sample(40:60, n, replace=TRUE)
nobs <- list(nobs1, nobs2)
truezs <- c( c(rep(1,10),rep(2,10),rep(3,10)), c(rep(1,5),rep(2,15),rep(3,10)), c(rep(2,15),rep(3,10),rep(4,5)), c(rep(2,10),rep(3,15),rep(5,5)) )


p <- 6
Xdat <- matrix(NA, p, sum(nts))


fval1 <- 2*log( te + 1/10 )
fval2 <- cos(5*te)+2
fval3 <- exp(te)
fval4 <- tan(te+1/4)
fval5 <- 3*sin( 3*te + 1/2 ) - 1
fval6 <- 3*te^3 + 1

Ymean1 <- Ymean2 <- matrix(NA, length(te), sum(nts))
for(i in 1:sum(nts)){ 
  
  zi <- truezs[i]
  
  tempu1 <- rnorm(1,    zi,  0.5)
  tempu2 <- rnorm(1, ifelse( (zi!=2&&zi!=5), zi, 0 ), 0.6) #rnorm(1, ifelse( (zi!=2||zi!=4), -zi, 0 ), 0.6)
  tempu3 <- rnorm(1, ifelse( (zi!=1&&zi!=4), zi, 0 ), 0.6)
  tempu4 <- rnorm(1,    -zi,  0.5)
  tempu5 <- rnorm(1, ifelse( (zi!=3&&zi!=4), zi, 0 ), 0.6)
  tempu6 <- rnorm(1,    zi*(-1)^zi,  1)#rnorm(1, ifelse( (zi!=1||zi!=4), -zi, 0 ), 0.6)
  
  Xdat[,i] <- c(tempu1, tempu2, tempu3, tempu4, tempu5, tempu6)
  
  Ymean1[,i] <- 3*te*zi + tempu1 * fval1 + tempu2 * fval2 + tempu3 * fval3
  Ymean2[,i] <- 2*zi + tempu4 * fval4 + tempu5 * fval5 + tempu6 * fval6
  
}


Ydat1 <- Ydat11 <- Ymean1 + matrix( rnorm(length(te)*sum(nts), 0, 1), length(te), sum(nts) )
Ydat2 <- Ydat22 <- Ymean2 + matrix( rnorm(length(te)*sum(nts), 0, 1), length(te), sum(nts) )

Ydat1 <- c(Ydat1)
Ydat2 <- c(Ydat2)
Ydat <- list(Ydat1, Ydat2)

tobs1 <- NULL; for(i in 1:sum(nts)){ tobs1 <- c( tobs1, te ) }
tobs2 <- NULL; for(i in 1:sum(nts)){ tobs2 <- c( tobs2, te ) }
tobs <- list(tobs1, tobs2)
cc <- 2





figurenum <- c("(a", "(b", "(c", "(d", "(e", "(f", "(g", "(h", "(i", "(j", "(k", "(l", "(m")

par(mfrow=c(2,3))
for(i in 1:6){

  DIM <- i
  YLIM <- quantile(Xdat[DIM,], c(0,1))

  plot(1, Xdat[DIM,1], ylim=YLIM, xlim=c(3/4,21/4), xlab='', ylab=paste("u",i,sep=''),
       main=paste(figurenum[i+4], ") Distribution of ", "U", i, sep=""), xaxt='n', cex.main=1.6, cex.lab=1.3)
  axis(1, at=c(1,2,3,4,5), labels=c("cluster 1", "cluster 2", "cluster 3","cluster 4","cluster 5"))
  points( 1+rnorm(sum(truezs == 1),0,1/15), Xdat[DIM,truezs == 1] )
  points( 2+rnorm(sum(truezs == 2),0,1/15), Xdat[DIM,truezs == 2] )
  points( 3+rnorm(sum(truezs == 3),0,1/15), Xdat[DIM,truezs == 3] )
  points( 4+rnorm(sum(truezs == 4),0,1/15), Xdat[DIM,truezs == 4] )
  points( 5+rnorm(sum(truezs == 5),0,1/15), Xdat[DIM,truezs == 5] )
}

dev.off()

par(mfrow=c(2,2))
YLIM <- quantile(Ymean1, c(0,1))
plot(te, Ymean1[,1], type='l', ylim=YLIM, ylab='y', xlab='t', 
     main = paste(figurenum[1], ") Mean Curves of dimension 1 for Simulation case (3)", sep=''), cex.main=1.4, cex.lab=1.3)
for(i in 1:sum(nts)){ lines(te, Ymean1[,i], col=truezs[i]) }


YLIM <- quantile(Ymean2, c(0,1))
plot(te, Ymean2[,1], type='l', ylim=YLIM, ylab='y', xlab='t', 
     main = paste(figurenum[2], ") Mean Curves of dimension 2 for Simulation case (3)", sep=''), cex.main=1.4, cex.lab=1.3)
for(i in 1:sum(nts)){ lines(te, Ymean2[,i], col=truezs[i]) }


YLIM <- quantile(Ydat11, c(0,1))
plot(te, Ydat11[,1], type='l', ylim=YLIM, ylab='y', xlab='t', 
     main = paste(figurenum[3], ") Observed Curves of dimension 1 for Simulation case (3)", sep=''), cex.main=1.4, cex.lab=1.3)
for(i in 1:sum(nts)){ lines(te, Ydat11[,i], col=truezs[i]) }


YLIM <- quantile(Ydat22, c(0,1))
plot(te, Ydat22[,1], type='l', ylim=YLIM, ylab='y', xlab='t', 
     main = paste(figurenum[4], ") Observed Curves of dimension 2 for Simulation case (3)", sep=''), cex.main=1.4, cex.lab=1.3)
for(i in 1:sum(nts)){ lines(te, Ydat22[,i], col=truezs[i]) }

dev.off()

######################
# radial basis 
######################
radialBasis <- function( t, Center, nu, p ){ rbasis <- rep(1,p); rbasis[2:p] <- exp( -nu * ( t-Center[-p] )^2); return(rbasis); }

# K2 <- 65
# fb <- create.fourier.basis(c(0,1), K2)
# B2 <- eval.basis(tobs[[2]], fb)

MakeBasis <- function(basis, TOBS, YDAT){
  
  if(basis == 1){
    
    K <- 65
    
    fb <- create.fourier.basis(c(0,1), K)
    B <- eval.basis(TOBS, fb)
    
  }else if(basis == 2){
    
    K <- 40
    
    # nu
    Nu <- 4
    center <- (0:(K-1))/(K-1);
    
    # B
    B <- matrix(0, length(TOBS), K)
    for(i in 1:length(TOBS)){ B[i,] <- radialBasis(TOBS[i], center, Nu, K) }
    
  }else if(basis == 3){
    
    K <- 20
    
    # B
    tempB <- getSplineInfo(TOBS, K)
    K <- ncol(tempB)
    
    Utobs <- sort(unique(TOBS))
    B <- matrix(0, length(TOBS), K)
    for(i in 1:length(TOBS)){ B[i,] <- tempB[which( Utobs %in% TOBS[i] ),] }
    
  }else{
    
    K <- 20
    
    Bbasis <- bs(YDAT, df = K, degree = 3)
    B <- Bbasis
  }
  
  kk <- K
  bb <- B
  
  return( list(kk, bb) )
}

################################################
BasisResult1 <- MakeBasis(3, tobs[[1]], Ydat[[1]]) # 1 2 3 4 
BasisResult2 <- MakeBasis(3, tobs[[2]], Ydat[[2]]) # 1 2 3 4
################################################

K1 <- BasisResult1[[1]]
B1 <- BasisResult1[[2]]

K2 <- BasisResult2[[1]]
B2 <- BasisResult2[[2]]

Bcs <- list(B1, B2)

#Be <- getSplineInfo(te, K-4)
n <- sum(nts)

Kcs <- c(K1,K2)
