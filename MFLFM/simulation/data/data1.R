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
n <- 40
n1 <- 15
n2 <- 15
n3 <- 10
te <- seq(0,1,length.out = 128) #seq(1,21,0.05)

# observation numbers
nobs <- rep(length(te),n) #sample(40:60, n, replace=TRUE)

p <- 3
Xdat <- matrix(NA, p, n)

q <-2
Zdat <- matrix(NA, q, n)

fval1 <- exp(te)
fval2 <- cos(5*te^2 + 1)
fval3 <- log( sin(10*te)+2 )

Ymean <- matrix(NA, length(te), n)
for(i in 1:n){ 
  
  if(i <= n1){
    
    tempu1 <- rnorm(1,1,0.3)
    tempu2 <- rnorm(1,0,0.4)
    tempu3 <- rnorm(1,0,0.6)
    
    tempg11 <- rnorm(1,1,1)
    
    tempg21 <- rnorm(1,-1,1)
    tempg22 <- rnorm(1,1,1)
    
    Xdat[,i] <- c(tempu1, tempu2, tempu3)
    Zdat[,i] <- c( which.max(c(tempg11, 0)), which.max(c(tempg21, tempg22, 0)) )
    
    Ymean[,i] <- (1*3/2 + tempg11/4) *te + (tempu1+tempg22/5) * fval1 + (tempu2+tempg21/4) * fval2 + tempu3 * fval3
    
  }else if(i <= n1+n2){ 
    
    tempu1 <- rnorm(1,2,0.6)
    tempu2 <- rnorm(1,2,0.4)
    tempu3 <- rnorm(1,0,0.6)
    
    tempg11 <- rnorm(1,-2,1)
    
    tempg21 <- rnorm(1,2,1)
    tempg22 <- rnorm(1,-1,1)
    
    Xdat[,i] <- c(tempu1, tempu2, tempu3)
    Zdat[,i] <- c( which.max(c(tempg11, 0)), which.max(c(tempg21, tempg22, 0)) )
    
    Ymean[,i] <- (2*3/2 + tempg11/4)*te + (tempu1+tempg22/5) * fval1 + (tempu2+tempg21/4) * fval2 + tempu3 * fval3
    
  }else{
    
    tempu1 <- rnorm(1,3,0.9)
    tempu2 <- rnorm(1,0,0.4)
    tempu3 <- rnorm(1,3,0.6)
    
    tempg11 <- rnorm(1,3,1)
    
    tempg21 <- rnorm(1,-3,1)
    tempg22 <- rnorm(1,-1,1)
    
    Xdat[,i] <- c(tempu1, tempu2, tempu3)
    Zdat[,i] <- c( which.max(c(tempg11, 0)), which.max(c(tempg21, tempg22, 0)) )
    
    Ymean[,i] <- (3*3/2 + tempg11/4)*te + (tempu1+tempg22/5) * fval1 + (tempu2+tempg21/4) * fval2 + tempu3 * fval3
  }
}

figurenum <- c("(a", "(b", "(c", "(d", "(e", "(f", "(g", "(h", "(i", "(j", "(k", "(l", "(m")
##### plot check
par(mfrow=c(3,2))

### 1
YLIM <- quantile(Ymean, c(0,1))
plot(te, Ymean[,1], type='l', ylim=YLIM, ylab='y', xlab='t', 
     main = paste(figurenum[1], ") Mean Curves for Simulation case (1)", sep=''), cex.main=1.6, cex.lab=1.3)
for(i in 1:n1){ lines(te, Ymean[,i]) }
for(i in (1:n2)+n1 ){ lines(te, Ymean[,i],col=2) }
for(i in (1:n3)+n1+n2 ){ lines(te, Ymean[,i], col=4) }

### 2
Ydat <- Ymean + matrix( rnorm(length(te)*n, 0, 0.5), length(te), n )

Ydat1 <- Ydat
YLIM <- quantile(Ydat1, c(0,1))
plot(te, Ydat1[,1], type='l', ylim=YLIM, ylab='y', xlab='t', 
     main = paste(figurenum[2], ") Observed Curves for Simulation case (1)", sep=''), cex.main=1.6, cex.lab=1.3)
for(i in 1:n1){ lines(te, Ydat1[,i]) }
for(i in (1:n2)+n1 ){ lines(te, Ydat1[,i],col=2) }
for(i in (1:n3)+n1+n2 ){ lines(te, Ydat1[,i], col=4) }

### 3

par(mfrow=c(2,3))
for(i in 1:3){
  
  DIM <- i
  YLIM <- quantile(Xdat[DIM,], c(0,1))
  
  plot(1, Xdat[DIM,1], ylim=YLIM, xlim=c(3/4,13/4), xlab='', ylab=paste("u",i,sep=''),
       main=paste(figurenum[2+i], ") Distribution of ", "U", i,  sep=""), xaxt='n', cex.main=1.6, cex.lab=1.3)
  axis(1, at=c(1,2,3), labels=c("cluster 1", "cluster 2", "cluster 3"))
  points( 1+rnorm(n1,0,1/15), Xdat[DIM,1:n1] )
  points( 2+rnorm(n2,0,1/15), Xdat[DIM,1:n2+n1] )
  points( 3+rnorm(n3,0,1/15), Xdat[DIM,1:n3+n1+n2] )
}

for(j in 1:q){
  
  x <- table( Zdat[j,], rep(c(1,2,3),c(15,15,10)) )
  
  colnames(x) <- c("cluster 1", "cluster 2", "cluster 3")
  barplot( x, beside=TRUE, ylab='Frequency', col=terrain.colors(length(unique(Zdat[j,]))), legend = rownames(x),
           args.legend = list(bty = "n", x = "top", ncol = 3), ylim=c(0,20), 
           main = paste(figurenum[5+j], ") Frequencies of ", "z", j, sep="") , cex.main=1.6, cex.lab=1.3)
}

dev.off()




Ydat <- c(Ydat)

tobs <- NULL; for(i in 1:n){ tobs <- c( tobs, te ) }

######################
# radial basis 
######################
radialBasis <- function( t, Center, nu, p ){ rbasis <- rep(1,p); rbasis[2:p] <- exp( -nu * ( t-Center[-p] )^2); return(rbasis); }


# ############ basis 1
# K <- 65
# 
# fb <- create.fourier.basis(c(0,1), K)
# B <- eval.basis(tobs, fb)
# 
# # te <- 0:300/300
# # Be <- eval.basis(te, fb)

# ############ basis 2
# # K: # of Gaussian radial bases
# K <- 40
# 
# # nu
# Nu <- 4
# center <- (0:(K-1))/(K-1);
# 
# # B
# B <- matrix(0, length(tobs), K)
# for(i in 1:length(tobs)){ B[i,] <- radialBasis(tobs[i], center, Nu, K) }
# 
# # te <- 0:300/300
# # Be <- matrix(0, length(te), K)
# # for(i in 1:length(te)){ Be[i,] <- radialBasis(te[i], center, Nu, K) }

############ basis 3
K <- 20

# B
tempB <- getSplineInfo(tobs, K)
K <- ncol(tempB)

Utobs <- sort(unique(tobs))
B <- matrix(0, length(tobs), K)
for(i in 1:length(tobs)){ B[i,] <- tempB[which( Utobs %in% tobs[i] ),] }

# te <- 0:300/300
# 


# ############# basis 4
# # K: # of Gaussian radial bases
# K <- 20
# 
# Bbasis <- bs(Ydat, df = K, degree = 3)
# B <- Bbasis

# te <- 0:300/300
# Be <- bs(te, df=K, degree=3, knots = attr(Bbasis,"knots"), Boundary.knots = attr(Bbasis,"Boundary.knots"))