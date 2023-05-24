library(R2WinBUGS)
library(R2jags)

# Data used in JAGS:
#    X , # treatment variable
#    M, # an N*TT matrix containing observed mediator scores at time point 1,2,...,TT for a sample of N individuals
#    Y, # an N*TT matrix containing observed outcome scores at time point 1,2,...,TT for a sample of N individuals
#    TimeM,  # a column vector of length TT containing the mediator slope loadings at time point 1, 2, ..., TT 
#    TimeY,  # a column vector of length TT containing the outcome slope loadings at time point 1, 2, ..., TT 
#    Z, # an N*numZ matrix containing a total of "numZ" observed pretreatment covariates, where "numZ" is the number of observed pretreatment covariates 
#    mZ, # a column vector containing the mean of Z; if Z is centered, mZ = rep(0, numZ)
#    varZ, # a numZ*numZ diagonal matrix containing the variance-covariance matrix of Z; if numZ=1 and Z is standardized, varZ = diag(1, nrow=1)  

# JAGS model with the diffuse priors         
jagsm <- function(){

  for (i in 1:N) {
    for (j in 1:TT) {
      M[i, j] ~ dnorm(mM[i, j], prec_eM)
      mM[i, j] <- ISM[i, 1] + ISM[i, 2] * TimeM[j]
      Y[i, j] ~ dnorm(mY[i, j], prec_eY)
      mY[i, j] <- ISY[i, 1] + ISY[i, 2] * TimeY[j]
    }
  }
  for (i in whichX0) {
    ISM[i, 1:2] ~ dmnorm(mISM[i, 1:2], prec_vMctrl[1:2, 1:2])
  }
  for (i in whichX1) {
    ISM[i, 1:2] ~ dmnorm(mISM[i, 1:2], prec_vMtrt[1:2, 1:2])
  }
  for (i in 1:N) {
    mISM[i, 1] <- a0[1] + aZ[1:numZ, 1] * Z[i, 1:numZ] + 
      aX[1] * X[i]
    mISM[i, 2] <- a0[2] + aZ[1:numZ, 2] * Z[i, 1:numZ] + 
      aX[2] * X[i]
  }
  for(i in 1:N){
    
    ISY[i,1:2] ~ dmnorm( mISY[i,1:2], prec_vY[1:2,1:2])
    mISY[i,1]<-b0[1]+bZ[1:numZ, 1]*Z[i,1:numZ]+ bX[1]*X[i] +
      bIM[1]*ISM[i,1] + bSM[1]*ISM[i,2] +
      bIMSM[1]*ISM[i,1]*ISM[i,2] + 
      bXIM[1]*X[i]*ISM[i,1] + bXSM[1]*X[i]*ISM[i,2] + bXIMSM[1]*X[i]*ISM[i,1]*ISM[i,2] 
    
    mISY[i,2]<-b0[2]+bZ[1:numZ, 2]*Z[i,1:numZ]+ bX[2]*X[i] +
      bIM[2]*ISM[i,1] + bSM[2]*ISM[i,2] +
      bIMSM[2]*ISM[i,1]*ISM[i,2] + 
      bXIM[2]*X[i]*ISM[i,1] + bXSM[2]*X[i]*ISM[i,2] + bXIMSM[2]*X[i]*ISM[i,1]*ISM[i,2]
    
  }
  
  prec_eM ~ dgamma(0.001, 0.001)
  var_eM <- 1/prec_eM
  prec_vMctrl[1:2, 1:2] ~ dwish(R, 3)
  prec_vMtrt[1:2, 1:2] ~ dwish(R, 3)
  Psi_vMtrt <- inverse(prec_vMtrt)
  Psi_vMctrl <- inverse(prec_vMctrl)
  prec_eY ~ dgamma(0.001, 0.001)
  prec_vY[1:2, 1:2] ~ dwish(R, 3)
  var_eY <- 1/prec_eY
  Psi_vY <- inverse(prec_vY)
  for (k in 1:2) {
    a0[k] ~ dnorm(0, 1e-06)
    aX[k] ~ dnorm(0, 1e-06)
    for (z in 1:numZ) {
      aZ[z, k] ~ dnorm(0, 1e-06)
    }
  }
  for (k in 1:2) {
    b0[k] ~ dnorm(0, 1e-06)
  }
  for (z in 1:numZ) {
    bZ[z, 1:2] ~ dmnorm(zero, prec_vY[1:2, 1:2])
  }
  bX[1:2] ~ dmnorm(zero, prec_vY[1:2, 1:2])
  bIM[1:2] ~ dmnorm(zero, prec_vY[1:2, 1:2])
  bSM[1:2] ~ dmnorm(zero, prec_vY[1:2, 1:2])
  bIMSM[1:2] ~ dmnorm(zero, prec_vY[1:2, 1:2])
  bXIM[1:2] ~ dmnorm(zero, prec_vY[1:2, 1:2])
  bXSM[1:2] ~ dmnorm(zero, prec_vY[1:2, 1:2])
  bXIMSM[1:2] ~ dmnorm(zero, prec_vY[1:2, 1:2])
  
  mIMt <- a0[1] + aX[1]+t(aZ[ ,1])%*%mZ
  mIMc <- a0[1] +t(aZ[ ,1])%*%mZ
  mSMt <- a0[2] + aX[2]+t(aZ[ ,2])%*%mZ
  mSMc <- a0[2] +t(aZ[ ,2])%*%mZ
  
  # the estimators of the interventional indirect and direct effects
  
  # IIE due to the mutual dependence on IY and on SY
  X_mu[1:2] <- (bIMSM[1:2]+bXIMSM[1:2]) * (Psi_vMtrt[1, 2] - Psi_vMctrl[1, 2])
  # pure IIE due to the mutual dependence on IY and on SY
  pX_mu[1:2] <- (bIMSM[1:2]) * (Psi_vMtrt[1, 2] - Psi_vMctrl[1, 2])
  # IIE via IM on IY and on SY
  X_IM[1:2] <- aX[1] * (bIM[1:2] +bXIM[1:2] + (bIMSM[1:2]+bXIMSM[1:2]) * mSMc[1] )
  # pure IIE via IM on IY and on SY
  pX_IM[1:2] <- aX[1] * (bIM[1:2] + (bIMSM[1:2]) * mSMc[1] )
  # IIE via SM on IY and on SY
  X_SM[1:2] <- aX[2] * (bSM[1:2] +bXSM[1:2] + (bIMSM[1:2]+bXIMSM[1:2]) * mIMt[1] )   # pure IIE via SM on IY and on SY
  pX_SM[1:2] <- aX[2] * (bSM[1:2] + (bIMSM[1:2]) * mIMt[1] )  
  # IIE joint on IY and on SY
  X_jo[1:2] <- X_mu[1:2] + X_IM[1:2] + X_SM[1:2]
  # pure IIE joint on IY and on SY
  pX_jo[1:2] <- pX_mu[1:2] + pX_IM[1:2] + pX_SM[1:2]
  
  # IDE on IY and on SY
  Xde[1:2] <- bX[1:2] +bXIM[1:2]*mIMc + bXSM[1:2]*mSMc + bXIMSM[1:2]*( mIMc[1]*mSMc[1] + t(aZ[ ,1])%*%varZ%*%aZ[ ,2] + Psi_vMctrl[1,2] )     
  # total IDE on IY and on SY
  tXde[1:2] <- bX[1:2] +bXIM[1:2]*mIMt + bXSM[1:2]*mSMt + bXIMSM[1:2]*( mIMt[1]*mSMt[1] + t(aZ[ ,1])%*%varZ%*%aZ[ ,2] + Psi_vMtrt[1,2] ) 
  
  # alternative IIE via IM on IY and on SY
  altX_IM[1:2] <- aX[1] * (bIM[1:2] +bXIM[1:2] + (bIMSM[1:2]+bXIMSM[1:2])  * mSMt[1] )
  # pure alternative IIE via IM on IY and on SY
  paltX_IM[1:2] <- aX[1] * (bIM[1:2] + (bIMSM[1:2])  * mSMt[1] )
  # alternative IIE via SM on IY and on SY
  altX_SM[1:2] <- aX[2] * (bSM[1:2] +bXSM[1:2] + (bIMSM[1:2]+bXIMSM[1:2]) * mIMc[1] )
  # pure alternative IIE via SM on IY and on SY
  paltX_SM[1:2] <- aX[2] * (bSM[1:2] + (bIMSM[1:2]) * mIMc[1] )
  
  # difference between two IIE versions 
  diffX_IM[1:2] <- altX_IM[1:2] - X_IM[1:2]
  diffX_SM[1:2] <- altX_SM[1:2] - X_SM[1:2]
  
}

# write the JAGS model to a file
write.model(jagsm, "jagsmod_interaction_XIMSM.txt")


# fitting the interaction model to the example data (https://github.com/xliu12/clgcmm/dat.RData)
load("dat.RData")
# observed mediator scores
M <- dat[ , grep("M",colnames(dat)) ]
# observed outcome scores
Y <- dat[ , grep("Y",colnames(dat)) ]
# treatment
X <- dat[, "X"]
# covariates
Z <- dat[ , grep("Z",colnames(dat)) ]

jagsdata <- list(
  M=M, Y=Y, X=X, 
  N = nrow(M),
  TT = ncol(M), 
  TimeM = ((1:ncol(M))-2), 
  TimeY = ((1:ncol(M))-ncol(M)),
  R=diag(1,nrow = 2),
  zero=c(0,0),
  Z=as.matrix(Z), numZ=1, 
  mZ=rep(0, 1), # column vector containing the mean of Z
  varZ=diag(1, nrow = 1), # matrix containing the variance-covariance matrix of Z
  whichX0=which(X==0), # row numbers of control group
  whichX1=which(X==1) # row numbers of treatment group
)

jagsout <- jags(data = jagsdata, jags.seed = 123,
      model.file = "jagsmod_interaction_XIMSM.txt", # the JAGS model script for fitting the interaction model (the default diffuse priors described in subsection "Bayesian estimation" are used)
      parameters.to.save = c("aX", "a0","aZ", "Psi_vMctrl","Psi_vMtrt","var_eM", # mediator model parameters 
      "bX","bIM","bSM","bIMSM","bXIM","bXSM","bXIMSM","b0","bZ", "Psi_vY", "var_eY", # outcome model parameters 
      "X_IM","altX_IM","X_SM", "altX_SM","X_mu", "Xde", 
      "pX_IM","paltX_IM","pX_SM", "paltX_SM","pX_mu", "tXde", # the estimators of the interventional (in)direct effects  
      "diffX_IM", "diffX_SM" # differences between the alternative and orignal versions of the interventional indirect effects via IM alone and those via SM alone  
      ), n.chains = 2,n.iter = 1e5, n.thin = 2
)

jagsout 

# if convergence is achieved, posterior means and 0.95 percential intervals can be saved as follows 
jagsres <- data.frame(jagsout$BUGSoutput$summary[ , c(1,3,7)])
round(jagsres, 3)


