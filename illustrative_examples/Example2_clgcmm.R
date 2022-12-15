# Using the Bayesian approach for
# estimating the interventional (in)direct effects in the LGCMM
# with the traditional (Trad) and interaction (Mint) models 

# load R packages needed
library(R2jags)

# Example 2 (dataset simulated under the interaction model Mint) ----
# read in the data
dat <- read.csv("data_simulated_under_Mint.csv")
# observed mediator scores
M <- dat[ , grep("M",colnames(dat)) ]
# observed outcome scores
Y <- dat[ , grep("Y",colnames(dat)) ]
# treatment
X <- dat[, "X"]
# covariates
Z <- dat[ , grep("Z",colnames(dat)) ]
# sample size 
N <- nrow(dat)
# number of time points
TT <- ncol(M)
# specify slope loadings
TimeM <- ((1:TT)-2) 
TimeY <- ((1:TT)-TT)



### The traditional model (Trad): Bayesian estimation via JAGS ####
# data for fitting Trad
jagsdat <- list(
  N=N,TT=TT, 
  TimeM=TimeM, TimeY=TimeY,
  M=M, Y=Y, X=X, 
  Z=as.matrix(Z), numZ=1, # numZ: number of covariates (i.e., number of columns of the covariate matrix Z)
  R=diag(1,nrow = 2), # parameter for the inverse-Wishart prior 
  zero=c(0,0) # parameter for the normal prior
)

# run jags
ex2_jagsout_Trad <- jags(data = jagsdat, jags.seed = 123,
                         model.file = "jagsmod_Trad.txt", # the JAGS model script for fitting Trad (the default diffuse priors described in subsection "Bayesian estimation" are used)
                         parameters.to.save = c( "aX","a0","aZ","Psi_vM","var_eM", # mediator model parameters 
                                                 "bX","bIM","bSM","b0","bZ","Psi_vY","var_eY", # outcome model parameters
                                                 'X_IM','X_SM','Xde' # the estimators of the interventional (in)direct effects 
                         ),
                         n.chains = 2,n.iter = 1e5, n.thin = 2
)
ex2_jagsout_Trad
# save posterior means and 0.95 percential intervals
ex2_jagsres_Trad <- data.frame(ex2_jagsout_Trad$BUGSoutput$summary[ , c(1,3,7)])
round(ex2_jagsres_Trad, 3)
# Note that with Trad, 
# the original and alternative versions of the IIEs via IM (or via SM) alone are constrained to be the same; 
# the IIEs due to the mutual dependence are constrained to be 0 (IIE_mu=0)

#                   mean     X2.5.    X97.5.
# X_IM[1]        -0.553    -0.676    -0.440  # the estimated IIE via IM on IY 
# X_SM[1]        -0.177    -0.374     0.020  # the estimated IIE via SM on IY 
# Xde[1]         -0.020    -0.239     0.200  # the estimated IDE on IY

# X_IM[2]        -0.059    -0.105    -0.015  # the estimated IIE via IM on SY 
# X_SM[2]        -0.228    -0.339    -0.118 # the estimated IIE via SM on SY
# Xde[2]          0.084    -0.037     0.206 # the estimated IDE on SY
write.csv(ex2_jagsres_Trad, "Example_2_jags_results_Trad.csv")


##### The interaction model (Mint): Bayesian estimation via JAGS  ######
# data for fitting Mint
jagsdat <- list(
  N=N,TT=TT, 
  TimeM=TimeM, TimeY=TimeY,
  M=M, Y=Y, X=X, Z=as.matrix(Z), numZ=1,
  R=diag(1,nrow = 2),
  zero=c(0,0),
  whichX0=which(X==0), # row numbers of control group
  whichX1=which(X==1) # row numbers of treatment group
)

ex2_jagsout_Mint <- jags(data = jagsdat, jags.seed = 123,
                         model.file = "jagsmod_Mint.txt", # the JAGS model script for fitting Mint (the default diffuse priors described in subsection "Bayesian estimation" are used)
                         parameters.to.save = c("aX", "a0","aZ", "Psi_vMctrl","Psi_vMtrt","var_eM", # mediator model parameters 
                                                "bX","bIM","bSM","bIMSM","b0","bZ", "Psi_vY", "var_eY", # outcome model parameters 
                                                "X_IM","altX_IM","X_SM", "altX_SM","X_mu", "Xde", # the estimators of the interventional (in)direct effects  
                                                "diffX_IM", "diffX_SM" # differences between the alternative and orignal versions of the interventional indirect effects via IM alone and those via SM alone  
                         ),
                         n.chains = 2,n.iter = 1e5, n.thin = 2
)
ex2_jagsout_Mint
# save posterior means and 0.95 percential intervals
ex2_jagsres_Mint<-data.frame(ex2_jagsout_Mint$BUGSoutput$summary[ , c(1,3,7)])
round(ex2_jagsres_Mint, 3)
# Note that with Mint, 
# the original and alternative versions of the IIEs via IM (or via SM) alone are allowed to be different; 
# the IIEs due to the mutual dependence are allowed to be nonzero 

#                   mean     X2.5.    X97.5.
# X_IM[1]           -0.379    -0.490    -0.278  # the estimated IIE via IM on IY 
# altX_IM[1]        -0.750    -0.909    -0.604  # the estimated alt.IIE via IM on IY 
# diffX_IM[1]       -0.371    -0.493    -0.262  # the estimated difference between the original (IIE) and alternative (alt.IIE) versions of the IIE via IM on IY 
# X_SM[1]           -0.211    -0.408    -0.016  # the estimated IIE via SM on IY 
# altX_SM[1]         0.160    -0.047     0.370  # the estimated alt.IIE via SM on IY 
# diffX_SM[1]        0.371     0.262     0.493  # the estimated difference between the original (IIE) and alternative (alt.IIE) versions of the IIE via SM on IY 
# X_mu[1]            -0.066    -0.096    -0.041  # the estimated IIE on IY due to the mutual dependent between IM on SM 
# Xde[1]             -0.090    -0.302     0.121  # the estimated IDE on IY

# X_IM[2]             0.033    -0.016     0.084 # the estimated IIE via IM on SY 
# altX_IM[2]         -0.164    -0.225    -0.109  # the estimated alt.IIE via IM on SY 
# diffX_IM[2]        -0.025    -0.081     0.029  # the estimated difference between the original (IIE) and alternative (alt.IIE) versions of the IIE via IM on SY 
# X_SM[2]            -0.246    -0.356    -0.137 # the estimated IIE via SM on SY 
# altX_SM[2]         -0.048    -0.165     0.070 # the estimated alt.IIE via SM on SY 
# diffX_SM[2]         0.025    -0.029     0.081  # the estimated difference between the original (IIE) and alternative (alt.IIE) versions of the IIE via SM on SY 
# X_mu[2]            -0.035    -0.052    -0.021  # the estimated IIE on SY due to the mutual dependent between IM on SM 
# Xde[2]              0.047    -0.071     0.165  # the estimated IDE on SY

write.csv(ex2_jagsres_Mint, "Example_2_jags_results_Mint.csv")

save.image("Example2_clgcmm.RData")


