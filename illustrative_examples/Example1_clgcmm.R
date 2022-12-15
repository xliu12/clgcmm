# Using the Bayesian approach for
# estimating the interventional (in)direct effects in the LGCMM
# with the traditional (Trad) and interaction (Mint) models 

# load R packages needed
library(R2jags)


# Example 1 (dataset simulated under the traditional model Trad) ----

# read in the data
dat <- read.csv("data_simulated_under_Trad.csv")
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
ex1_jagsout_Trad <- jags(data = jagsdat, jags.seed = 123,
                    model.file = "jagsmod_Trad.txt", # the JAGS model script for fitting Trad (the default diffuse priors described in subsection "Bayesian estimation" are used)
                    parameters.to.save = c( "aX","a0","aZ","Psi_vM","var_eM", # mediator model parameters 
                                           "bX","bIM","bSM","b0","bZ","Psi_vY","var_eY", # outcome model parameters
                                           'X_IM','X_SM','Xde' # the estimators of the interventional (in)direct effects 
                                           ),
                    n.chains = 2,n.iter = 1e5, n.thin = 2
)
ex1_jagsout_Trad
# save posterior means and 0.95 percential intervals
ex1_jagsres_Trad <- data.frame(ex1_jagsout_Trad$BUGSoutput$summary[ , c(1,3,7)])
round(ex1_jagsres_Trad, 3)
# Note that with Trad, 
# the original and alternative versions of the IIEs via IM (or via SM) alone are constrained to be the same; 
# the IIEs due to the mutual dependence are constrained to be 0 (IIE_mu=0)

#                   mean     X2.5.    X97.5.
# X_IM[1]        -0.312    -0.400    -0.233  # the estimated IIE via IM on IY 
# X_SM[1]        -0.113    -0.290     0.063  # the estimated IIE via SM on IY 
# Xde[1]         -0.003    -0.216     0.211  # the estimated IDE on IY

# X_IM[2]         0.002    -0.035     0.039  # the estimated IIE via IM on SY 
# X_SM[2]         0.039    -0.060     0.139  # the estimated IIE via SM on SY
# Xde[2]         -0.079    -0.200     0.042  # the estimated IDE on SY
write.csv(ex1_jagsres_Trad, "Example_1_jags_results_Trad.csv")


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

ex1_jagsout_Mint <- jags(data = jagsdat, jags.seed = 123,
                     model.file = "jagsmod_Mint.txt", # the JAGS model script for fitting Mint (the default diffuse priors described in subsection "Bayesian estimation" are used)
                     parameters.to.save = c("aX", "a0","aZ", "Psi_vMctrl","Psi_vMtrt","var_eM", # mediator model parameters 
                                            "bX","bIM","bSM","bIMSM","b0","bZ", "Psi_vY", "var_eY", # outcome model parameters 
                                            "X_IM","altX_IM","X_SM", "altX_SM","X_mu", "Xde", # the estimators of the interventional (in)direct effects  
                                            "diffX_IM", "diffX_SM" # differences between the alternative and orignal versions of the interventional indirect effects via IM alone and those via SM alone  
                     ),
                     n.chains = 2,n.iter = 1e5, n.thin = 2
)
ex1_jagsout_Mint
# save posterior means and 0.95 percential intervals
ex1_jagsres_Mint<-data.frame(ex1_jagsout_Mint$BUGSoutput$summary[ , c(1,3,7)])
round(ex1_jagsres_Mint, 3)
# Note that with Mint, 
# the original and alternative versions of the IIEs via IM (or via SM) alone are allowed to be different; 
# the IIEs due to the mutual dependence are allowed to be nonzero 

#                   mean     X2.5.    X97.5.
# X_IM[1]            -0.338    -0.442    -0.245  # the estimated IIE via IM on IY 
# altX_IM[1]         -0.285    -0.384    -0.196  # the estimated alt.IIE via IM on IY 
# diffX_IM[1]         0.053    -0.042     0.152  # the estimated difference between the original (IIE) and alternative (alt.IIE) versions of the IIE via IM on IY 
# X_SM[1]            -0.095    -0.273     0.084  # the estimated IIE via SM on IY 
# altX_SM[1]         -0.148    -0.334     0.039  # the estimated alt.IIE via SM on IY 
# diffX_SM[1]        -0.053    -0.152     0.042  # the estimated difference between the original (IIE) and alternative (alt.IIE) versions of the IIE via SM on IY 
# X_mu[1]             0.004    -0.003     0.013  # the estimated IIE on IY due to the mutual dependent between IM on SM 
# Xde[1]              0.001    -0.212     0.213  # the estimated IDE on IY

# X_IM[2]             0.016    -0.029     0.061  # the estimated IIE via IM on SY 
# altX_IM[2]         -0.010    -0.057     0.037  # the estimated alt.IIE via IM on SY 
# diffX_IM[2]        -0.025    -0.081     0.029  # the estimated difference between the original (IIE) and alternative (alt.IIE) versions of the IIE via IM on SY 
# X_SM[2]             0.030    -0.071     0.131  # the estimated IIE via SM on SY 
# altX_SM[2]          0.056    -0.049     0.162  # the estimated alt.IIE via SM on SY 
# diffX_SM[2]         0.025    -0.029     0.081  # the estimated difference between the original (IIE) and alternative (alt.IIE) versions of the IIE via SM on SY 
# X_mu[2]            -0.002    -0.007     0.002  # the estimated IIE on SY due to the mutual dependent between IM on SM 
# Xde[2]             -0.082    -0.201     0.037  # the estimated IDE on SY

write.csv(ex1_jagsres_Mint, "Example_1_jags_results_Mint.csv")

save.image("Example1_clgcmm.RData")


