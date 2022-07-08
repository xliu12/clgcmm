
library(parallel)
library(R2jags)
# library(MplusAutomation)
library(mvtnorm)

OneSimu <- function(seed=12,
                    N=200, # number of persons
                    TT=6, # number of time points nested within each person
                    rho12_vM1=0, # difference in the correlations of vIMi and vSMi between treatment vs. control groups
                    f2IMSM = c(0,0) # f^2 of the interaction IM*SM on IY SY
){
  set.seed(seed)
  
  ############ Data Generation ################ 
  #mediator intercept effects on IY SY 
  f2bIM= c(0.15,0.15);
  #mediator slope effects on IY SY 
  f2bSM= c(0.15,0.15); 
  #mediator intercept-slope interaction effects
  f2IMSM = c(0.15,0.15)
  #X on IY SY (direct effects) #Cohen's d
  dbX=c(0,0)
  
  TimeM = ((1:TT)-2)
  TimeY = ((1:TT)-TT)
  
  ## variance parameters
  vSM=16; 
  vIM=90; 
  vSY=9 #vIY=vSY*abs(TimeY[1])^2
  
  # sigma_eM=sqrt(50) ##when TT=3, GRR too small (around .3), over 1e5 iterations still nonconvergence ##when TT=6, 5e4 iters no convergence 
  if(TT==3) {sigma_eM=sqrt(50) } 
  if(TT>3) {sigma_eM=sqrt(50) } 
  
  if(TT==3){ sigma_eY=sqrt(16) }
  if(TT==6){ sigma_eY=sqrt(100) }
  if(TT==12){ sigma_eY=sqrt(400) }
  # pretreatment confounder 
  Z=rnorm(N)
  # binary randomized treatment 
  X=rbinom(N, 1, 0.5)
  
  # mediator 
  # structural models for IM SM 
  # with randomi effects of X
  # IMi = a_IM0 + (a_IM.X)*Xi + a_IM.Z*Zi + v_IMtrti*X + v_IMctrli*(1-X)
  # SMi = a_SM0 + (a_SM.X)*Xi + a_SM.Z*Zi + v_SMtrti*X + v_SMctrli*(1-X)
  
  a0=c(0,0)
  ##medium size f^2 for Z; a_IM.Z^2=f^2*vIM/((f^2+1)*Var(Z))
  aZ=rep(NA,2)
  aZ[1]=sqrt(0.02*vIM/((0.15+1)*1))
  aZ[2]= sqrt(0.02*vSM/((0.15+1)*1))
  ##medium size d=0.5
  aX=rep(NA,2)
  # aX[2]= 0.7*sqrt(vSM)
  ###aX set under the constraint: \alpha_{IM.X}+\alpha_{SM.X}Time_{M1}=0 
  aX[1]= 0.32*sqrt(vIM)
  aX[2]= -aX[1]/TimeM[1]
  # covariance matrix of the residuals of IM SM 
  # var(SM|X=0)*0.5+var(SM|X=1)*0.5 + Var(E(SM|X))
  # Psi_vMctrl = diag(c(vIM,vSM)-aZ^2*1, 2, 2)
  Psi_vMctrl = diag(c(vIM,vSM)-aX^2*0.25-aZ^2*1, 2, 2)
  Psi_vMtrt = Psi_vMctrl
  # rho12_vM1 = 0.3 # 0.1 0.3 0.5 small medium large
  rho12_vM1 =rho12_vM1
  Psi_vMtrt[1,2]=Psi_vMtrt[2,1]=rho12_vM1*sqrt(Psi_vMtrt[1,1]*Psi_vMtrt[2,2])
  
  vMctrl=rmvnorm(N, c(0,0), Psi_vMctrl)
  vMtrt=rmvnorm(N, c(0,0), Psi_vMtrt)
  
  IM=a0[1]+aX[1]*X+aZ[1]*Z + vMtrt[ ,1]*X + vMctrl[ ,1]*(1-X)
  SM=a0[2]+aX[2]*X+aZ[2]*Z + vMtrt[ ,2]*X + vMctrl[ ,2]*(1-X)
  
  # #check
  # ##f^2 for Z
  # (var(aX[1]*X+aZ[1]*Z)-var(aX[1]*X))/ (var(IM)-var(aX[1]*X+aZ[1]*Z))
  # (var(aX[2]*X+aZ[2]*Z)-var(aX[2]*X))/ (var(SM)-var(aX[2]*X+aZ[2]*Z))
  # ##d for X
  # (mean(IM[X==1])-mean(IM[X==0])) / sd(IM)
  # (mean(SM[X==1])-mean(SM[X==0])) / sd(SM)
  
  # outcome
  TimeY = ((1:TT)-TT)
  b0=c(0,0)
  vIY=vSY*abs(TimeY[1])^2
  ##### TT=12,6,3
  # if(TT==12){
  if( rho12_vM1==0 ){
    
    # f2IMSM=c(.15, .15) ##in this file where bIM=bSM=0, f2IMSM=c(.15, .15) means the nonzero indirect effects altX_IM, X_SM are entirely from the IM*SM interaction ; X_IM , altX_SM are 0, because a0=c(0,0), see below X_IM <- aX[1]*( bIM[1:2]+  (bIMSM[1:2] )*a0[2] ) ;
    if (f2IMSM[1]==0.15 & f2IMSM[2]==0.15){
      bZ=bIM=bSM=bX=bIMSM=rep(NA, 2)
      bIMSM[1]=sqrt( f2IMSM[1]*vIY/((0.15+1)*Psi_vMctrl[1,1]*Psi_vMctrl[2,2]))*0.6
      bIMSM[2]= -bIMSM[1]/TimeY[1]
      bZ[1]=sqrt(0.02*vIY/((0.15+1)*1))*0.1
      # bZ[2]= sqrt(0.02*vSY/((0.15+1)*1))*0.8
      bZ[2]=-bZ[1]/TimeY[1]
      #f2bIM=c(0,0)
      bIM[1]=sqrt(f2bIM[1] *vIY/((0.15+1)*Psi_vMctrl[1,1]))*0.8##diag(Psi_vMctrl)=diag(Psi_vMTrt)
      # bIM[2]= sqrt(0.15*vSY/((0.15+1)*Psi_vMctrl[1,1]))*0.9
      bIM[2]=-bIM[1]/TimeY[1]
      #f2bSM=c(0,0)
      bSM[1]=sqrt(f2bSM[1]*vIY/((0.15+1)*Psi_vMctrl[2,2]))*0.8 ##diag(Psi_vMctrl)=diag(Psi_vMTrt)
      # bSM[2]= sqrt(0.15*vSY/((0.15+1)*Psi_vMctrl[2,2]))*0.9
      bSM[2]=-bSM[1]/TimeY[1]
      bX[1]=dbX[1]*sqrt(vIY)*1
      bX[2]= -bX[1]/TimeY[1]
      
      # ##residual variance
      # varexplained=rep(NA,2)
      # ###calculate with N=100000
      # varexplained[1]=var(b0[1] + bX[1]*X + bIM[1]*IM + bSM[1]*SM + bIMSM[1]*IM*SM + bZ[1]*Z )
      # # varexplained[2]= var(b0[2] + bX[2]*X + bIM[2]*IM + bSM[2]*SM + bIMSM[2]*IM*SM + bZ[2]*Z)
      # varexplained[2]=varexplained[1]/(TimeY[1]^2)
      # Psi_vY= diag(c(vIY,vSY)-varexplained, 2, 2)
      
      
      Psi_vY=diag(c(NA, 5.7),2,2)
      Psi_vY[1,1]=(TimeY[1]^2)*Psi_vY[2,2]
    }
    
  }
  
  if( rho12_vM1==0.3 ){
    # f2IMSM=c(0,0)
    
    # f2IMSM=c(.15,.15)
    if (f2IMSM[1]==0.15 & f2IMSM[2]==0.15){
      bZ=bIM=bSM=bX=bIMSM=rep(NA, 2)
      bIMSM[1]=sqrt( f2IMSM[1]*vIY/((0.15+1)*Psi_vMctrl[1,1]*Psi_vMctrl[2,2]))*0.55
      bIMSM[2]= -bIMSM[1]/TimeY[1]
      bZ[1]=sqrt(0.02*vIY/((0.15+1)*1))*0.5
      # bZ[2]= sqrt(0.02*vSY/((0.15+1)*1))*0.8
      bZ[2]=-bZ[1]/TimeY[1]
      ##f2bIM= c(0,0);
      bIM[1]=sqrt(f2bIM[1]*vIY/((0.15+1)*Psi_vMctrl[1,1]))*0.7 ##diag(Psi_vMctrl)=diag(Psi_vMTrt)
      # bIM[2]= sqrt(0.15*vSY/((0.15+1)*Psi_vMctrl[1,1]))*0.9
      bIM[2]=-bIM[1]/TimeY[1]
      ##f2bSM= c(0,0); 
      bSM[1]=sqrt(f2bSM[1]*vIY/((0.15+1)*Psi_vMctrl[2,2]))*0.7 ##diag(Psi_vMctrl)=diag(Psi_vMTrt)
      # bSM[2]= sqrt(0.15*vSY/((0.15+1)*Psi_vMctrl[2,2]))*0.9
      bSM[2]=-bSM[1]/TimeY[1]
      ##in this file, dbX=c(.5,.5)
      bX[1]=dbX[1]*sqrt(vIY)*1.2
      bX[2]= -bX[1]/TimeY[1]
      
      # ##residual variance
      # varexplained=rep(NA,2)
      # ###calculate with N=100000
      # varexplained[1]=var(b0[1] + bX[1]*X + bIM[1]*IM + bSM[1]*SM + bIMSM[1]*IM*SM + bZ[1]*Z )
      # # varexplained[2]= var(b0[2] + bX[2]*X + bIM[2]*IM + bSM[2]*SM + bIMSM[2]*IM*SM + bZ[2]*Z)
      # varexplained[2]=varexplained[1]/(TimeY[1]^2)
      # Psi_vY= diag(c(vIY,vSY)-varexplained, 2, 2)
      
      
      Psi_vY=diag(c(NA, 5.8),2,2)
      Psi_vY[1,1]=(TimeY[1]^2)*Psi_vY[2,2]
    }
    
  }
  
  # structural models for IY SY
  # with interactions
  vY=rmvnorm(N, c(0,0), Psi_vY)
  
  IY= b0[1] + bX[1]*X + bIM[1]*IM + bSM[1]*SM + bIMSM[1]*IM*SM + bZ[1]*Z +vY[ ,1]
  SY= b0[2] + bX[2]*X + bIM[2]*IM + bSM[2]*SM + bIMSM[2]*IM*SM + bZ[2]*Z +vY[ ,2]
  
  # # #check with N=100000
  # ###IY
  # #f^2 for Z
  # (var(bX[1]*X +bZ[1]*Z +bIM[1]*IM + bSM[1]*SM +bIMSM[1]*IM*SM)-
  #     var(bX[1]*X +bIM[1]*IM + bSM[1]*SM +bIMSM[1]*IM*SM ) ) / Psi_vY[1,1]
  # #d for X
  # (mean(IY[X==1])-mean(IY[X==0])) / sd(IY)
  # #f^2 for IM
  # (var(bX[1]*X +bZ[1]*Z +bIM[1]*IM + bSM[1]*SM +bIMSM[1]*IM*SM)-
  #     var(bX[1]*X +bZ[1]*Z + bSM[1]*SM +bIMSM[1]*IM*SM) ) / Psi_vY[1,1]
  # #f^2 for SM
  # (var(bX[1]*X +bZ[1]*Z +bIM[1]*IM + bSM[1]*SM +bIMSM[1]*IM*SM)-
  #     var(bX[1]*X +bZ[1]*Z +bIM[1]*IM + bIMSM[1]*IM*SM) ) / Psi_vY[1,1]
  # #f^2 for IMSM
  # (var(bX[1]*X +bZ[1]*Z +bIM[1]*IM + bSM[1]*SM +bIMSM[1]*IM*SM)-
  #     var(bX[1]*X +bZ[1]*Z +bIM[1]*IM + bSM[1]*SM  ) ) / Psi_vY[1,1]
  # 
  # ###SY
  # # f^2 for Z
  # (var(bX[2]*X +bZ[2]*Z +bIM[2]*IM + bSM[2]*SM +bIMSM[2]*IM*SM)-
  #     var(bX[2]*X +bIM[2]*IM + bSM[2]*SM +bIMSM[2]*IM*SM) ) / Psi_vY[2,2]
  # #f^2 for X
  # (mean(SY[X==1])-mean(SY[X==0])) / sd(SY)
  # #f^2 for IM
  # (var(bX[2]*X +bZ[2]*Z +bIM[2]*IM + bSM[2]*SM +bIMSM[2]*IM*SM)-
  #     var(bX[2]*X +bZ[2]*Z + bSM[2]*SM +bIMSM[2]*IM*SM) ) / Psi_vY[2,2]
  # #f^2 for SM
  # (var(bX[2]*X +bZ[2]*Z +bIM[2]*IM + bSM[2]*SM +bIMSM[2]*IM*SM)-
  #     var(bX[2]*X +bZ[2]*Z +bIM[2]*IM + bIMSM[2]*IM*SM) ) / Psi_vY[2,2]
  # #f^2 for IMSM
  # (var(bX[2]*X +bZ[2]*Z +bIM[2]*IM + bSM[2]*SM +bIMSM[2]*IM*SM)-
  #     var(bX[2]*X +bZ[2]*Z +bIM[2]*IM + bSM[2]*SM  ) ) / Psi_vY[2,2]
  
  
  # LGC models
  # Mit= IMi + SMi*TimeM_it + e_Mit
  # Yit= IYi + SYi*TimeY_it + e_Yit
  TimeM = ((1:TT)-2) 
  TimeY = ((1:TT)-TT)
  
  M = matrix(NA, N, TT)
  Y = matrix(NA, N, TT)
  
  for(i in 1:N) {
    M[i, ]= IM[i] + TimeM*SM[i] + rnorm(TT,0,sigma_eM)
    Y[i, ]= IY[i] + TimeY*SY[i] + rnorm(TT,0,sigma_eY)
  }
  # check 
  # var(M[,3]) / var(M[,1])
  # datComp = cbind(Y, M, X, Z)
  # colnames(datComp) = c(paste0("Y",1:TT), paste0("M",1:TT), "X","Z")
  # datComp = as.data.frame(datComp)
  
  # dat=datComp
  # simudata=mget( ls() )
  # print( environment())
  
  
  ################# DATA ANALYSIS ##########################
  
  ### traditional LGCMM ####
  
  # jags data
  jagsdat=list(
    N=N,TT=TT, 
    TimeM=TimeM, TimeY=TimeY,
    M=M, Y=Y, X=X, Z=as.matrix(Z), numZ=1,
    R=diag(1,nrow = 2),
    zero=c(0,0)
  )
  
  # run jags
  if( TT==3 ){ niter=1.2e5 }
  if( TT==6 ){ niter=1e5 }
  if( TT==12 ){ niter=5e4 }
  jagsout_linear=jags(data = jagsdat, 
                      model.file = "jagsm_linear_conjugate.txt",
                      parameters.to.save = c("aX","aZ","Psi_vM","var_eM","a0",
                                             "b0","bX","bIM","bSM","bZ","Psi_vY","var_eY",
                                             'X_IM','X_SM','X_jo','Xde'),
                      n.chains = 2,n.iter = niter, n.thin = 2
  )
  jagsres_linear=data.frame(jagsout_linear$BUGSoutput$summary[ , c(1,2,3,7,8)])
  
  
  
  ##### Interaction LGCMM ##############
  jagsdat=list(
    N=N,TT=TT, 
    TimeM=TimeM, TimeY=TimeY,
    M=M, Y=Y, X=X, Z=as.matrix(Z), numZ=1,
    R=diag(1,nrow = 2),
    zero=c(0,0),
    whichX0=which(X==0),
    whichX1=which(X==1)
  )
  
  jagsout_Mint=jags(data = jagsdat, 
                    model.file = "jagsm_2way_binaryX_conjugate.txt",
                    parameters.to.save = c("Psi_vMtrt",'a0','aZ',"aX",'Psi_vMctrl','Mcovdiff',
                                           'Psi_vY',"b0","bX","bIMSM","bIM","bSM", 'var_eM','var_eY',
                                           'X_IM','X_SM','X_mu','X_jo','Xde','altX_IM','altX_SM'),
                    n.chains = 2,n.iter = niter, n.thin = 5
  )
  jagsres_Mint=data.frame(jagsout_Mint$BUGSoutput$summary[ , c(1,2,3,7,8)])
  
  
  ########### ONE Data RESULT ###################################
  
  #simu cond and true values
  
  X_mu<- (bIMSM[1:2] )*(Psi_vMtrt[1,2]-Psi_vMctrl[1,2])
  X_IM <- aX[1]*( bIM[1:2]+  (bIMSM[1:2] )*a0[2] ) ;
  X_SM <- aX[2]*(  bSM[1:2]+ (bIMSM[1:2] )*(a0[1]+aX[1]) ) ;
  X_jo = X_mu[1:2]+X_IM[1:2]+X_SM[1:2];
  Xde = bX[1:2]
  
  altX_IM <- aX[1]*( bIM[1:2]+  (bIMSM[1:2] )*(a0[2]+aX[2]) ) ;
  altX_SM <- aX[2]*(  bSM[1:2]+ (bIMSM[1:2] )*(a0[1] ) ) ;
  
  true_value=c(
    ##(in)direct effects
    X_IM=X_IM,X_SM=X_SM,X_mu=X_mu,X_jo=X_jo,Xde=Xde,altX_IM=altX_IM,altX_SM=altX_SM,
    ##key params that may be of interest
    aX=aX, Mcovdiff=Psi_vMtrt[1,2]-Psi_vMctrl[1,2],bX=bX,bIM=bIM,bSM=bSM,bIMSM=bIMSM, 
    # a0=a0,aZ=aZ,b0=b0,bZ=bZ,
    ##variance/cov components
    Psi_vMtrt=Psi_vMtrt[c(1,2,4)], Psi_vMctrl=Psi_vMctrl[c(1,2,4)], var_eM=sigma_eM^2, 
    Psi_vY=Psi_vY[c(1,2,4)],var_eY=sigma_eY^2
  )
  
  # out = data.frame(N=N, TT=TT,rho12_vM1=rho12_vM1, f2IMSM1=f2IMSM[1],f2IMSM2=f2IMSM[2], param=names(true_value), true_value=true_value) 
  out = data.frame(N=N, TT=TT,rho12_vM1=rho12_vM1, f2IMSM1=f2IMSM[1],f2IMSM2=f2IMSM[2], 
                   f2bIM1=f2bIM[1],f2bIM2=f2bIM[2],f2bSM1=f2bSM[1],f2bSM2=f2bSM[2], dbX1=dbX[1],dbX2=dbX[2],
                   param=names(true_value), true_value=true_value) 
  
  ### for jagsres_linear
  colnames(jagsres_linear)[1:5]=c('Estimate','PostSD','CIlow','CIup','PSR')
  jagsres_linear$param= gsub("\\]",'', gsub("\\[",'', row.names(jagsres_linear)))
  
  out1=merge(out, jagsres_linear, by='param', all.x = T)
  
  out1[out1$param %in% c('altX_IM1','altX_IM2','altX_SM1','altX_SM2'), c('Estimate','PostSD','CIlow','CIup','PSR')]=jagsres_linear[jagsres_linear$param%in%c('X_IM1','X_IM2','X_SM1','X_SM2'), 1:5]
  out1[out1$param %in% c('bIMSM1','bIMSM2','Mcovdiff','X_mu1','X_mu2'), c('Estimate','PostSD','CIlow','CIup','PSR')]=0
  out1[out1$param %in% paste0('Psi_vMtrt',1:3), c('Estimate','PostSD','CIlow','CIup','PSR')]=jagsres_linear[jagsres_linear$param%in%c("Psi_vM1,1", "Psi_vM2,1","Psi_vM2,2"), 1:5]
  out1[out1$param %in% paste0('Psi_vMctrl',1:3), c('Estimate','PostSD','CIlow','CIup','PSR')]=jagsres_linear[jagsres_linear$param%in%c("Psi_vM1,1", "Psi_vM2,1","Psi_vM2,2"), 1:5]
  out1[out1$param %in% paste0('Psi_vY',1:3), c('Estimate','PostSD','CIlow','CIup','PSR')]=jagsres_linear[jagsres_linear$param%in%c("Psi_vY1,1", "Psi_vY2,1","Psi_vY2,2"), 1:5]
  
  #DIC
  out1$DIC=jagsout_linear$BUGSoutput$DIC
  out1$method='linear'
  
  
  ###for jagsres_Mint
  colnames(jagsres_Mint)[1:5]=c('Estimate','PostSD','CIlow','CIup','PSR')
  jagsres_Mint$param= gsub("\\]",'', gsub("\\[",'', row.names(jagsres_Mint)))
  
  out2=merge(out, jagsres_Mint, by='param', all.x = T)
  out2[out2$param %in% paste0('Psi_vMtrt',1:3), c('Estimate','PostSD','CIlow','CIup','PSR')]=jagsres_Mint[jagsres_Mint$param%in%c("Psi_vMtrt1,1", "Psi_vMtrt2,1","Psi_vMtrt2,2"), 1:5]
  out2[out2$param %in% paste0('Psi_vMctrl',1:3), c('Estimate','PostSD','CIlow','CIup','PSR')]=jagsres_Mint[jagsres_Mint$param%in%c("Psi_vMctrl1,1", "Psi_vMctrl2,1","Psi_vMctrl2,2"), 1:5]
  out2[out2$param %in% paste0('Psi_vY',1:3), c('Estimate','PostSD','CIlow','CIup','PSR')]=jagsres_Mint[jagsres_Mint$param%in%c("Psi_vY1,1", "Psi_vY2,1","Psi_vY2,2"), 1:5]
  
  #DIC
  out2$DIC=jagsout_Mint$BUGSoutput$DIC
  out2$method='Mint'
  
  out12 = rbind(out1, out2)
  # outlist=list(out.onedata=out12, jagsout_linear=jagsout_linear, jagsout_Mint=jagsout_Mint)
  outlist=out12
  
  return(outlist)
}

