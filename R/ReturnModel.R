#' Performs the DIC-tau_g procedure and returns the posterior quantities of the optimal model.
#'
#' Performs the DIC-tau_g procedure by first running the function SCRSELECTRUN with 60% burnin, which performs SVSS on two disperse starting values for beta1,beta2,beta3. Afterwards, the function DICTAUG is used to extract the DIC values for unique models visited by the grid search and the optimal model is determined as the one with the lowest DIC which is the most parsimonius. After the optimal model is determined, one final MCMC is performed to obtain posterior beta1,beta2 and beta3 quantities for this model, returning summary values for each hazard.
#'
#' @importFrom graphics par plot
#' @importFrom stats dgamma dnorm dpois rgamma rnorm runif median quantile
#' @importFrom utils write.table
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @param Y1   Vector Containing non-terminal event times (or censoring time due to death/censoring)
#' @param I1   Vector Containing non-terminal event indicators (1 if non-terminal event for a patient, 0 otherwise)
#' @param Y2   Vector Containing Terminal Event times (or censoring)
#' @param I2   Vector Containing Terminal event indicators (1 if a patients experiences a non-ternminal event, 0 if censored)
#' @param X  Matrix of Patient Covariates. The last inc will be left out of variable selection.
#' @param hyperparameters  List containing 29 hyperparameters and four starting values. In order they are: psi-the swap rate of the SVSS algorithm.
#'  c-parameter involved in Sigma matrix for selection. z1a, z1b, z2a, z2b, z3a, z3b - beta hyper parameters on probability of inclusion for each of the three hazard functions.
#'  a1,b1,a2,b2,a3,b3- hyperparameters on sigma_lambda_1, sigma_lambda_2, and sigma_lambda_3.
#'   clam1, clam2, clam3 - spatial dependency of baseline hazard (between 0 and 1) for the three hazard functions.
#'    Alpha1, Alpha2, Alpha3 - The parameter for the number of split points in hazards 1,2 and 3 (must be whole number).
#'    J1max, J2max, J3max - Maximum number of split points allowed (must be whole number).
#'    J1, J2, J3- Starting number of split points. w, psi1- hyperparameters on theta^{-1}. cep=Tuning Parameter for theta^{-1} sampler.
#'    epstart-Starting value for theta^{-1}. cl1,cl2,cl3-Tuning parameters for log baseline hazard height sampler.
#' @param c sparsity parameter involved in Sigma matrix for selection. This should be the same c as that used in the hyperparameters vector.
#' @param BSVSS Number of iterations to perform during the SVSS procedure. 100,000 is a reccomended value to achieve convergence.
#' @param BDIC Number of iterations to perform during the DIC-tau_g grid search. 10,000 is a reccomended value to achieve convergence in a reasonable amount of time.
#' @param inc Number of variables left out of selection.
#' @param Path Where to save posterior coefficient samples for the optimal model.
#' @return Returns the optimal model determined by the DIC-Tau_g procedure and it's DIC along with summaries of these posterior quantities. Additionally, this function saves these posterior samples to a desired path.
#' @examples
#' ####Randomly Generate Semicompeting Risks Data
#' ####Generates random patient time, indicator and covariates.
#' n=100
#' Y1=runif(n,0,100)
#' I1=rbinom(n,1,.5)
#' Y2=Y1
#' I2=I1
#' for(i in 1:n){if(I1[i]==0){Y2[i]=Y1[i]}else{Y2[i]=Y1[i]+runif(1,0,100)}}
#' I2=rbinom(n,1,.5)
#' library(mvtnorm)
#' X=rmvnorm(n,rep(0,7),diag(7))
#' ####Read in Hyperparameters
#' ##Swap Rate
#' psi=.5
#' c=5
#' ###Eta Beta function probabilities
#' z1a=.4
#' z1b=1.6
#' z2a=.4
#' z2b=1.6
#' z3a=.4
#' z3b=1.6
#' ####Hierarchical lam params
#' ###Sigma^2 lambda_g hyperparameters
#' a1=.7
#' b1=.7
#' a2=a1
#' b2=b1
#' a3=a1
#' b3=b1
#' ##Spacing dependence c in [0,1]
#' clam1=1
#' clam2=1
#' clam3=1
#' #####NumSplit
#' alpha1=3
#' alpha2=3
#' alpha3=3
#' J1max=10
#' J2max=10
#' J3max=10
#' ####Split Point Starting Value ###
#' J1=3
#' J2=3
#' J3=3
#' ###epsilon starting values/hyperparameters###
#' w=.7
#' psi1=.7
#' cep=2.4
#' #############
#' epstart=1.5
#' cl1=.25
#' cl2=.25
#' cl3=.25
#' ###Beta Starting Values
#' hyper1=c(psi,c,z1a,z1b,z2a,z2b,z3a,z3b,a1,b1,a2,b2,a3,b3,clam1,clam2,clam3)
#' hyper2=c(alpha1,alpha2,alpha3,J1max,J2max,J3max,J1,J2,J3,w,psi1,cep,epstart,cl1,cl2,cl3)
#' hyper=c(hyper1,hyper2)
#' ###Number of iterations and output location
#' BSVSS=10
#' BDIC=4
#'Path=tempdir()
#' ###Number of variables to exclude from selection and burnin percent
#'inc=2
#' ReturnModel(Y1,I1,Y2,I2,X,hyper,inc,c,BSVSS,BDIC,Path)
#' @export
ReturnModel=function(Y1,I1,Y2,I2,X,hyperparameters,inc,c,BSVSS,BDIC,Path){
  cat("Running SVSS MCMC

      ")



  B=BSVSS
  burn=.6








  beta1start=rep(1,ncol(X))
  beta2start=beta1start
  beta3start=beta1start
  z1=SCRSELECTRUN(Y1,I1,Y2,I2,X,hyperparameters,beta1start,beta2start,beta3start,B,inc,burn)

  cat("

      Chain 1 Finished

      ")

  beta1start=c(rep(0,ncol(X)-inc),-1,-1)
  beta2start=beta1start
  beta3start=beta1start
  z2=SCRSELECTRUN(Y1,I1,Y2,I2,X,hyperparameters,beta1start,beta2start,beta3start,B,inc,burn)

  cat("

      Chain 2 Finished


      ")
  PCT1=colMeans(rbind(z1[[1]],z2[[1]]))
  PCT2=colMeans(rbind(z1[[2]],z2[[2]]))
  PCT3=colMeans(rbind(z1[[3]],z2[[3]]))
















  gam=colMeans(rbind(z1[[4]],z2[[4]]))


  J1=c(z1[[5]],z2[[5]])
  J2=c(z1[[6]],z2[[6]])
  J3=c(z1[[7]],z2[[7]])

  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }

  J1=Mode(J1)
  J2=Mode(J2)
  J3=Mode(J3)

  s1=rbind(z1[[8]],z2[[8]])
  s2=rbind(z1[[9]],z2[[9]])
  s3=rbind(z1[[10]],z2[[10]])

  lam1=rbind(z1[[11]],z2[[11]])
  lam2=rbind(z1[[12]],z2[[12]])
  lam3=rbind(z1[[13]],z2[[13]])


  s=matrix(rep(NA,length(s1)),nrow=nrow(s1))
  l=matrix(rep(NA,length(lam1)),nrow=nrow(lam1))

  for(b in 1:nrow(s)){
    if(is.na(s1[b,(J1+2)])==FALSE && is.na(s1[b,(J1+3)])==TRUE){
      s[b,1:(J1+2)]=s1[b,1:(J1+2)]
      l[b,1:(J1+1)]=lam1[b,1:(J1+1)]
    }
  }


  s1=colMeans(s,na.rm=TRUE)[1:(J1+2)]
  lam1=colMeans(l,na.rm=TRUE)[1:(J1+1)]


  s=matrix(rep(NA,length(s2)),nrow=nrow(s2))
  l=matrix(rep(NA,length(lam2)),nrow=nrow(lam2))

  for(b in 1:nrow(s)){
    if(is.na(s2[b,(J2+2)])==FALSE && is.na(s2[b,(J2+3)])==TRUE){
      s[b,1:(J2+2)]=s2[b,1:(J2+2)]
      l[b,1:(J2+1)]=lam2[b,1:(J2+1)]
    }
  }

  s2=colMeans(s,na.rm=TRUE)[1:(J2+2)]
  lam2=colMeans(l,na.rm=TRUE)[1:(J2+1)]


  s=matrix(rep(NA,length(s3)),nrow=nrow(s3))
  l=matrix(rep(NA,length(lam3)),nrow=nrow(lam3))

  for(b in 1:nrow(s)){
    if(is.na(s3[b,(J3+2)])==FALSE && is.na(s3[b,(J3+3)])==TRUE){
      s[b,1:(J3+2)]=s3[b,1:(J3+2)]
      l[b,1:(J3+1)]=lam3[b,1:(J3+1)]
    }
  }


  s3=colMeans(s,na.rm=TRUE)[1:(J3+2)]
  lam3=colMeans(l,na.rm=TRUE)[1:(J3+1)]


  Z=list(PCT1,PCT2,PCT3,gam,s1,s2,s3,lam1,lam2,lam3)


  cat("

      Marginal Posterior Probabilities of Inclusion

      ")


  PCT1=Z[[1]]
  PCT2=Z[[2]]
  PCT3=Z[[3]]




  cat("

      Hazard 1: Non-terminal Event

      ")

  print(PCT1)

  cat("

      Hazard 2: Death Before Non-terminal Event

      ")


  print(PCT2)

  cat("

      Hazard 3: Death After Non-terminal Event

      ")


  print(PCT3)



  gam=Z[[4]]
  s1=Z[[5]]
  s2=Z[[6]]
  s3=Z[[7]]
  lam1=Z[[8]]
  lam2=Z[[9]]
  lam3=Z[[10]]

  B=BDIC

  cat("

      Starting DIC-tau_g Procedure

      ")

  X1=DICTAUG(PCT1,PCT2,PCT3,X,Y1,Y2,I1,I2,s1,lam1,s2,lam2,s3,lam3,gam,c,B,inc)


  Inc=inc

  ext=function(X){

    m1=rep(NA,length(X))

    for(b in 1:length(X)){
      if(is.null(dim(X[[b]]))==FALSE){
        m1[b]=suppressWarnings(min(as.numeric(X[[b]]),na.rm=TRUE))
      }

    }

    return(m1)
  }


  Vec=ext(X1)

  VEC1=sort(unique(Vec))
  VEC1=VEC1[!is.na(VEC1)]

  ## How many different TAU1 candidates are there?
  ##Min DIC
  G1=min(Vec,na.rm=TRUE)
  t1=1
  for(b in 2:length(VEC1)){
    if(VEC1[b]<(G1+1)){
      t1=t1+1
    }
  }

  tau1=rep(NA,t1)
  tau12=rep(NA,t1)
  tau13=rep(NA,t1)

  for(k in 1:t1){
    for(b in 1:length(Vec)){
      if(is.na(Vec[b])==FALSE){
        if(VEC1[k]==Vec[b]){
          tau1[k]=b
        }
      }
    }
  }

  minDIC=rep(NA,t1)

  for( k in 1:t1){
    G2=as.numeric(min(X1[[tau1[k]]],na.rm=TRUE))
    NumCand=sum(X1[[tau1[k]]]<(G2+1))

    H=suppressWarnings(as.numeric(sort(unique(X1[[tau1[k]]]))))
    H=H[!is.na(H)]
    H=H[1:NumCand]
    j1=1
    tau2=rep(NA,NumCand)
    tau3=rep(NA,NumCand)
    for(j in 1:NumCand){
      for(l in 1:18){
        for(m in 1:18){
          if(X1[[tau1[k]]][l,m]==H[j]){
            tau2[j]=l/20
            tau3[j]=m/20
          }
        }
      }
    }

    if(NumCand==1){
      tau12[k]=tau2[1]
      tau13[k]=tau3[1]

      minDIC[k]=min(H[1])


    }else{
      sumINC=rep(NA,NumCand)

      for(j in 1:NumCand){
        sumINC[j]=sum(PCT2>tau2[j])+sum(PCT3>tau3[j])
      }

      minINC = min(sumINC)
      H1=sumINC==minINC

      spot1=rep(NA,length(sumINC))

      for(j in 1:NumCand){
        if(H1[j]==TRUE){
          spot1[j]=j
        }
      }

      spot1=spot1[!is.na(spot1)]


      minDIC[k]=min(H[spot1])

      for(j in 1:length(spot1)){
        if(H[spot1[j]]==minDIC[k]){
          spot2=spot1[j]
        }

      }



      tau12[k]=tau2[spot2]
      tau13[k]=tau3[spot2]





    }










  }

  taumark=tau1

  tau1=tau1/20



  sumINC=rep(NA,length(tau1))

  for(j in 1:t1){
    sumINC[j]=sum(PCT1>tau1[j])+sum(PCT2>tau12[j])+sum(PCT3>tau13[j])
  }

  minINC = min(sumINC)
  H1=sumINC==minINC




  spot1=rep(NA,length(sumINC))

  for(j in 1:t1){
    if(H1[j]==TRUE){
      spot1[j]=j
    }
  }

  spot1=spot1[!is.na(spot1)]


  minDIC1=min(minDIC)

  for(j in 1:length(spot1)){
    if(Vec[taumark[spot1[j]]]==minDIC1){
      spot2=spot1[j]
    }

  }



  tauG=c(tau1[spot2],tau12[spot2],tau13[spot2])


  cat("Grid Search Finished

      ")

  cat("DIC of best model: ",minDIC)
  cat("

      Tau_g value for optimum model

      ")

  cat(tauG)

  cat("

      Variables Included

      ")





  cat("

      Hazard 1: Non-terminal Event

      ")







  E1=PCT1>tauG[1]
  print(E1)
  cat("

      Hazard 2: Death Before Non-terminal Event

      ")

  E2=PCT2>tauG[2]
  print(E2)
  cat("

      Hazard 3: Death After Non-terminal Event

      ")

  E3=PCT3>tauG[3]
  print(E3)


  cat("

      Performing MCMC for Final Model

      ")

  COV=X


  eta1=rep(0,length(E1))
  eta2=rep(0,length(E2))
  eta3=rep(0,length(E3))


  for(i in 1:(ncol(COV)-inc)){
    if(E1[i]==1){
      eta1[i]=i
    }
  }


  for(i in 1:(ncol(COV)-inc)){
    if(E2[i]==1){
      eta2[i]=i
    }
  }

  for(i in 1:(ncol(COV)-inc)){
    if(E3[i]==1){
      eta3[i]=i
    }
  }


  eta1=eta1[!(eta1==0)]
  eta2=eta2[!(eta2==0)]
  eta3=eta3[!(eta3==0)]


  B=BSVSS




  Include=(ncol(COV)-Inc+1):ncol(COV)

  ###Covariate Matrices
  COV1=as.matrix(COV[,c(eta1,Include)])
  COV2=as.matrix(COV[,c(eta2,Include)])
  COV3=as.matrix(COV[,c(eta3,Include)])


  p1=ncol(COV1)
  p2=ncol(COV2)
  p3=ncol(COV3)




  ### Sets up Storage Matrices ##
  ##Hundred Thousand

  ###Beta/Eta###
  beta1=matrix(rep(0,B*(p1)),nrow=B)
  beta2=matrix(rep(0,B*(p2)),nrow=B)
  beta3=matrix(rep(0,B*(p3)),nrow=B)


  n=length(Y1)
  Like=rep(0,B)



  Indcond1=matrix(rep(0,p1*B),nrow=B)
  Indcond2=matrix(rep(0,p2*B),nrow=B)
  Indcond3=matrix(rep(0,p3*B),nrow=B)


  Sigma1=c*solve(t(COV1)%*%COV1)
  Sigma2=c*solve(t(COV2)%*%COV2)
  Sigma3=c*solve(t(COV3)%*%COV3)

  Indmix1=rep(0,B)
  Indmix2=rep(0,B)
  Indmix3=rep(0,B)




  G1=length(s1)-1
  G2=length(s2)-1
  G3=length(s3)-1



  LK1=function(Y1,Y2,I1,I2,Beta1){

    LOGBH=0
    et1=COV1%*%Beta1


    LOGBH=LOGBH+sum(I1*et1)

    for(k in 1:G1){


      Del=pmax(0,pmin(Y1,s1[k+1])-s1[k])




      LOGBH=LOGBH-sum(gam*Del*exp(lam1[k])*exp(et1))

    }



    return(LOGBH)
  }


  ##






  LK2=function(Y1,Y2,I1,I2,Beta2){

    LOGBH=0
    et1=COV2%*%Beta2
    LOGBH=LOGBH+sum(I2*(1-I1)*et1)

    for(k in 1:G2){


      Del=pmax(0,pmin(Y1,s2[k+1])-s2[k])




      LOGBH=LOGBH-sum(gam*Del*exp(lam2[k])*exp(et1))

    }

    return(LOGBH)

  }

  ###


  ##

  LK3=function(Y1,Y2,I1,I2,Beta3){

    LOGBH=0
    et1=COV3%*%Beta3

    LOGBH=LOGBH+sum(I2*(I1)*et1)

    for(k in 1:G3){


      Del=pmax(0,pmin(Y2-Y1,s3[k+1])-s3[k])




      LOGBH=LOGBH-sum(gam*Del*exp(lam3[k])*exp(et1))

    }

    return(LOGBH)

  }




  iter=0

  for(b in 2:B){

    iter="haz1"


    ##Print iteration
    if(b%%10000==0){cat(b, "iterations ")}

    beta1[b,]=beta1[b-1,]

    for(m in 1:p1){


      V1 = Sigma1[m,m]
      V2 = as.matrix(Sigma1[-m,-m])
      V12 = as.matrix(Sigma1[m,-m])
      thetab=beta1[b,]
      thetano = as.matrix(thetab[-m])
      meannew = t(V12)%*%solve(V2)%*%thetano
      varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
      ##################

      beta1[b,m]=rnorm(1,meannew,varnew)

      #beta1[b,m]=beta1[b-1,m] + runif(1,-clb,clb)
      dn=log(dnorm(beta1[b,m],meannew,varnew))
      ###density old
      do=log(dnorm(thetab[m],meannew,varnew))



      Likeo=LK1(Y1,Y2,I1,I2,thetab)

      Liken=LK1(Y1,Y2,I1,I2,beta1[b,])

      alpha=Liken-Likeo+dn-do
      U=log(runif(1,0,1))

      if(U>alpha){
        Indcond1[b,m]=0
        beta1[b,]=thetab
      }else{
        Indcond1[b,m]=1
      }

    }

    iter="haz2"




    beta2[b,]=beta2[b-1,]

    for(m in 1:p2){


      V1 = Sigma2[m,m]
      V2 = as.matrix(Sigma2[-m,-m])
      V12 = as.matrix(Sigma2[m,-m])
      thetab=beta2[b,]
      thetano = as.matrix(thetab[-m])
      meannew = t(V12)%*%solve(V2)%*%thetano
      varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
      ##################

      beta2[b,m]=rnorm(1,meannew,varnew)
      dn=log(dnorm(beta2[b,m],meannew,varnew))
      ###density old
      do=log(dnorm(thetab[m],meannew,varnew))



      Likeo=LK2(Y1,Y2,I1,I2,thetab)

      Liken=LK2(Y1,Y2,I1,I2,beta2[b,])

      alpha=Liken-Likeo+dn-do
      U=log(runif(1,0,1))

      if(U>alpha){
        Indcond2[b,m]=0
        beta2[b,]=thetab
      }else{
        Indcond2[b,m]=1
      }

    }





    iter="haz3"


    ##Print iteration


    beta3[b,]=beta3[b-1,]

    for(m in 1:p3){


      V1 = Sigma3[m,m]
      V2 = as.matrix(Sigma3[-m,-m])
      V12 = as.matrix(Sigma3[m,-m])
      thetab=beta3[b,]
      thetano = as.matrix(thetab[-m])
      meannew = t(V12)%*%solve(V2)%*%thetano
      varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
      ##################

      beta3[b,m]=rnorm(1,meannew,varnew)
      dn=log(dnorm(beta3[b,m],meannew,varnew))
      ###density old
      do=log(dnorm(thetab[m],meannew,varnew))



      Likeo=LK3(Y1,Y2,I1,I2,thetab)

      Liken=LK3(Y1,Y2,I1,I2,beta3[b,])

      alpha=Liken-Likeo+dn-do
      U=log(runif(1,0,1))

      if(U>alpha){
        Indcond3[b,m]=0
        beta3[b,]=thetab
      }else{
        Indcond3[b,m]=1
      }

    }



    ###End
  }


  cat("MCMC is finished, returning posterior quantities

      ")

  ##burnin

  beta1=beta1[(B*burn):B,]
  beta2=beta2[(B*burn):B,]
  beta3=beta3[(B*burn):B,]


  ##Re-parse out vector and return 0 for non-included variables.

  BETA1=matrix(rep(0,nrow(beta1)*ncol(COV)),nrow=nrow(beta1))
  BETA2=matrix(rep(0,nrow(beta1)*ncol(COV)),nrow=nrow(beta1))
  BETA3=matrix(rep(0,nrow(beta1)*ncol(COV)),nrow=nrow(beta1))



  BETA1[,c(eta1,Include)]=beta1
  BETA2[,c(eta2,Include)]=beta2
  BETA3[,c(eta3,Include)]=beta3



  Path1= paste0(Path,"/beta1.txt")

  write.table(BETA1, Path1, sep="\t")


  Path1= paste0(Path,"/beta2.txt")

  write.table(BETA2, Path1, sep="\t")



  Path1= paste0(Path,"/beta3.txt")

  write.table(BETA3, Path1, sep="\t")



  haz1=colMeans(beta1>0)
  haz2=colMeans(beta2>0)
  haz3=colMeans(beta3>0)

  HAZ1=rep(NA,ncol(COV))
  HAZ2=HAZ1
  HAZ3=HAZ1


  HAZ1[c(eta1,Include)]=haz1
  HAZ2[c(eta2,Include)]=haz2
  HAZ3[c(eta3,Include)]=haz3




  cat("

      Posterior Probability of Increasing Hazard

      ")




  cat("Hazard 1: Non-terminal Event

      ")

  print(HAZ1)

  cat("

      Hazard 2: Death Before Non-terminal Event

      ")


  print(HAZ2)

  cat("

      Hazard 3: Death After Non-terminal Event

      ")


  print(HAZ3)



  cat("

      Posterior Median Hazard Ratio and Credible Interval

      ")


  med1=exp(apply(beta1,2,median))

  med2=exp(apply(beta2,2,median))
  med3=exp(apply(beta3,2,median))

  ql1=exp(apply(beta1,2,quantile,probs=.025))
  ql2=exp(apply(beta2,2,quantile,probs=.025))
  ql3=exp(apply(beta3,2,quantile,probs=.025))


  qu1=exp(apply(beta1,2,quantile,probs=.975))
  qu2=exp(apply(beta2,2,quantile,probs=.975))
  qu3=exp(apply(beta3,2,quantile,probs=.975))


  MED1=rep(NA,ncol(COV))
  MED2=MED1
  MED3=MED1
  QL1=MED1
  QL2=MED1
  QL3=MED1
  QU1=MED1
  QU2=MED1
  QU3=MED1

  MED1[c(eta1,Include)]=med1
  MED2[c(eta2,Include)]=med2
  MED3[c(eta3,Include)]=med3


  QL1[c(eta1,Include)]=ql1
  QL2[c(eta2,Include)]=ql2
  QL3[c(eta3,Include)]=ql3


  QU1[c(eta1,Include)]=qu1
  QU2[c(eta2,Include)]=qu2
  QU3[c(eta3,Include)]=qu3







  cat("

      Hazard 1: Non-terminal Event

      ")

  print(MED1)
  print(QL1)
  print(QU1)

  cat("

      Hazard 2: Death Before Non-terminal Event

      ")


  print(MED2)
  print(QL2)
  print(QU2)


  cat("

      Hazard 3: Death After Non-terminal Event

      ")


  print(MED3)
  print(QL3)
  print(QU3)










}
