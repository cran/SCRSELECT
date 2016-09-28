#' Performs a grid search over the marginal posterior probabilities of inclusion and returns a list of DIC values corresponding to each grid point. This is used in the ReturnModel function.
#'
#' @importFrom stats dgamma dnorm dpois rgamma rnorm runif
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @param PCT1 Vector Containing posterior probabilities of inclusion for the hazard of a non-terminal event. This must be of length ncol(COV)-inc.
#' @param PCT2 Vector Containing posterior probabilities of inclusion for the hazard of death without a non-terminal event. This must be of length ncol(COV)-inc.
#' @param PCT3 Vector Containing posterior probabilities of inclusion for the hazard of death after a non-terminal event. This must be of length ncol(COV)-inc.
#' @param Y1   Vector Containing non-terminal event times (or censoring time due to death/censoring).
#' @param I1   Vector Containing non-terminal event indicators (1 if non-terminal event for a patient, 0 otherwise).
#' @param Y2   Vector Containing Terminal Event times (or censoring).
#' @param I2   Vector Containing Terminal event indicators (1 if a patients experiences a non-ternminal event, 0 if censored).
#' @param COV  Matrix of Patient Covariates. The last inc will be left out of variable selection.
#' @param s1 Vector containing the posterior locations of the split points in the hazard of a non-terminal event.
#' @param lam1 Vector containing the posterior log hazard heights on the split point intervals in the hazard of a non-terminal event.
#' @param s2 Vector containing the posterior locations of the split points in the hazard of death without a non-terminal event.
#' @param lam2 Vector containing the posterior log hazard heights on the split point intervals in the hazard of death without a non-terminal event.
#' @param s3 Vector containing the posterior locations of the split points in the hazard of death after a non-terminal event.
#' @param lam3 Vector containing the posterior log hazard heights on the split point intervals in the hazard of death after a non-terminal event.
#' @param gam Vector of length n containing the posterior mean frailties of the patients.
#' @param c Hyperparameter involved in the sampling of hazard coefficients. This should be the same value that controls the degree of sparsity achieved by the SVSS.
#' @param B Number of iterations
#' @param inc Number of variables left out of selection
#' @return Returns a list of size 18 containing 18x18 matrices of DIC values and skipped entries.
#'
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
#' ###Read in Posterior mean quantities from SCRSELECTRETURN
#' PCT1=c(.2,.4,.7,.8,.5)
#' PCT2=c(.02,.06,.1,.5,.7)
#' PCT3=c(.85,.87,.3,.45,.51)
#' gam=rgamma(n,1,1)
#' s1=c(0,3,5,max(Y1[I1==1]))
#' lam1=c(-1,-3,0)
#' s2=c(0,1,max(Y2[I1==0]))
#' lam2=c(0,-2)
#' s3=c(0,max(Y2[I1==1]))
#' lam3=-2
#' ####Read in Hyperparameters
#' c=5
#' ###Number of iterations and output location
#' B=4
#' ###Number of variables to exclude from selection and burnin percent
#'inc=2
#' DICTAUG(PCT1,PCT2,PCT3,X,Y1,Y2,I1,I2,s1,lam1,s2,lam2,s3,lam3,gam,c,B,inc)
#' @export
DICTAUG=function(PCT1,PCT2,PCT3,COV,Y1,Y2,I1,I2,s1,lam1,s2,lam2,s3,lam3,gam,c,B,inc){
  Inc=inc

  sum1=-1
  sum2=-1
  sum3=-1

  tau1=seq(.05,.9,by=.05)
  tau2=tau1
  tau3=tau1

  E1=rep(0,ncol(COV))
  E2=rep(0,ncol(COV))
  E3=rep(0,ncol(COV))


  Include=(ncol(COV)-Inc+1):ncol(COV)


  DICMAT=matrix(rep(NA,length(tau1)^2),nrow=length(tau1))

  DICLIST=list(rep(0,length(tau1)))

  for(g in 1:length(tau1)){

    E1=PCT1>tau1[g]

    cat("

        Grid Search on tau_1=", tau1[g], "Started

        ")

    if(sum(E1)==sum1){
      DICLIST[[g]]="skip"
    }else{

      sum1=sum(E1)
      for(h in 1:length(tau2)){

        E2=PCT2>tau2[h]

        if(sum(E2)==sum2){
          DICMAT[h,]=rep("skip",length(tau2))
        }else{

          sum2=sum(E2)
          for(l in 1:length(tau3)){



            J1=length(s1)-2
            J2=length(s2)-2
            J3=length(s3)-2
            G1=J1+1
            G2=J2+1
            G3=J3+1

            E3=PCT3>tau3[l]

            ###########################

            if( sum(E3)==sum3){

              DICMAT[h,l]="Skip"




            }else{



              sum2=sum(E2)
              sum3=sum(E3)

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





              ###Start DIC



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


                Like[b]=LK3(Y1,Y2,I1,I2,beta3[b,])+LK2(Y1,Y2,I1,I2,beta2[b,])+LK1(Y1,Y2,I1,I2,beta1[b,])

                ###End
              }








              mbeta1=colMeans(beta1[(B/2):B,])
              mbeta2=colMeans(beta2[(B/2):B,])
              mbeta3=colMeans(beta3[(B/2):B,])






              A=LK1(Y1,Y2,I1,I2,mbeta1)+
                LK2(Y1,Y2,I1,I2,mbeta2)+
                LK3(Y1,Y2,I1,I2,mbeta3)

              pdic=(-2*mean(Like[(B/2):B])+2*A)


              DIC=-2*A+2*pdic



              DICMAT[h,l]=DIC

            }

          }
        }
      }
      DICLIST[[g]]= DICMAT


    }




  }

  return(DICLIST)

}

