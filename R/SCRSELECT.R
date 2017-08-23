#' Performs Bayesian Variable Selection on the covariates in a semi-competing risks model
#' @importFrom graphics par plot
#' @importFrom stats dgamma dnorm dpois rgamma rnorm runif
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
#' @param beta1start Starting Values for Beta1
#' @param beta2start Starting Values for Beta2
#' @param beta3start Starting Values for Beta3
#' @param B Number of iterations
#' @param inc Number of variables left out of selection
#' @param Path Where to save posterior samples
#' @param burn  percent of posterior sample to burn in (burn*B must be a whole number)
#' @return Returns marginal posterior probability of inclusion (post burn-in) for each hazard function along with acceptance rates for the various Metropolis-Hastings (and Metropolis-Hastings-Green) samplers.
#'
#' @references
#' [1] Lee, K. H., Haneuse, S., Schrag, D. and Dominici, F. (2015), Bayesian semi-parametric analysis of semi-competing risks data: investigating hospital readmission after a pancreatic cancer diagnosis. Journal of the Royal Statistical Society: Series C (Applied Statistics), 64: 253-273. doi: 10.1111/rssc.12078
#' [2] Chapple, A.C., Vannucci, M., Thall, P.F., Lin, S.(2017), Bayesian Variable selection for a semi-competing risks model with three hazard functions. Journal of Computational Statistics & Data Analysis, Volume 112, August 2017, Pages 170-185
#' [3] https://adventuresinstatistics.wordpress.com/2017/04/10/package-scrselect-using-returnmodel/

#'
#' @examples
#' ####Randomly Generate Semicompeting Risks Data
#' set.seed(1)
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
#' beta1start=c(1,1,1,1,1,-1,-1)
#' beta2start=c(1,1,1,1,1,-1,-1)
#' beta3start=c(-1,1,1,1,1,-1,-1)
#' hyper1=c(psi,c,z1a,z1b,z2a,z2b,z3a,z3b,a1,b1,a2,b2,a3,b3,clam1,clam2,clam3)
#' hyper2=c(alpha1,alpha2,alpha3,J1max,J2max,J3max,J1,J2,J3,w,psi1,cep,epstart,cl1,cl2,cl3)
#' hyper=c(hyper1,hyper2)
#' ###Number of iterations and output location
#' B=100
#'Path=tempdir()
#' ###Number of variables to exclude from selection and burnin percent
#'inc=2
#'burn=.1
#'SCRSELECT(Y1,I1,Y2,I2,X,hyper,beta1start,beta2start,beta3start,B,inc,Path,burn)
#' @export
SCRSELECT=function(Y1,I1,Y2,I2,X,hyperparameters,beta1start,beta2start,beta3start,B,inc,Path,burn){


  Inc=inc
  inc=inc

  iter=c(0,0)

  if(inc%%1>0){
    cat("inc must be a natural number")
  }else{




    ####Hyperparameters##
    ##Swap Rate
    psi=hyperparameters[1]
    ##
    c=hyperparameters[2]
    ###Eta Beta function probabilities
    z1a=hyperparameters[3]
    z1b=hyperparameters[4]
    z2a=hyperparameters[5]
    z2b=hyperparameters[6]
    z3a=hyperparameters[7]
    z3b=hyperparameters[8]
    ####Hierarchical lam params
    ###Siglam
    a1=hyperparameters[9]
    b1=hyperparameters[10]
    a2=hyperparameters[11]
    b2=hyperparameters[12]
    a3=hyperparameters[13]
    b3=hyperparameters[14]
    ##Spacing dependence c in [0,1]
    clam1=hyperparameters[15]
    clam2=hyperparameters[16]
    clam3=hyperparameters[17]
    ##Lamsampler params
    #####NumSplit
    alpha1=hyperparameters[18]
    alpha2=hyperparameters[19]
    alpha3=hyperparameters[20]
    J1max=hyperparameters[21]
    J2max=hyperparameters[22]
    J3max=hyperparameters[23]
    ####Split Points###
    J1=hyperparameters[24]
    J2=hyperparameters[25]
    J3=hyperparameters[26]
    ###epsilon functions###
    ###hyperparams###
    w=hyperparameters[27]
    psi1=hyperparameters[28]
    cep=hyperparameters[29]

    epstart=hyperparameters[30]

    cl1=hyperparameters[31]

    cl2=hyperparameters[32]
    cl3=hyperparameters[33]



    p1=ncol(X)-inc

    n=length(Y1)











    #####In program
    ###Make Acceptance Matrices
    ###Beta/Eta###
    beta1=matrix(rep(1,B*(p1+inc)),nrow=B)
    beta2=matrix(rep(1,B*(p1+inc)),nrow=B)
    beta3=matrix(rep(1,B*(p1+inc)),nrow=B)
    eta1=matrix(rep(1,B*p1),nrow=B)
    eta2=matrix(rep(1,B*p1),nrow=B)
    eta3=matrix(rep(1,B*p1),nrow=B)
    ####Frailty Matrix###
    ###
    Mulam1=rep(0,B)
    Siglam1=rep(1,B)
    Mulam2=rep(0,B)
    Siglam2=rep(1,B)
    Mulam3=rep(0,B)
    Siglam3=rep(1,B)
    gam = matrix(rep(1,n*B),nrow=B)
    theta=rep(1,B)
    epsilon=rep(epstart,B)
    Indepsilon=rep(0,B)
    gam[1,]=rgamma(n,1/epsilon[1],1/epsilon[1])





    ###Make Eta1Start
    beta1[1,]=beta1start
    beta2[1,]=beta2start
    beta3[1,]=beta3start
    ##
    eta1start=rep(1,p1)
    eta2start=eta1start
    eta3start=eta1start

    for(i in 1:p1){
      if(beta1start[i]==0){
        eta1start[i]=0
      }

      if(beta2start[i]==0){
        eta2start[i]=0
      }

      if(beta3start[i]==0){
        eta3start[i]=0
      }
    }



    eta1[1,]=eta1start
    eta2[1,]=eta2start
    eta3[1,]=eta3start


    m1 = max(Y1[I1==1])+.001
    m2 = max(Y2[I1==0 & I2==1])+.001
    m3 = max(Y2[I1==1 & I2==1])+.001




    ####Acceptance Matrices


    Acceptlam1=matrix(rep(NA,B*(J1max+1)),nrow=B)
    Acceptlam2=matrix(rep(NA,B*(J2max+1)),nrow=B)
    Acceptlam3=matrix(rep(NA,B*(J3max+1)),nrow=B)
    accepts1=rep(0,B)
    accepts2=rep(0,B)
    accepts3=rep(0,B)
    Indmix1=rep(0,B)
    Indmix2=rep(0,B)
    Indmix3=rep(0,B)
    sum1=rep(0,B)
    sum2=rep(0,B)
    sum3=rep(0,B)
    split1=rep(0,B)
    split2=rep(0,B)
    split3=rep(0,B)

    Indcond1=matrix(rep(NA,p1*B),nrow=B)
    Indcond2=matrix(rep(NA,p1*B),nrow=B)
    Indcond3=matrix(rep(NA,p1*B),nrow=B)




    #########################S Matrices!!!
    #Reset up lam and S1 matrices
    s1=matrix(rep(NA,B*(J1max+2)),nrow=B)
    s1[1,1:(J1+2)]=sort(seq(0,m1,length.out = J1+2))
    s2=matrix(rep(NA,B*(J2max+2)),nrow=B)
    s2[1,1:(J2+2)]=sort(seq(0,m2,length.out = J2+2))
    s3=matrix(rep(NA,B*(J3max+2)),nrow=B)
    s3[1,1:(J3+2)]=sort(seq(0,m3,length.out = J3+2))


    lam1=matrix(rep(NA,B*(J1max+1)),nrow=B)
    lam2=matrix(rep(NA,B*(J2max+1)),nrow=B)
    lam3=matrix(rep(NA,B*(J3max+1)),nrow=B)
    lam1[1,1:(J1+1)]=rep(0,J1+1)
    lam2[1,1:(J2+1)]=rep(0,J2+1)
    lam3[1,1:(J3+1)]=rep(0,J3+1)

    ###Acceptance
    split1=rep(0,B)
    split2=rep(0,B)
    split3=rep(0,B)
    ####Birth
    IndB1=rep(0,B)
    IndB2=rep(0,B)
    IndB3=rep(0,B)
    ###Death
    IndD1=rep(0,B)
    IndD2=rep(0,B)
    IndD3=rep(0,B)
    Indeta1=rep(0,B)
    Indeta2=rep(0,B)
    Indeta3=rep(0,B)

    Ind1s=rep(0,B)
    Ind2s=rep(0,B)
    Ind3s=rep(0,B)
    Indcor1=rep(0,B)
    Indcor2=rep(0,B)
    Indcor3=rep(0,B)



    ######################################
    ### Function########################33
    ####################################

    ## MAtrices for storing



    ######################################
    ### Function########################33
    ####################################

    ## MAtrices for storing



    n=length(Y1)
    G1=J1+1
    G2=J2+1
    G3=J3+1





    ###Haz 2





    ###








    #####
    LK1L=function(Y1,Y2,I1,I2,X,Beta1,Beta2,Beta3,s1,s2,s3,lam1,lam2,lam3,gam){

      LOGBH=0
      et1=X%*%Beta1


      for(k in 1:G1){


        Del=pmax(0,pmin(Y1,s1[k+1])-s1[k])



        LOGBH=LOGBH-sum(gam*Del*exp(lam1[k])*exp(et1))

        zu=Y1<=s1[k+1]
        zl=Y1>s1[k]
        LOGBH=LOGBH+sum(zu*zl*I1)*lam1[k]
      }



      return(LOGBH)

    }
    ###Haz 2


    LK2L=function(Y1,Y2,I1,I2,X,Beta1,Beta2,Beta3,s1,s2,s3,lam1,lam2,lam3,gam){

      LOGBH=0
      et1=X%*%Beta2

      Y=Y1
      Y[I1==0]=Y2[I1==0]

      for(k in 1:G2){


        Del=pmax(0,pmin(Y,s2[k+1])-s2[k])



        LOGBH=LOGBH-sum(gam*Del*exp(lam2[k])*exp(et1))


        zu=Y1<=s2[k+1]
        zl=Y1>s2[k]
        LOGBH=LOGBH+sum(zu*zl*I2*(1-I1))*lam2[k]

      }

      return(LOGBH)
    }



    ###

    LK3L=function(Y1,Y2,I1,I2,X,Beta1,Beta2,Beta3,s1,s2,s3,lam1,lam2,lam3,gam){

      LOGBH=0
      et1=X%*%Beta3
      for(k in 1:G3){


        Del=pmax(0,pmin(Y2[I1==1]-Y1[I1==1],s3[k+1])-s3[k])



        LOGBH=LOGBH-sum(gam[I1==1]*Del*exp(lam3[k])*exp(et1[I1==1]))

        zu=Y2-Y1<=s3[k+1]
        zl=Y2-Y1>s3[k]
        LOGBH=LOGBH+sum(zu*zl*I2*I1)*lam3[k]

      }

      return(LOGBH)
    }


    ###



    #####
    LK1=function(Y1,Y2,I1,I2,X,Beta1,Beta2,Beta3,s1,s2,s3,lam1,lam2,lam3,gam){

      LOGBH=0
      et1=X%*%Beta1


      for(k in 1:G1){


        Del=pmax(0,pmin(Y1,s1[k+1])-s1[k])



        LOGBH=LOGBH-sum(gam*Del*exp(lam1[k])*exp(et1))

      }

      LOGBH=LOGBH+sum(I1*et1)



      return(LOGBH)

    }
    ###Haz 2


    LK2=function(Y1,Y2,I1,I2,X,Beta1,Beta2,Beta3,s1,s2,s3,lam1,lam2,lam3,gam){

      LOGBH=0
      et1=X%*%Beta2

      Y=Y1
      Y[I1==0]=Y2[I1==0]


      for(k in 1:G2){


        Del=pmax(0,pmin(Y,s2[k+1])-s2[k])



        LOGBH=LOGBH-sum(gam*Del*exp(lam2[k])*exp(et1))

      }

      LOGBH=LOGBH+sum(I2*(1-I1)*et1)


      return(LOGBH)
    }



    ###

    LK3=function(Y1,Y2,I1,I2,X,Beta1,Beta2,Beta3,s1,s2,s3,lam1,lam2,lam3,gam){

      LOGBH=0
      et1=X%*%Beta3
      for(k in 1:G3){


        Del=pmax(0,pmin(Y2[I1==1]-Y1[I1==1],s3[k+1])-s3[k])



        LOGBH=LOGBH-sum(gam[I1==1]*Del*exp(lam3[k])*exp(et1[I1==1]))

      }

      LOGBH=LOGBH+sum(I1*I2*et1)


      return(LOGBH)
    }






    D1=function(ep,gamma){
      D=n*log(ep)+(n*ep+psi1-1)/ep-sum(gamma)-w+n*dgamma(ep,1,1)+sum(log(gamma))
      return(D)
    }

    D2=function(ep){
      D=n/ep + (psi1-1)/(ep^2)-n*trigamma(ep)
      return(D)
    }


    mep=function(ep,gamma){
      D=ep-min(0,D1(ep,gamma)/D2(ep))
      return(D)
    }

    vep=function(ep){
      D=-(cep^2)/D2(ep)
      return(D)
    }


    ###Phifunction



    phifun=function(Y1,Y2,I1,I2,B1,B2,B3,S1,S2,S3,Lam1,Lam2,Lam3,Ep,X){
      Ep=1/Ep
      et1=exp(X%*%B1)
      et2=exp(X%*%B2)
      et3=exp(X%*%B3)
      n=length(Y1)
      phi=rep(0,n)
      In2=rep(0,G3)
      In3=rep(0,G3)
      for(i in 1:n){
        if(I1[i]==0 & I2[i]==0){
          delta1=rep(0,G1)
          delta2=rep(0,G2)


          for(m in 1:G1){delta1[m]=max(0,min(Y1[i],S1[m+1])-S1[m])}
          for(m in 1:G2){delta2[m]=max(0,min(Y2[i],S2[m+1])-S2[m])}


          phi[i]=et1[i]*(t(exp(Lam1))%*%as.matrix(delta1)) + et2[i]*(t(exp(Lam2))%*%as.matrix(delta2))+Ep

        }




        ###Case2###
        if(I1[i]==1 & I2[i]==0){
          delta1=rep(0,G1)
          delta2=rep(0,G2)
          delta3=rep(0,G3)

          for(m in 1:G3){delta3[m]=max(0,min(Y2[i]-Y1[i],S3[m+1])-S3[m])
          In3[m]=Y1[i]<=S3[m+1]}
          for(m in 1:G1){delta1[m]=max(0,min(Y1[i],S1[m+1])-S1[m])}
          for(m in 1:G2){delta2[m]=max(0,min(Y1[i],S2[m+1])-S2[m])}

          delta3=delta3*In3

          phi[i]=et1[i]*(t(exp(Lam1))%*%as.matrix(delta1))+et2[i]*(t(exp(Lam2))%*%as.matrix(delta2))+
            et3[i]*(t(exp(Lam3))%*%as.matrix(delta3))+Ep
        }
        ###Case 3###
        if(I1[i]==0 & I2[i]==1){
          delta1=rep(0,G1)
          delta2=rep(0,G2)

          for(m in 1:G1){delta1[m]=max(0,min(Y1[i],S1[m+1])-S1[m])}
          for(m in 1:G2){delta2[m]=max(0,min(Y2[i],S2[m+1])-S2[m])}

          phi[i]=et1[i]*(t(exp(Lam1))%*%as.matrix(delta1)) + et2[i]*(t(exp(Lam2))%*%as.matrix(delta2))+Ep
        }
        ###Case 4###
        if(I1[i]==1 & I2[i]==1){


          delta1=rep(0,G1)
          delta2=rep(0,G2)
          delta3=rep(0,G3)
          In3=rep(0,G3)
          In2=rep(0,G3)


          for(m in 1:G3){delta3[m]=max(0,min(Y2[i]-Y1[i],S3[m+1])-S3[m])
          In3[m]=Y1[i]<=S3[m+1]
          In2[m]=Y2[i]<=S3[m+1]
          }
          for(m in 1:G1){delta1[m]=max(0,min(Y1[i],S1[m+1])-S1[m])}
          for(m in 1:G2){delta2[m]=max(0,min(Y2[i],S2[m+1])-S2[m])}


          delta3=delta3*In3*In2

          phi[i]=et1[i]*(t(exp(Lam1))%*%as.matrix(delta1))+et2[i]*(t(exp(Lam2))%*%as.matrix(delta2))+
            et3[i]*(t(exp(Lam3))%*%as.matrix(delta3))+Ep
        }

      }
      return(phi)
    }


    if(inc>1){
      cat("More than One Variable Included", "


          ")

      ###Set Up Additional Acceptance Matrix

      IncCond1=matrix(rep(0,B*inc),nrow=B)
      IncCond2=matrix(rep(0,B*inc),nrow=B)
      IncCond3=matrix(rep(0,B*inc),nrow=B)


      iter=c(0,0)


      ##Sampler


      for(b in 2:B){





        if(b%%10000==0){cat(b, "iterations",date(), "  ")}else{
          if(b%%5000==0){cat(b, " iterations ")}}

        U=runif(1,0,1)


        iter[1]="etabeta1"

        ###eta1,beta1
        eta1[b,]=eta1[b-1,]
        beta1[b,]=beta1[b-1,]

        if(sum(eta1[b-1,])==0|sum(eta1[b-1,])==p1){
          if(sum(eta1[b-1,])==0){
            ###Add Automatically
            iter[2]="Add"
            Ind=sample(1:p1,1)
            eta1[b,Ind]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ####
            V1 = Sigmanew[spot2,spot2]
            V2 = as.matrix(Sigmanew[-spot2,-spot2])
            V12 = as.matrix(Sigmanew[spot2,-spot2])
            thetab=beta1[b-1,c(includednew, (p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot2])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            beta1[b,Ind]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta1[b,Ind],meannew,varnew))
            ######Accept reject###

            Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b-1,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo+dn + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
            U=log(runif(1,0,1))

            if(is.finite(alphab1)==FALSE){
              eta1[b,]=eta1[b-1,]
              beta1[b,]=beta1[b-1,]
              Indeta1[b]=0
            }else{

              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}
            }
          }
          if(sum(eta1[b-1,])==p1){
            ###Delete Automatically
            Ind=sample(1:p1,1)
            iter[2]="delete"
            eta1[b,Ind]=0
            beta1[b,Ind]=0
            includedold=rep(0,p1)
            for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(includedold[k]==Ind){spotold[k]=k}}
            spot1=max(spotold)


            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])

            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta1[b-1,c(includedold,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta1[b-1,Ind],meanold,varold))
            ######Accept reject###
            Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b-1,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo-do + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta1[b,]=eta1[b-1,]
              beta1[b,]=beta1[b-1,]
              Indeta1[b]=0
            }else{

              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}
            }}
        }else{

          U=runif(1,0,1)

          if(U<psi){
            ###Swapper
            includedold=rep(0,p1)
            iter[2]="swap"
            for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            ones=includedold
            zeros=rep(0,p1)
            for(k in 1:p1){if(eta1[b-1,k]==0){zeros[k]=k}}
            zeros=zeros[zeros != 0]
            ###Sample swap indices###
            if(length(ones)==1){
              Indone=ones}else{
                Indone=sample(ones,1)}
            if(length(zeros)==1){Indzero=zeros}else{
              Indzero=sample(zeros,1)}
            ####Change Beta/eta
            eta1[b,Indone]=0
            eta1[b,Indzero]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(Indone==includedold[k]){spotold[k]=k}}
            spot1=max(spotold)
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Indzero==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ###Generate new vector##
            beta1[b,Indone]=0

            ##meannew,varnew##
            V1 = Sigmanew[spot2,spot2]
            V2 = as.matrix(Sigmanew[-spot2,-spot2])
            V12 = as.matrix(Sigmanew[spot2,-spot2])
            thetab=beta1[b-1,c(includednew,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot2])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            beta1[b,Indzero]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta1[b,Indzero],meannew,varnew))
            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta1[b-1,c(includedold,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta1[b-1,Indone],meanold,varold))
            ######Accept reject###
            Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b-1,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo+dn-do
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta1[b,]=eta1[b-1,]
              beta1[b,]=beta1[b-1,]
              Indeta1[b]=0
            }else{
              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}

            }

          }else{
            ###Add/Delete
            Ind=sample(1:p1,1)
            if(eta1[b-1,Ind]==1){
              ##delete##
              iter[2]="delete"
              eta1[b,Ind]=0
              beta1[b,Ind]=0
              includedold=rep(0,p1)
              for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
              includedold=includedold[includedold != 0]
              spotold=rep(0,length(includedold))
              for(k in 1:length(includedold)){if(Ind==includedold[k]){spotold[k]=k}}
              spot1=max(spotold)


              ###Make sigma matrices##
              Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])

              ###Old density###
              V1 = Sigmaold[spot1,spot1]
              V2 = as.matrix(Sigmaold[-spot1,-spot1])
              V12 = as.matrix(Sigmaold[spot1,-spot1])
              thetab=beta1[b-1,c(includedold,(p1+1):(p1+inc))]
              thetano = as.matrix(thetab[-spot1])
              meanold = t(V12)%*%solve(V2)%*%thetano
              varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              do=log(dnorm(beta1[b-1,Ind],meanold,varold))
              ######Accept reject###
              Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b-1,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo-do + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
              U=log(runif(1,0,1))
              if(is.finite(alphab1)==FALSE){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{
                if(U>alphab1){
                  eta1[b,]=eta1[b-1,]
                  beta1[b,]=beta1[b-1,]
                  Indeta1[b]=0
                }else{Indeta1[b]=1}

              }

            }else{
              ###Add###
              eta1[b,Ind]=1

              iter[2]="add"
              includednew=rep(0,p1)
              for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
              includednew=includednew[includednew != 0]
              spotnew=rep(0,length(includednew))
              for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
              spot2=max(spotnew)
              ###Make sigma matrices##
              Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
              ####
              V1 = Sigmanew[spot2,spot2]
              V2 = as.matrix(Sigmanew[-spot2,-spot2])
              V12 = as.matrix(Sigmanew[spot2,-spot2])
              thetab=beta1[b-1,c(includednew,(p1+1):(p1+inc))]
              thetano = as.matrix(thetab[-spot2])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              beta1[b,Ind]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta1[b,Ind],meannew,varnew))
              ######Accept reject###
              Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b-1,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo+dn + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
              U=log(runif(1,0,1))
              if(is.finite(alphab1)==FALSE){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{
                if(U>alphab1){
                  eta1[b,]=eta1[b-1,]
                  beta1[b,]=beta1[b-1,]
                  Indeta1[b]=0
                }else{Indeta1[b]=1}
              }


            }

          }}



        ####ETABETA 2


        iter[1]="etabeta2"

        ###eta1,beta1
        eta2[b,]=eta2[b-1,]
        beta2[b,]=beta2[b-1,]

        if(sum(eta2[b-1,])==0|sum(eta2[b-1,])==p1){
          if(sum(eta2[b-1,])==0){
            ###Add Automatically
            iter[2]="Add"
            Ind=sample(1:p1,1)
            eta2[b,Ind]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta2[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ####
            V1 = Sigmanew[spot2,spot2]
            V2 = as.matrix(Sigmanew[-spot2,-spot2])
            V12 = as.matrix(Sigmanew[spot2,-spot2])
            thetab=beta2[b-1,c(includednew, (p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot2])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            beta2[b,Ind]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta2[b,Ind],meannew,varnew))
            ######Accept reject###

            Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo+dn + log(beta(sum(eta2[b,])+z1a,p1-sum(eta2[b,])+z1b)) - log(beta(sum(eta2[b-1,])+z1a,p1-sum(eta2[b-1,])+z1b))
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta2[b,]=eta2[b-1,]
              beta2[b,]=beta2[b-1,]
              Indeta2[b]=0
            }else{
              if(U>alphab1){
                eta2[b,]=eta2[b-1,]
                beta2[b,]=beta2[b-1,]
                Indeta2[b]=0
              }else{Indeta2[b]=1}
            }}
          if(sum(eta2[b-1,])==p1){
            ###Delete Automatically
            Ind=sample(1:p1,1)
            iter[2]="delete"
            eta2[b,Ind]=0
            beta2[b,Ind]=0
            includedold=rep(0,p1)
            for(k in 1:p1){if(eta2[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(includedold[k]==Ind){spotold[k]=k}}
            spot1=max(spotold)


            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])

            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta2[b-1,c(includedold,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta2[b-1,Ind],meanold,varold))
            ######Accept reject###
            Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo-do + log(beta(sum(eta2[b,])+z1a,p1-sum(eta2[b,])+z1b)) - log(beta(sum(eta2[b-1,])+z1a,p1-sum(eta2[b-1,])+z1b))
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta2[b,]=eta2[b-1,]
              beta2[b,]=beta2[b-1,]
              Indeta2[b]=0
            }else{
              if(U>alphab1){
                eta2[b,]=eta2[b-1,]
                beta2[b,]=beta2[b-1,]
                Indeta2[b]=0
              }else{Indeta2[b]=1}
            }}
        }else{

          U=runif(1,0,1)

          if(U<psi){
            ###Swapper
            includedold=rep(0,p1)
            iter[2]="swap"
            for(k in 1:p1){if(eta2[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            ones=includedold
            zeros=rep(0,p1)
            for(k in 1:p1){if(eta2[b-1,k]==0){zeros[k]=k}}
            zeros=zeros[zeros != 0]
            ###Sample swap indices###
            if(length(ones)==1){
              Indone=ones}else{
                Indone=sample(ones,1)}
            if(length(zeros)==1){Indzero=zeros}else{
              Indzero=sample(zeros,1)}
            ####Change Beta/eta
            eta2[b,Indone]=0
            eta2[b,Indzero]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta2[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(Indone==includedold[k]){spotold[k]=k}}
            spot1=max(spotold)
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Indzero==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ###Generate new vector##
            beta2[b,Indone]=0

            ##meannew,varnew##
            V1 = Sigmanew[spot2,spot2]
            V2 = as.matrix(Sigmanew[-spot2,-spot2])
            V12 = as.matrix(Sigmanew[spot2,-spot2])
            thetab=beta2[b-1,c(includednew,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot2])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            beta2[b,Indzero]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta2[b,Indzero],meannew,varnew))
            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta2[b-1,c(includedold,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta2[b-1,Indone],meanold,varold))
            ######Accept reject###
            Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo+dn-do
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta2[b,]=eta2[b-1,]
              beta2[b,]=beta2[b-1,]
              Indeta2[b]=0
            }else{
              if(U>alphab1){
                eta2[b,]=eta2[b-1,]
                beta2[b,]=beta2[b-1,]
                Indeta2[b]=0
              }else{Indeta2[b]=1}
            }


          }else{
            ###Add/Delete
            Ind=sample(1:p1,1)
            if(eta2[b-1,Ind]==1){
              ##delete##
              iter[2]="delete"
              eta2[b,Ind]=0
              beta2[b,Ind]=0
              includedold=rep(0,p1)
              for(k in 1:p1){if(eta2[b-1,k]==1){includedold[k]=k}}
              includedold=includedold[includedold != 0]
              spotold=rep(0,length(includedold))
              for(k in 1:length(includedold)){if(Ind==includedold[k]){spotold[k]=k}}
              spot1=max(spotold)


              ###Make sigma matrices##
              Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])

              ###Old density###
              V1 = Sigmaold[spot1,spot1]
              V2 = as.matrix(Sigmaold[-spot1,-spot1])
              V12 = as.matrix(Sigmaold[spot1,-spot1])
              thetab=beta2[b-1,c(includedold,(p1+1):(p1+inc))]
              thetano = as.matrix(thetab[-spot1])
              meanold = t(V12)%*%solve(V2)%*%thetano
              varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              do=log(dnorm(beta2[b-1,Ind],meanold,varold))
              ######Accept reject###
              Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo-do + log(beta(sum(eta2[b,])+z1a,p1-sum(eta2[b,])+z1b)) - log(beta(sum(eta2[b-1,])+z1a,p1-sum(eta2[b-1,])+z1b))
              U=log(runif(1,0,1))
              if(is.finite(alphab1)==FALSE){
                eta2[b,]=eta2[b-1,]
                beta2[b,]=beta2[b-1,]
                Indeta2[b]=0
              }else{
                if(U>alphab1){
                  eta2[b,]=eta2[b-1,]
                  beta2[b,]=beta2[b-1,]
                  Indeta2[b]=0
                }else{Indeta2[b]=1}
              }


            }else{
              ###Add###
              eta2[b,Ind]=1

              iter[2]="add"
              includednew=rep(0,p1)
              for(k in 1:p1){if(eta2[b,k]==1){includednew[k]=k}}
              includednew=includednew[includednew != 0]
              spotnew=rep(0,length(includednew))
              for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
              spot2=max(spotnew)
              ###Make sigma matrices##
              Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
              ####
              V1 = Sigmanew[spot2,spot2]
              V2 = as.matrix(Sigmanew[-spot2,-spot2])
              V12 = as.matrix(Sigmanew[spot2,-spot2])
              thetab=beta2[b-1,c(includednew,(p1+1):(p1+inc))]
              thetano = as.matrix(thetab[-spot2])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              beta2[b,Ind]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta2[b,Ind],meannew,varnew))
              ######Accept reject###
              Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo+dn + log(beta(sum(eta2[b,])+z1a,p1-sum(eta2[b,])+z1b)) - log(beta(sum(eta2[b-1,])+z1a,p1-sum(eta2[b-1,])+z1b))
              U=log(runif(1,0,1))
              if(is.finite(alphab1)==FALSE){
                eta2[b,]=eta2[b-1,]
                beta2[b,]=beta2[b-1,]
                Indeta2[b]=0
              }else{
                if(U>alphab1){
                  eta2[b,]=eta2[b-1,]
                  beta2[b,]=beta2[b-1,]
                  Indeta2[b]=0
                }else{Indeta2[b]=1}
              }


            }

          }}



        #####ETA3###


        ####

        iter[1]="etabeta3"

        ###eta1,beta1
        eta3[b,]=eta3[b-1,]
        beta3[b,]=beta3[b-1,]

        if(sum(eta3[b-1,])==0|sum(eta3[b-1,])==p1){
          if(sum(eta3[b-1,])==0){
            ###Add Automatically
            iter[2]="Add"
            Ind=sample(1:p1,1)
            eta3[b,Ind]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta3[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ####
            V1 = Sigmanew[spot2,spot2]
            V2 = as.matrix(Sigmanew[-spot2,-spot2])
            V12 = as.matrix(Sigmanew[spot2,-spot2])
            thetab=beta3[b-1,c(includednew, (p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot2])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            beta3[b,Ind]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta3[b,Ind],meannew,varnew))
            ######Accept reject###

            Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo+dn + log(beta(sum(eta3[b,])+z1a,p1-sum(eta3[b,])+z1b)) - log(beta(sum(eta3[b-1,])+z1a,p1-sum(eta3[b-1,])+z1b))
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta3[b,]=eta3[b-1,]
              beta3[b,]=beta3[b-1,]
              Indeta3[b]=0
            }else{
              if(U>alphab1){
                eta3[b,]=eta3[b-1,]
                beta3[b,]=beta3[b-1,]
                Indeta3[b]=0
              }else{Indeta3[b]=1}
            }}
          if(sum(eta3[b-1,])==p1){
            ###Delete Automatically
            Ind=sample(1:p1,1)
            iter[2]="delete"
            eta3[b,Ind]=0
            beta3[b,Ind]=0
            includedold=rep(0,p1)
            for(k in 1:p1){if(eta3[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(includedold[k]==Ind){spotold[k]=k}}
            spot1=max(spotold)


            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])

            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta3[b-1,c(includedold,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta3[b-1,Ind],meanold,varold))
            ######Accept reject###
            Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo-do + log(beta(sum(eta3[b,])+z1a,p1-sum(eta3[b,])+z1b)) - log(beta(sum(eta3[b-1,])+z1a,p1-sum(eta3[b-1,])+z1b))
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta3[b,]=eta3[b-1,]
              beta3[b,]=beta3[b-1,]
              Indeta3[b]=0
            }else{
              if(U>alphab1){
                eta3[b,]=eta3[b-1,]
                beta3[b,]=beta3[b-1,]
                Indeta3[b]=0
              }else{Indeta3[b]=1}
            }}
        }else{

          U=runif(1,0,1)

          if(U<psi){
            ###Swapper
            includedold=rep(0,p1)
            iter[2]="swap"
            for(k in 1:p1){if(eta3[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            ones=includedold
            zeros=rep(0,p1)
            for(k in 1:p1){if(eta3[b-1,k]==0){zeros[k]=k}}
            zeros=zeros[zeros != 0]
            ###Sample swap indices###
            if(length(ones)==1){
              Indone=ones}else{
                Indone=sample(ones,1)}
            if(length(zeros)==1){Indzero=zeros}else{
              Indzero=sample(zeros,1)}
            ####Change Beta/eta
            eta3[b,Indone]=0
            eta3[b,Indzero]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta3[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(Indone==includedold[k]){spotold[k]=k}}
            spot1=max(spotold)
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Indzero==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ###Generate new vector##
            beta3[b,Indone]=0

            ##meannew,varnew##
            V1 = Sigmanew[spot2,spot2]
            V2 = as.matrix(Sigmanew[-spot2,-spot2])
            V12 = as.matrix(Sigmanew[spot2,-spot2])
            thetab=beta3[b-1,c(includednew,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot2])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            beta3[b,Indzero]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta3[b,Indzero],meannew,varnew))
            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta3[b-1,c(includedold,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta3[b-1,Indone],meanold,varold))
            ######Accept reject###
            Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo+dn-do
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta3[b,]=eta3[b-1,]
              beta3[b,]=beta3[b-1,]
              Indeta3[b]=0
            }else{
              if(U>alphab1){
                eta3[b,]=eta3[b-1,]
                beta3[b,]=beta3[b-1,]
                Indeta3[b]=0
              }else{Indeta3[b]=1}}



          }else{
            ###Add/Delete
            Ind=sample(1:p1,1)
            if(eta3[b-1,Ind]==1){
              ##delete##
              iter[2]="delete"
              eta3[b,Ind]=0
              beta3[b,Ind]=0
              includedold=rep(0,p1)
              for(k in 1:p1){if(eta3[b-1,k]==1){includedold[k]=k}}
              includedold=includedold[includedold != 0]
              spotold=rep(0,length(includedold))
              for(k in 1:length(includedold)){if(Ind==includedold[k]){spotold[k]=k}}
              spot1=max(spotold)


              ###Make sigma matrices##
              Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])

              ###Old density###
              V1 = Sigmaold[spot1,spot1]
              V2 = as.matrix(Sigmaold[-spot1,-spot1])
              V12 = as.matrix(Sigmaold[spot1,-spot1])
              thetab=beta3[b-1,c(includedold,(p1+1):(p1+inc))]
              thetano = as.matrix(thetab[-spot1])
              meanold = t(V12)%*%solve(V2)%*%thetano
              varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              do=log(dnorm(beta3[b-1,Ind],meanold,varold))
              ######Accept reject###
              Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo-do + log(beta(sum(eta3[b,])+z1a,p1-sum(eta3[b,])+z1b)) - log(beta(sum(eta3[b-1,])+z1a,p1-sum(eta3[b-1,])+z1b))
              U=log(runif(1,0,1))
              if(is.finite(alphab1)==FALSE){
                eta3[b,]=eta3[b-1,]
                beta3[b,]=beta3[b-1,]
                Indeta3[b]=0
              }else{
                if(U>alphab1){
                  eta3[b,]=eta3[b-1,]
                  beta3[b,]=beta3[b-1,]
                  Indeta3[b]=0
                }else{Indeta3[b]=1}}



            }else{
              ###Add###
              eta3[b,Ind]=1

              iter[2]="add"
              includednew=rep(0,p1)
              for(k in 1:p1){if(eta3[b,k]==1){includednew[k]=k}}
              includednew=includednew[includednew != 0]
              spotnew=rep(0,length(includednew))
              for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
              spot2=max(spotnew)
              ###Make sigma matrices##
              Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
              ####
              V1 = Sigmanew[spot2,spot2]
              V2 = as.matrix(Sigmanew[-spot2,-spot2])
              V12 = as.matrix(Sigmanew[spot2,-spot2])
              thetab=beta3[b-1,c(includednew,(p1+1):(p1+inc))]
              thetano = as.matrix(thetab[-spot2])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              beta3[b,Ind]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta3[b,Ind],meannew,varnew))
              ######Accept reject###
              Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo+dn + log(beta(sum(eta3[b,])+z1a,p1-sum(eta3[b,])+z1b)) - log(beta(sum(eta3[b-1,])+z1a,p1-sum(eta3[b-1,])+z1b))
              U=log(runif(1,0,1))
              if(is.finite(alphab1)==FALSE){
                eta3[b,]=eta3[b-1,]
                beta3[b,]=beta3[b-1,]
                Indeta3[b]=0
              }else{
                if(U>alphab1){
                  eta3[b,]=eta3[b-1,]
                  beta3[b,]=beta3[b-1,]
                  Indeta3[b]=0
                }else{Indeta3[b]=1}}



            }

          }}








        ###INCLUDED SAMPLERS

        iter[1]="Beta1"
        iter[2]="Included"

        if(sum(eta1[b,])==0){
          ##Sample Included
          Sigmanew= c*solve(t(X[,(p1+1):(p1+inc)])%*%X[,(p1+1):(p1+inc)])
          zeta1n=beta1[b,(p1+1):(p1+inc)]
          for(k in 1:inc){
            zeta1=zeta1n
            V1 = Sigmanew[k,k]
            V2 = as.matrix(Sigmanew[-k,-k])
            V12 = as.matrix(Sigmanew[k,-k])
            thetano=zeta1[-k]
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            zeta1[k]=rnorm(1,meannew,varnew)
            dn=log(dnorm(zeta1[k],meannew,varnew))
            ###density old
            do=log(dnorm(zeta1n[k],meannew,varnew))
            beta=beta1[b,]
            beta[(p1+1):(p1+inc)]=zeta1
            Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK1(Y1,Y2,I1,I2,X,beta,beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1m=Liken-Likeo+dn -do
            U=log(runif(1,0,1))

            if(is.finite(alphab1m)==FALSE){
              IncCond1[b,k]=0
            }else{
              if(U>alphab1m){
                IncCond1[b,k]=0
              }else{IncCond1[b,k]=1
              beta1[b,]=beta
              zeta1n=zeta1
              }}
            ##End Inc Sampler
          } }else{
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            zeta1n=beta1[b,c(includednew,(p1+1):(p1+inc))]

            ###Make sigma matrices##
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ####
            p=length(includednew)+inc
            ####Update All included variables

            for(k in (length(includednew)+1):(length(includednew)+inc)){
              zeta1=zeta1n
              V1 = Sigmanew[k,k]
              V2 = as.matrix(Sigmanew[-k,-k])
              V12 = as.matrix(Sigmanew[k,-k])
              thetano = as.matrix(zeta1[-k])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              zeta1[k]=rnorm(1,meannew,varnew)
              dn=log(dnorm(zeta1[k],meannew,varnew))
              ###density old
              do=log(dnorm(beta1[b,(p1+k-length(includednew))],meannew,varnew))


              ######Accept reject###
              Likeo=LK1(Y1,Y2,I1,I2,X,c(beta1[b,1:p1],zeta1n[(length(zeta1n)-inc+1):length(zeta1n)]),beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK1(Y1,Y2,I1,I2,X,c(beta1[b,1:p1],zeta1[(length(zeta1n)-inc+1):length(zeta1n)]),beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1s=Liken-Likeo+dn -do
              U=log(runif(1,0,1))

              if(is.finite(alphab1s)==FALSE){
                IncCond1[b,(k-p1)]=0

              }else{

                if(U>alphab1s){



                  IncCond1[b,(k-p1)]=0

                }else{IncCond1[b,(k-p1)]=1
                zeta1n=zeta1
                beta1[b,]=c(beta1[b,1:p1],zeta1[(length(zeta1)-inc+1):length(zeta1)])

                }
              }

            }

            ###End included sampler###
          }



        #####Conditional Sampler for Included!###


        if(sum(eta1[b,])>0){

          iter[2]="Conditional Inclusion"
          ##Jointly Update nonzero betas
          zeta1=beta1[b,]
          zeta1=zeta1[zeta1!=0]
          zeta1n=zeta1
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])



          ###############
          ####

          for(k in 1:length(includednew)){


            V1 = Sigmanew[k,k]
            V2 = as.matrix(Sigmanew[-k,-k])
            V12 = as.matrix(Sigmanew[k,-k])
            thetab=beta1[b,c(includednew,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-k])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            zeta1n[k]=rnorm(1,meannew,varnew)
            dn=log(dnorm(zeta1n[k],meannew,varnew))
            ###density old
            do=log(dnorm(zeta1[k],meannew,varnew))






            beta=beta1[b,]
            beta[c(includednew,(p1+1):(p1+inc))]=zeta1n


            Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK1(Y1,Y2,I1,I2,X,beta,beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1m=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphab1m)==FALSE){
              Indcond1[b,k]=0
            }else{
              if(U>alphab1m){
                Indcond1[b,includednew[k]]=0
                zeta1n[k]=zeta1[k]
              }else{Indcond1[b,includednew[k]]=1
              beta1[b,]=beta
              zeta1[k]=zeta1n[k]
              }}

          }



          ##Jointly Update nonzero betas
          iter[2]="mixing"
          zeta1n=beta1[b,]
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])


          zeta1n[c(includednew,(p1+1):(p1+inc))]=rmvnorm(1,rep(0,length(includednew)+inc),Sigmanew)

          beta=beta1[b,]
          beta=beta[beta!=0]

          dn=log(dmvnorm(zeta1n[c(includednew,(p1+1):(p1+inc))],rep(0,length(includednew)+inc),Sigmanew))
          ###density old
          do=log(dmvnorm(beta,rep(0,length(includednew)+inc),Sigmanew))

          ######Accept reject###
          Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
          Liken=LK1(Y1,Y2,I1,I2,X,zeta1n,beta2[b-1,],beta3[b-1,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

          alphamix1=Liken-Likeo+dn -do
          U=log(runif(1,0,1))
          if(is.finite(alphamix1)==FALSE){
            Indmix1[b]=0
          }else{
            if(U>alphamix1){

              Indmix1[b]=0
            }else{Indmix1[b]=1
            beta1[b,]=zeta1n
            }}

        }else{
          ##Jointly Update nonzero betas
          iter[2]="mixing No eta"
          zeta1n=beta1[b,]
          Sigmanew=c*solve(t(X[,(p1+1):(p1+inc)])%*%X[,(p1+1):(p1+inc)])


          zeta1n[(p1+1):(p1+inc)]=rmvnorm(1,rep(0,inc),Sigmanew)

          beta=beta1[b,]
          beta=beta[beta!=0]

          dn=log(dmvnorm(zeta1n[(p1+1):(p1+inc)],rep(0,inc),Sigmanew))
          ###density old
          do=log(dmvnorm(beta,rep(0,inc),Sigmanew))

          ######Accept reject###
          Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
          Liken=LK1(Y1,Y2,I1,I2,X,zeta1n,beta2[b-1,],beta3[b-1,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

          alphamix1=Liken-Likeo+dn -do
          U=log(runif(1,0,1))

          if(is.finite(alphamix1)==FALSE){
            Indmix1[b]=0}else{

              if(U>alphamix1){



                Indmix1[b]=0
              }else{Indmix1[b]=1
              beta1[b,]=zeta1n
              }}

        }





        iter[1]="Beta2"
        iter[2]="Included"

        if(sum(eta2[b,])==0){
          ##Sample Included
          Sigmanew= c*solve(t(X[,(p1+1):(p1+inc)])%*%X[,(p1+1):(p1+inc)])
          zeta1n=beta2[b,(p1+1):(p1+inc)]
          for(k in 1:inc){
            zeta1=zeta1n
            V1 = Sigmanew[k,k]
            V2 = as.matrix(Sigmanew[-k,-k])
            V12 = as.matrix(Sigmanew[k,-k])
            thetano=zeta1[-k]
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            zeta1[k]=rnorm(1,meannew,varnew)
            dn=log(dnorm(zeta1[k],meannew,varnew))
            ###density old
            do=log(dnorm(zeta1n[k],meannew,varnew))
            beta=beta2[b,]
            beta[(p1+1):(p1+inc)]=zeta1
            Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta,beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1m=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphab1m)==FALSE){
              IncCond2[b,k]=0
            }else{
              if(U>alphab1m){
                IncCond2[b,k]=0
              }else{IncCond2[b,k]=1
              beta2[b,]=beta
              zeta1n=zeta1
              }}
            ##End Inc Sampler
          } }else{
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta2[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            zeta1n=beta2[b,c(includednew,(p1+1):(p1+inc))]

            ###Make sigma matrices##
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ####
            p=length(includednew)+inc
            ####Update All included variables

            for(k in (length(includednew)+1):(length(includednew)+inc)){
              zeta1=zeta1n
              V1 = Sigmanew[k,k]
              V2 = as.matrix(Sigmanew[-k,-k])
              V12 = as.matrix(Sigmanew[k,-k])
              thetano = as.matrix(zeta1[-k])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              zeta1[k]=rnorm(1,meannew,varnew)
              dn=log(dnorm(zeta1[k],meannew,varnew))
              ###density old
              do=log(dnorm(beta2[b,(p1+k-length(includednew))],meannew,varnew))



              ######Accept reject###
              Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],c(beta2[b,1:p1],zeta1n[(length(zeta1n)-inc+1):length(zeta1n)]),beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],c(beta2[b,1:p1],zeta1[(length(zeta1n)-inc+1):length(zeta1n)]),beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1s=Liken-Likeo+dn -do
              U=log(runif(1,0,1))

              if(is.finite(alphab1s)==FALSE){
                IncCond2[b,(k-p1)]=0

              }else{

                if(U>alphab1s){



                  IncCond2[b,(k-p1)]=0

                }else{IncCond2[b,(k-p1)]=1
                zeta1n=zeta1
                beta2[b,]=c(beta2[b,1:p1],zeta1[(length(zeta1)-inc+1):length(zeta1)])

                }}

            }

            ###End included sampler###
          }



        #####Conditional Sampler for Included!###

        if(sum(eta2[b,])>0){


          iter[2]="Conditional Inclusion"
          ##Jointly Update nonzero betas
          zeta1=beta2[b,]
          zeta1=zeta1[zeta1!=0]
          zeta1n=zeta1
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])



          ###############
          ####

          for(k in 1:length(includednew)){


            V1 = Sigmanew[k,k]
            V2 = as.matrix(Sigmanew[-k,-k])
            V12 = as.matrix(Sigmanew[k,-k])
            thetab=beta1[b,c(includednew,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-k])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            zeta1n[k]=rnorm(1,meannew,varnew)
            dn=log(dnorm(zeta1n[k],meannew,varnew))
            ###density old
            do=log(dnorm(zeta1[k],meannew,varnew))

            beta=beta2[b,]
            beta[c(includednew,(p1+1):(p1+inc))]=zeta1n


            Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta,beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1m=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphab1m)==FALSE){
              Indcond2[b,k]=0
            }else{
              if(U>alphab1m){
                Indcond2[b,includednew[k]]=0
                zeta1n[k]=zeta1[k]
              }else{Indcond2[b,includednew[k]]=1
              beta2[b,]=beta
              zeta1[k]=zeta1n[k]
              }}

          }



          ##Jointly Update nonzero betas
          iter[2]="mixing"
          zeta1n=beta2[b,]
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])


          zeta1n[c(includednew,(p1+1):(p1+inc))]=rmvnorm(1,rep(0,length(includednew)+inc),Sigmanew)

          beta=beta2[b,]
          beta=beta[beta!=0]

          dn=log(dmvnorm(zeta1n[c(includednew,(p1+1):(p1+inc))],rep(0,length(includednew)+inc),Sigmanew))
          ###density old
          do=log(dmvnorm(beta,rep(0,length(includednew)+inc),Sigmanew))

          ######Accept reject###
          Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
          Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],zeta1n,beta3[b-1,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

          alphamix1=Liken-Likeo+dn -do
          U=log(runif(1,0,1))
          if(is.finite(alphamix1)==FALSE){
            Indmix2[b]=0
          }else{

            if(U>alphamix1){

              Indmix2[b]=0
            }else{Indmix2[b]=1
            beta2[b,]=zeta1n
            }
          }


        }else{
          ##Jointly Update nonzero betas
          iter[2]="mixing no eta"
          zeta1n=beta2[b,]
          Sigmanew=c*solve(t(X[,(p1+1):(p1+inc)])%*%X[,(p1+1):(p1+inc)])


          zeta1n[(p1+1):(p1+inc)]=rmvnorm(1,rep(0,inc),Sigmanew)

          beta=beta2[b,]
          beta=beta[beta!=0]

          dn=log(dmvnorm(zeta1n[(p1+1):(p1+inc)],rep(0,inc),Sigmanew))
          ###density old
          do=log(dmvnorm(beta,rep(0,inc),Sigmanew))

          ######Accept reject###
          Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
          Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],zeta1n,beta3[b-1,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

          alphamix1=Liken-Likeo+dn -do
          U=log(runif(1,0,1))

          if(is.finite(alphamix1)==FALSE){
            Indmix2[b]=0
          }else{

            if(U>alphamix1){

              Indmix2[b]=0
            }else{Indmix2[b]=1
            beta2[b,]=zeta1n
            }}
        }



        iter[1]="Beta3"
        iter[2]="Included"

        if(sum(eta3[b,])==0){
          ##Sample Included
          Sigmanew= c*solve(t(X[,(p1+1):(p1+inc)])%*%X[,(p1+1):(p1+inc)])
          zeta1n=beta3[b,(p1+1):(p1+inc)]
          for(k in 1:inc){
            zeta1=zeta1n
            V1 = Sigmanew[k,k]
            V2 = as.matrix(Sigmanew[-k,-k])
            V12 = as.matrix(Sigmanew[k,-k])
            thetano=zeta1[-k]
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            zeta1[k]=rnorm(1,meannew,varnew)
            dn=log(dnorm(zeta1[k],meannew,varnew))
            ###density old
            do=log(dnorm(zeta1n[k],meannew,varnew))
            beta=beta3[b,]
            beta[(p1+1):(p1+inc)]=zeta1
            Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta,
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1m=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphab1m)==FALSE){
              IncCond3[b,k]=0
            }else{
              if(U>alphab1m){
                IncCond3[b,k]=0
              }else{IncCond2[b,k]=1
              beta3[b,]=beta
              zeta1n=zeta1
              }}
            ##End Inc Sampler
          } }else{
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta3[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            zeta1n=beta3[b-1,c(includednew,(p1+1):(p1+inc))]

            ###Make sigma matrices##
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ####
            p=length(includednew)+inc
            ####Update All included variables

            for(k in (length(includednew)+1):(length(includednew)+inc)){
              zeta1=zeta1n
              V1 = Sigmanew[k,k]
              V2 = as.matrix(Sigmanew[-k,-k])
              V12 = as.matrix(Sigmanew[k,-k])
              thetano = as.matrix(zeta1[-k])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              zeta1[k]=rnorm(1,meannew,varnew)
              dn=log(dnorm(zeta1[k],meannew,varnew))
              ###density old
              do=log(dnorm(beta3[b-1,(p1+k-length(includednew))],meannew,varnew))



              ######Accept reject###
              Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],c(beta3[b,1:p1],zeta1n[(length(zeta1n)-inc+1):length(zeta1n)]),
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],c(beta3[b,1:p1],zeta1[(length(zeta1n)-inc+1):length(zeta1n)]),
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1s=Liken-Likeo+dn -do
              U=log(runif(1,0,1))
              if(is.finite(alphab1s)==FALSE){
                IncCond3[b,(k-p1)]=0

              }else{
                if(U>alphab1s){



                  IncCond3[b,(k-p1)]=0

                }else{IncCond3[b,(k-p1)]=1
                zeta1n=zeta1
                beta3[b,]=c(beta3[b,1:p1],zeta1[(length(zeta1)-inc+1):length(zeta1)])

                }}

            }

            ###End included sampler###
          }



        #####Conditional Sampler for Included!###
        if(sum(eta3[b,])>0){


          iter[2]="Conditional Inclusion"
          ##Jointly Update nonzero betas
          zeta1=beta3[b,]
          zeta1=zeta1[zeta1!=0]
          zeta1n=zeta1
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])



          ###############
          ####

          for(k in 1:length(includednew)){


            V1 = Sigmanew[k,k]
            V2 = as.matrix(Sigmanew[-k,-k])
            V12 = as.matrix(Sigmanew[k,-k])
            thetab=beta1[b,c(includednew,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-k])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            zeta1n[k]=rnorm(1,meannew,varnew)
            dn=log(dnorm(zeta1n[k],meannew,varnew))
            ###density old
            do=log(dnorm(zeta1[k],meannew,varnew))

            beta=beta3[b,]
            beta[c(includednew,(p1+1):(p1+inc))]=zeta1n


            Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta,
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1m=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphab1m)==FALSE){
              Indcond3[b,k]=0
            }else{
              if(U>alphab1m){
                Indcond3[b,includednew[k]]=0
                zeta1n[k]=zeta1[k]
              }else{
                Indcond3[b,includednew[k]]=1
                beta3[b,]=beta
                zeta1[k]=zeta1n[k]
              }}

          }



          ##Jointly Update nonzero betas
          iter[2]="mixing"
          zeta1n=beta3[b,]
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])


          zeta1n[c(includednew,(p1+1):(p1+inc))]=rmvnorm(1,rep(0,length(includednew)+inc),Sigmanew)

          beta=beta3[b,]
          beta=beta[beta!=0]

          dn=log(dmvnorm(zeta1n[c(includednew,(p1+1):(p1+inc))],rep(0,length(includednew)+inc),Sigmanew))
          ###density old
          do=log(dmvnorm(beta,rep(0,length(includednew)+inc),Sigmanew))

          ######Accept reject###
          Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
          Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],zeta1n,
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

          alphamix1=Liken-Likeo+dn -do


          U=log(runif(1,0,1))

          if(is.finite(alphamix1)==FALSE){
            Indmix3[b]=0

          }else{

            if(U>alphamix1){

              Indmix3[b]=0
            }else{Indmix3[b]=1
            beta3[b,]=zeta1n
            }
          }

        }else{
          ##Jointly Update nonzero betas
          iter[2]="mixing no eta"
          zeta1n=beta3[b,]
          Sigmanew=c*solve(t(X[,(p1+1):(p1+inc)])%*%X[,(p1+1):(p1+inc)])


          zeta1n[(p1+1):(p1+inc)]=rmvnorm(1,rep(0,inc),Sigmanew)

          beta=beta3[b,]
          beta=beta[beta!=0]

          dn=log(dmvnorm(zeta1n[(p1+1):(p1+inc)],rep(0,inc),Sigmanew))
          ###density old
          do=log(dmvnorm(beta,rep(0,inc),Sigmanew))

          ######Accept reject###
          Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
          Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],zeta1n,
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

          alphamix1=Liken-Likeo+dn -do
          U=log(runif(1,0,1))

          if(is.finite(alphamix1)==0){
            Indmix3[b]=0}else{


              if(U>alphamix1){

                Indmix3[b]=0
              }else{Indmix3[b]=1
              beta3[b,]=zeta1n
              }}

        }







        ########################
        ###Frailty Samplers#####
        ########################
        ############Epsilon Sampler#####

        iter[1]="frailty"
        iter[2]="hier"
        Der1o=D1(epsilon[b-1],gam[b-1,])
        Der2o=D2(epsilon[b-1])

        epsilon[b]=rgamma(1,((epsilon[b-1]-min(0,Der1o/Der2o))^2)/(-(cep^2)/Der2o),rate=(epsilon[b-1]-min(0,Der1o/Der2o))/(-(cep^2)/Der2o))

        Der1n=D1(epsilon[b],gam[b-1,])
        Der2n=D2(epsilon[b])

        dn=dgamma(epsilon[b-1],((epsilon[b]-min(0,Der1n/Der2n))^2)/(-(cep^2)/Der2n),rate=(epsilon[b]-min(0,Der1n/Der2n))/(-(cep^2)/Der2n))
        do=dgamma(epsilon[b],((epsilon[b-1]-min(0,Der1o/Der2o))^2)/(-(cep^2)/Der2o),rate=(epsilon[b-1]-min(0,Der1o/Der2o))/(-(cep^2)/Der2o))
        pn=(n*epsilon[b]+psi1-1)*log(epsilon[b])-epsilon[b]*(sum(gam[b-1,])+w)+(epsilon[b]-1)*sum(log(gam[b-1,]))-n*log(gamma(epsilon[b]))
        po=(n*epsilon[b-1]+psi1-1)*log(epsilon[b-1])-epsilon[b-1]*(sum(gam[b-1,])+w)+(epsilon[b-1]-1)*sum(log(gam[b-1,]))-n*log(gamma(epsilon[b-1]))




        alphaep=log(dn)-log(do)+pn-po

        if(is.nan(alphaep)==TRUE){
          epsilon[b]=epsilon[b-1]
          Indepsilon[b]=0
        }else{
          U=log(runif(1,0,1))

          if(U>alphaep){
            epsilon[b]=epsilon[b-1]
            Indepsilon[b]=0
          }else{Indepsilon[b]=1}
        }

        ####Frailty Sampler here
        ####Gam here is not how it's done
        iter[2]="gamma"

        S1=s1[b-1,]
        S1=S1[!is.na(S1)]
        S2=s2[b-1,]
        S2=S2[!is.na(S2)]
        S3=s3[b-1,]
        S3=S3[!is.na(S3)]

        L1=lam1[b-1,]
        L1=as.matrix(L1[!is.na(L1)])
        L2=lam2[b-1,]
        L2=as.matrix(L2[!is.na(L2)])
        L3=lam3[b-1,]
        L3=as.matrix(L3[!is.na(L3)])

        phi1=phifun(Y1,Y1,I1,I2,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),S1,S2,S3,
                    L1,L2,L3,epsilon[b],X)
        ##Sample
        for(i in 1:n){
          gam[b,i]=rgamma(1,1/epsilon[b]+I1[i]+I2[i],rate=phi1[i])
        }


        ############################################
        #####Start LogBH Samplers###################
        ############################################
        ####Lam1####

        iter[1]="LogBH1"
        iter[2]="matrixsetup"

        W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
        Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


        length1=rep(0,J1+1)




        for(j in 1:length(length1)){
          length1[j]=s1[b-1,j+1]-s1[b-1,j]
        }


        if(J1<2){
          if(J1==1){
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            SigLam1=solve(diag(J1+1)-W1)%*%Q1





          }else{

            Q1=as.matrix(2/(m1))
            SigLam1=Q1
          }
        }else{


          for(j in 2:J1){
            W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
          }


          Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
          Q1[1,1]=2/(2*length1[1]+length1[2])
          W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
          W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


          SigLam1=solve(diag(J1+1)-W1)%*%Q1

        }



        iter[2]="Mu"
        ##Lambda1 Hierarchical Sampler
        ##Mulam

        if(J1>0){

          Mulam1[b]=rnorm(1,(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%L1)/(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%as.matrix(rep(1,J1+1))),sqrt(Siglam1[b-1]/(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%as.matrix(rep(1,J1+1)))))


          Siglam1[b]=1/rgamma(1,a1+(J1+1)/2,rate=b1+.5*(t(as.matrix(rep(Mulam1[b],J1+1))-L1)%*%solve(SigLam1)%*%(as.matrix(rep(Mulam1[b],J1+1))-L1)))


          ##Siglam

          iter[2]="Sigma"
        }else{



          Mulam1[b]=rnorm(1,lam1[b-1,1],sqrt(Siglam1[b-1]))


          Siglam1[b]=1/rgamma(1,a1+1/2,rate=b1+.5*(Mulam1[b]-lam1[b-1,1])^2)



        }

        #if(is.finite(Mulam1[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}
        #if(is.finite(Siglam1[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}



        #lambda1
        iter[2]="lam1"
        lam1[b,]=lam1[b-1,]
        #######

        for(m in 1:(J1+1)){



          lam=lam1[b,]
          lam=lam[is.na(lam)==FALSE]
          lambda=lam

          lam[m]=lambda[m]+runif(1,-cl1,cl1)


          if(J1==0){
            do=log(dnorm(lambda[m],Mulam1[b],sqrt(Siglam1[b])))
            dn=log(dnorm(lam[m],Mulam1[b],sqrt(Siglam1[b])))
          }else{


            #do=-(t(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])

            #dn=-(t(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])


            do=dmvnorm(lambda,rep(Mulam1[b],J1+1),Siglam1[b]*SigLam1)
            do=dmvnorm(lam,rep(Mulam1[b],J1+1),Siglam1[b]*SigLam1)
          }

          Likeo=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b-1,],lam3[b-1,],gam[b,])

          Liken=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam,lam2[b-1,],lam3[b-1,],gam[b,])



          U=log(runif(1,0,1))
          alphalam=Liken-Likeo+dn-do

          if(is.nan(alphalam)==TRUE){
            lam1[b,m]=lam1[b-1,m]
            Acceptlam1[b,m]=0
          }else{

            if(U<alphalam){
              Acceptlam1[b,m]=1
              lam1[b,m]=lam[m]
            }else{Acceptlam1[b,m]=0}
          }


        }




        ####Lam2####

        iter[1]="LogBH2"
        iter[2]="matrixsetup"

        W1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)
        Q1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)


        length1=diff(s2[b-1,])







        if(J2<2){
          if(J2==1){
            W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
            SigLam2=solve(diag(J2+1)-W1)%*%Q1





          }else{

            Q1=as.matrix(2/(m2))
            SigLam2=Q1
          }
        }else{


          for(j in 2:J2){
            W1[j,j-1]=(clam2*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            W1[j,j+1]=(clam2*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
          }


          Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
          Q1[1,1]=2/(2*length1[1]+length1[2])
          W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
          W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])


          SigLam2=solve(diag(J2+1)-W1)%*%Q1

        }

        iter[2]="Mu"
        ##Lambda1 Hierarchical Sampler
        ##Mulam
        if(J2>0){

          Mulam2[b]=rnorm(1,(t(as.matrix(rep(1,J2+1)))%*%solve(SigLam2)%*%L2)/(t(as.matrix(rep(1,J2+1)))%*%solve(SigLam2)%*%as.matrix(rep(1,J2+1))),sqrt(Siglam2[b-1]/(t(as.matrix(rep(1,J2+1)))%*%solve(SigLam2)%*%as.matrix(rep(1,J2+1)))))


          Siglam2[b]=1/rgamma(1,a2+(J2+1)/2,rate=b2+.5*(t(as.matrix(rep(Mulam2[b],J2+1))-L2)%*%solve(SigLam2)%*%(as.matrix(rep(Mulam2[b],J2+1))-L2)))


          ##Siglam
          iter[2]="Sigma"

        }else{



          Mulam2[b]=rnorm(1,lam2[b-1,1],sqrt(Siglam2[b-1]))

          Siglam2[b]=1/rgamma(1,a2+1/2,rate=b2+.5*(Mulam2[b]-lam2[b-1,1])^2)




        }

        #if(is.finite(Mulam2[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}
        #if(is.finite(Siglam2[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}
        #lambda1
        iter[2]="lam2"
        lam2[b,]=lam2[b-1,]
        #######


        for(m in 1:(J2+1)){



          lam=lam2[b,]
          lam=lam[is.na(lam)==FALSE]
          lambda=lam

          lam[m]=lambda[m]+runif(1,-cl2,cl2)






          if(J2==0){
            do=log(dnorm(lambda[m],Mulam2[b],sqrt(Siglam2[b])))
            dn=log(dnorm(lam[m],Mulam2[b],sqrt(Siglam2[b])))
          }else{


            #do=-(t(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])

            #dn=-(t(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])


            do=dmvnorm(lambda,rep(Mulam2[b],J2+1),Siglam2[b]*SigLam2)
            do=dmvnorm(lam,rep(Mulam2[b],J2+1),Siglam2[b]*SigLam2)
          }



          #do=-(t(as.matrix(lambda)-as.matrix(rep(Mulam2[b],J2+1)))%*%solve(SigLam2)%*%(as.matrix(lambda)-as.matrix(rep(Mulam2[b],J2+1))))/(2*Siglam2[b])

          #dn=-(t(as.matrix(lam)-as.matrix(rep(Mulam2[b],J2+1)))%*%solve(SigLam2)%*%(as.matrix(lam)-as.matrix(rep(Mulam2[b],J2+1))))/(2*Siglam2[b])


          Likeo=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam3[b-1,],gam[b,])

          Liken=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam,lam3[b-1,],gam[b,])



          U=log(runif(1,0,1))
          alphalam=Liken-Likeo+dn-do

          if(is.nan(alphalam)==TRUE){
            lam2[b,m]=lam2[b-1,m]
            Acceptlam2[b,m]=0
          }else{

            if(U<alphalam){
              Acceptlam2[b,m]=1
              lam2[b,m]=lam[m]
            }else{Acceptlam2[b,m]=0}
          }


        }


        ####Lam2####

        iter[1]="LogBH3"
        iter[2]="matrixsetup"

        W1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)
        Q1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)


        length1=diff(s3[b-1,])







        if(J3<2){
          if(J3==1){
            W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
            SigLam3=solve(diag(J3+1)-W1)%*%Q1





          }else{

            Q1=as.matrix(2/(m3))
            SigLam3=Q1
          }
        }else{


          for(j in 2:J3){
            W1[j,j-1]=(clam3*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            W1[j,j+1]=(clam3*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
          }


          Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
          Q1[1,1]=2/(2*length1[1]+length1[2])
          W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
          W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])


          SigLam3=solve(diag(J3+1)-W1)%*%Q1

        }

        iter[2]="Mu"
        ##Lambda1 Hierarchical Sampler
        ##Mulam

        if(J3>0){

          iter[2]="Sigma"


          Mulam3[b]=rnorm(1,(t(as.matrix(rep(1,J3+1)))%*%solve(SigLam3)%*%L3)/(t(as.matrix(rep(1,J3+1)))%*%solve(SigLam3)%*%as.matrix(rep(1,J3+1))),sqrt(Siglam3[b-1]/(t(as.matrix(rep(1,J3+1)))%*%solve(SigLam3)%*%as.matrix(rep(1,J3+1)))))
          ##Siglam

          Siglam3[b]=1/rgamma(1,a3+(J3+1)/2,rate=b3+.5*(t(as.matrix(rep(Mulam3[b],J3+1))-L3)%*%solve(SigLam3)%*%(as.matrix(rep(Mulam3[b],J3+1))-L3)))


        }else{



          Mulam3[b]=rnorm(1,lam3[b-1,1],sqrt(Siglam3[b-1]))

          Siglam3[b]=1/rgamma(1,a3+1/2,rate=b3+.5*(Mulam3[b]-lam3[b-1,1])^2)



        }

        #if(is.finite(Mulam3[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}
        #if(is.finite(Siglam3[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}

        #lambda3
        iter[2]="lam3"
        lam3[b,]=lam3[b-1,]
        #######

        for(m in 1:(J3+1)){


          lam=lam3[b,]
          lam=lam[is.na(lam)==FALSE]
          lambda=lam

          lam[m]=lambda[m]+runif(1,-cl3,cl3)



          if(J3==0){
            do=log(dnorm(lambda[m],Mulam3[b],sqrt(Siglam3[b])))
            dn=log(dnorm(lam[m],Mulam3[b],sqrt(Siglam3[b])))
          }else{


            #do=-(t(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])

            #dn=-(t(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])


            do=dmvnorm(lambda,rep(Mulam3[b],J3+1),Siglam3[b]*SigLam3)
            do=dmvnorm(lam,rep(Mulam3[b],J3+1),Siglam3[b]*SigLam3)
          }



          Likeo=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])

          Liken=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam,gam[b,])


          U=log(runif(1,0,1))
          alphalam=Liken-Likeo+dn-do

          if(is.nan(alphalam)==TRUE){
            lam3[b,m]=lam3[b-1,m]
            Acceptlam3[b,m]=0
          }else{

            if(U<alphalam){
              Acceptlam3[b,m]=1
              lam3[b,m]=lam[m]
            }else{Acceptlam3[b,m]=0}
          }


        }


        ##############################################
        ######## PUT BACK LAMBDA SAMPLERS HERE!!! ###

        ###Delete these later
        s2[b,]=s2[b-1,]
        s3[b,]=s3[b-1,]

        #####################################################
        ###################################################

        iter[1]="Haz1"
        iter[2]="Birth"

        ###Random Perturbation###
        U1=runif(1,0,1)
        #####

        s=s1[b-1,]
        s=s[!is.na(s)]

        if(length(s)<J1max){
          Birth=runif(1,0,m1)

          s1[b,1:(J1+3)]=sort(c(s,Birth))

          for(k in 2:(J1+2)){
            if(Birth>s1[b-1,k-1] & Birth<s1[b-1,k]){
              Ind=k-1
            }
          }

          lam=rep(0,J1+2)

          if(Ind==1 | Ind==J1+1){
            if(Ind==1){
              lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[(Ind+2):length(lam)]=lam1[b,(Ind+1):(J1+1)]
            }else{
              lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[1:(Ind-1)]=lam1[b,1:(Ind-1)]
            }
          }else{
            lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
            lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
            lam[1:(Ind-1)]=lam1[b,1:(Ind-1)]
            lam[(Ind+2):length(lam)]=lam1[b,(Ind+1):(J1+1)]
          }

          lam=lam[!is.na(lam)]

          lambda=lam1[b,]
          lambda=lambda[!is.na(lambda)]

          Lo=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])
          if(J1>0){
            do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))
          }else{
            do=log(dpois(J1,alpha1))+log(dnorm(lambda,Mulam1[b],Siglam1[b]))
          }

          prior=((2*J1+3)*(2*J1+2)*(Birth-s1[b-1,Ind])*(s1[b-1,Ind+1]-Birth))/((m1^2)*(s1[b-1,Ind+1]-s1[b-1,Ind]))

          G1=G1+1
          J1=J1+1

          Ln=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b-1,],s3[b-1,],lam,lam2[b,],lam3[b,],gam[b,])


          ##Make SigLam1



          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


          length1=diff(s1[b,])


          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1n=solve(diag(J1+1)-W1)%*%Q1





            }else{

              SigLam1n=2/m1
            }
          }else{


            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1n=solve(diag(J1+1)-W1)%*%Q1

          }



          dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),Siglam1[b]*SigLam1n))




          alpha=Ln-Lo+dn-do-log(U1*(1-U1)) + log(prior)

          if(is.nan(alpha)==TRUE){
            IndB1[b]=0
            s1[b,]=s1[b-1,]
            J1=J1-1
            G1=G1-1
          }else{

            U=log(runif(1,0,1))

            if(U<alpha){
              IndB1[b]=1
              lam1[b,1:(J1+1)]=lam
            }else{
              s1[b,]=s1[b-1,]
              IndB1[b]=0
              J1=J1-1
              G1=G1-1
            }

          }


        }else{
          s1[b,]=s1[b-1,]
          IndB1[b]=0
        }


        #########################################################
        ###################Death Sampler#########################
        ##########################################################
        iter[2]="Death"

        U1=runif(1,0,1)

        if(J1==0){
          IndD1[b]=0
          s1[b,]=s1[b-1,]
        }else{

          if(J1==1){
            Ind=2
          }else{

            Ind=sample(2:(J1+1),1)
          }


          s=s1[b,]
          s=s[-Ind]

          lam=lam1[b,]
          lambda=lam[!is.na(lam)]

          lam=lam[!is.na(lam)]
          lam=lam[-Ind]

          lam[Ind-1]=((s1[b,Ind]-s1[b,Ind-1])*lam1[b,Ind-1]+(s1[b,Ind+1]-s1[b,Ind])*lam1[b,Ind])/(s1[b,Ind+1]-s1[b,Ind-1])



          #############################################
          ####Sets up SigLam1 matrix for old density###
          #############################################


          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


          length1=diff(s1[b,])



          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1=solve(diag(J1+1)-W1)%*%Q1

              do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))




            }else{


              do=log(dpois(J1,alpha1))+log(dnorm(lambda,Mulam1[b],Siglam1[b]))

            }
          }else{

            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1=solve(diag(J1+1)-W1)%*%Q1

            do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))


          }
          #############################################
          #############################################

          Lo=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])


          prior=((m1^2)*(s1[b,Ind+1]-s1[b,Ind-1]))/((2*J1+1)*(2*J1)*(s1[b,Ind]-s1[b,Ind-1])*(s1[b,Ind+1]-s1[b,Ind]))


          G1=G1-1
          J1=J1-1


          Ln=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s,s2[b-1,],s3[b-1,],lam,lam2[b,],lam3[b,],gam[b,])

          ###Make siglam matrix



          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)



          length1=rep(0,J1+1)




          for(j in 1:length(length1)){
            length1[j]=s[j+1]-s[j]
          }


          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1n=solve(diag(J1+1)-W1)%*%Q1


              dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),SigLam1n*Siglam1[b]))



            }else{

              SigLam1n=2/m1
              dn=log(dpois(J1,alpha1))+log(dnorm(lam,Mulam1[b],Siglam1[b]))

            }
          }else{


            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1n=solve(diag(J1+1)-W1)%*%Q1

            dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),SigLam1n*Siglam1[b]))


          }
          ####




          alpha=Ln-Lo+dn-do+log(prior)+log(U1*(1-U1))

          if(is.nan(alpha)==TRUE){
            IndD1[b]=0
            J1=J1+1
            G1=G1+1
          }else{

            U=log(runif(1,0,1))

            iter[2]="AcceptRejDeath"

            if(U<alpha){
              s1[b,]=c(s,NA)
              IndD1[b]=1
              lam1[b,1:(J1+1)]=lam
              lam1[b,(J1+2):J1max]=rep(NA,J1max-J1-1)
            }else{
              IndD1[b]=0
              J1=J1+1
              G1=G1+1
            }
          }

          ####End else
        }
        ##






        #######################
        #####End of Death sampler
        ######################


        #####################################################
        ###################################################

        iter[1]="Haz2"
        iter[2]="Birth"

        ###Random Perturbation###
        U2=runif(1,0,1)
        #####

        s=s2[b-1,]
        s=s[!is.na(s)]

        if(length(s)<J2max){
          Birth=runif(1,0,m2)

          s2[b,1:(J2+3)]=sort(c(s,Birth))

          for(k in 2:(J2+2)){
            if(Birth>s2[b-1,k-1] & Birth<s2[b-1,k]){
              Ind=k-1
            }
          }

          lam=rep(0,J2+2)

          if(Ind==1 | Ind==J2+1){
            if(Ind==1){
              lam[Ind]=lam2[b,Ind] - ((s2[b-1,Ind+1]-Birth)/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
              lam[Ind+1]=lam2[b,Ind] + ((Birth-s2[b-1,Ind])/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
              lam[(Ind+2):length(lam)]=lam2[b,(Ind+1):(J2+1)]
            }else{
              lam[Ind]=lam2[b,Ind] - ((s2[b-1,Ind+1]-Birth)/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
              lam[Ind+1]=lam2[b,Ind] + ((Birth-s2[b-1,Ind])/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
              lam[1:(Ind-1)]=lam2[b,1:(Ind-1)]
            }
          }else{
            lam[Ind]=lam2[b,Ind] - ((s2[b-1,Ind+1]-Birth)/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
            lam[Ind+1]=lam2[b,Ind] + ((Birth-s2[b-1,Ind])/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
            lam[1:(Ind-1)]=lam2[b,1:(Ind-1)]
            lam[(Ind+2):length(lam)]=lam2[b,(Ind+1):(J2+1)]
          }

          lam=lam[!is.na(lam)]

          lambda=lam2[b,]
          lambda=lambda[!is.na(lambda)]

          Lo=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])
          if(J2>0){
            do=log(dpois(J2,alpha2))+log(dmvnorm(lambda,rep(Mulam2[b],length(lambda)),SigLam2*Siglam2[b]))
          }else{
            do=log(dpois(J2,alpha2))+log(dnorm(lambda,Mulam2[b],Siglam2[b]))
          }

          prior=((2*J2+3)*(2*J2+2)*(Birth-s2[b-1,Ind])*(s2[b-1,Ind+1]-Birth))/((m2^2)*(s2[b-1,Ind+1]-s2[b-1,Ind]))

          G2=G2+1
          J2=J2+1

          Ln=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s3[b-1,],lam1[b,],lam,lam3[b,],gam[b,])


          ##Make SigLam1



          W1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)
          Q1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)


          length1=diff(s2[b,])




          if(J2<2){
            if(J2==1){
              W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
              SigLam2n=solve(diag(J2+1)-W1)%*%Q1





            }else{

              SigLam2n=2/m2
            }
          }else{


            for(j in 2:J2){
              W1[j,j-1]=(clam2*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam2*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])


            SigLam2n=solve(diag(J2+1)-W1)%*%Q1

          }



          dn=log(dpois(J2,alpha2))+log(dmvnorm(lam,rep(Mulam2[b],length(lam)),Siglam2[b]*SigLam2n))




          alpha=Ln-Lo+dn-do-log(U2*(1-U2)) + log(prior)

          if(is.nan(alpha)==TRUE){
            IndB2[b]=0
            s2[b,]=s2[b-1,]
            J2=J2-1
            G2=G2-1
          }else{

            U=log(runif(1,0,1))

            if(U<alpha){
              IndB2[b]=1
              lam2[b,1:(J2+1)]=lam
            }else{
              s2[b,]=s2[b-1,]
              IndB2[b]=0
              J2=J2-1
              G2=G2-1
            }

          }


        }else{
          s2[b,]=s2[b-1,]
          IndB2[b]=0
        }


        #########################################################
        ###################Death Sampler#########################
        ##########################################################
        iter[2]="Death"

        U2=runif(1,0,1)
        if(J2==0){
          IndD2[b]=0
          s2[b,]=s2[b-1,]
        }else{

          if(J2==1){
            Ind=2
          }else{

            Ind=sample(2:(J2+1),1)
          }


          s=s2[b,]
          s=s[-Ind]

          lam=lam2[b,]
          lambda=lam[!is.na(lam)]

          lam=lam[!is.na(lam)]
          lam=lam[-Ind]

          lam[Ind-1]=((s2[b,Ind]-s2[b,Ind-1])*lam2[b,Ind-1]+(s2[b,Ind+1]-s2[b,Ind])*lam2[b,Ind])/(s2[b,Ind+1]-s2[b,Ind-1])



          #############################################
          ####Sets up SigLam1 matrix for old density###
          #############################################


          W1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)
          Q1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)


          length1=diff(s2[b,])







          if(J2<2){
            if(J2==1){
              W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
              SigLam2=solve(diag(J2+1)-W1)%*%Q1

              do=log(dpois(J2,alpha2))+log(dmvnorm(lambda,rep(Mulam2[b],length(lambda)),SigLam2*Siglam2[b]))




            }else{

              SigLam2=2/m2

              do=log(dpois(J2,alpha2))+log(dnorm(lambda,Mulam2[b],Siglam2[b]))

            }
          }else{

            for(j in 2:J2){
              W1[j,j-1]=(clam2*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam2*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])


            SigLam2=solve(diag(J2+1)-W1)%*%Q1

            do=log(dpois(J2,alpha2))+log(dmvnorm(lambda,rep(Mulam2[b],length(lambda)),SigLam2*Siglam2[b]))


          }
          #############################################
          #############################################

          Lo=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])


          prior=((m2^2)*(s2[b,Ind+1]-s2[b,Ind-1]))/((2*J2+1)*(2*J2)*(s2[b,Ind]-s2[b,Ind-1])*(s2[b,Ind+1]-s2[b,Ind]))


          G2=G2-1
          J2=J2-1


          Ln=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s,s3[b-1,],lam1[b,],lam,lam3[b,],gam[b,])

          ###Make siglam matrix



          W1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)
          Q1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)



          length1=rep(0,J2+1)




          for(j in 1:length(length1)){
            length1[j]=s[j+1]-s[j]
          }


          if(J2<2){
            if(J2==1){
              W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
              SigLam2n=solve(diag(J2+1)-W1)%*%Q1


              dn=log(dpois(J2,alpha2))+log(dmvnorm(lam,rep(Mulam2[b],length(lam)),SigLam2n*Siglam2[b]))



            }else{

              dn=log(dpois(J2,alpha2))+log(dnorm(lam,Mulam2[b],Siglam2[b]))

            }
          }else{


            for(j in 2:J2){
              W1[j,j-1]=(clam2*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam2*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])


            SigLam2n=solve(diag(J2+1)-W1)%*%Q1

            dn=log(dpois(J2,alpha2))+log(dmvnorm(lam,rep(Mulam2[b],length(lam)),SigLam2n*Siglam2[b]))


          }
          ####




          alpha=Ln-Lo+dn-do+log(prior)+log(U2*(1-U2))

          if(is.nan(alpha)==TRUE){
            IndD2[b]=0
            J2=J2+1
            G2=G2+1
          }else{

            U=log(runif(1,0,1))

            iter[2]="AcceptRejDeath"

            if(U<alpha){
              s2[b,]=c(s,NA)
              IndD2[b]=1
              lam2[b,1:(J2+1)]=lam
              lam2[b,(J2+2):J2max]=rep(NA,J2max-J2-1)
            }else{
              IndD2[b]=0
              J2=J2+1
              G2=G2+1
            }
          }

          ####End else
        }
        ##






        #######################
        #####End of Death sampler
        ######################


        #####################################################
        ###################################################

        iter[1]="Haz3"
        iter[2]="Birth"

        ###Random Perturbation###
        U3=runif(1,0,1)
        #####

        s=s3[b-1,]
        s=s[!is.na(s)]

        if(length(s)<J3max){
          Birth=runif(1,0,m3)

          s3[b,1:(J3+3)]=sort(c(s,Birth))

          for(k in 2:(J3+2)){
            if(Birth>s3[b-1,k-1] & Birth<s3[b-1,k]){
              Ind=k-1
            }
          }

          lam=rep(0,J3+2)

          if(Ind==1 | Ind==J3+1){
            if(Ind==1){
              lam[Ind]=lam3[b,Ind] - ((s3[b-1,Ind+1]-Birth)/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
              lam[Ind+1]=lam3[b,Ind] + ((Birth-s3[b-1,Ind])/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
              lam[(Ind+2):length(lam)]=lam3[b,(Ind+1):(J3+1)]
            }else{
              lam[Ind]=lam3[b,Ind] - ((s3[b-1,Ind+1]-Birth)/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
              lam[Ind+1]=lam3[b,Ind] + ((Birth-s3[b-1,Ind])/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
              lam[1:(Ind-1)]=lam3[b,1:(Ind-1)]
            }
          }else{
            lam[Ind]=lam3[b,Ind] - ((s3[b-1,Ind+1]-Birth)/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
            lam[Ind+1]=lam3[b,Ind] + ((Birth-s3[b-1,Ind])/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
            lam[1:(Ind-1)]=lam3[b,1:(Ind-1)]
            lam[(Ind+2):length(lam)]=lam3[b,(Ind+1):(J3+1)]
          }

          lam=lam[!is.na(lam)]

          lambda=lam3[b,]
          lambda=lambda[!is.na(lambda)]

          Lo=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])
          if(J3>0){
            do=log(dpois(J3,alpha3))+log(dmvnorm(lambda,rep(Mulam3[b],length(lambda)),SigLam3*Siglam3[b]))
          }else{
            do=log(dpois(J3,alpha3))+log(dnorm(lambda,Mulam3[b],Siglam3[b]))
          }

          prior=((2*J3+3)*(2*J3+2)*(Birth-s3[b-1,Ind])*(s3[b-1,Ind+1]-Birth))/((m3^2)*(s3[b-1,Ind+1]-s3[b-1,Ind]))

          G3=G3+1
          J3=J3+1

          Ln=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s3[b,],lam1[b,],lam2[b,],lam,gam[b,])


          ##Make SigLam1



          W1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)
          Q1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)


          length1=diff(s3[b,])



          if(J3<2){
            if(J3==1){
              W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
              SigLam3n=solve(diag(J3+1)-W1)%*%Q1





            }else{

              SigLam3n=2/m3
            }
          }else{


            for(j in 2:J3){
              W1[j,j-1]=(clam3*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam3*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])


            SigLam3n=solve(diag(J3+1)-W1)%*%Q1

          }



          dn=log(dpois(J3,alpha3))+log(dmvnorm(lam,rep(Mulam3[b],length(lam)),Siglam3[b]*SigLam3n))




          alpha=Ln-Lo+dn-do-log(U3*(1-U3)) + log(prior)

          if(is.nan(alpha)==TRUE){
            IndB3[b]=0
            s3[b,]=s3[b-1,]
            J3=J3-1
            G3=G3-1
          }else{

            U=log(runif(1,0,1))

            if(U<alpha){
              IndB3[b]=1
              lam3[b,1:(J3+1)]=lam
            }else{
              s3[b,]=s3[b-1,]
              IndB3[b]=0
              J3=J3-1
              G3=G3-1
            }

          }


        }else{
          s3[b,]=s3[b-1,]
          IndB3[b]=0
        }


        #########################################################
        ###################Death Sampler#########################
        ##########################################################
        iter[2]="Death"

        U3=runif(1,0,1)

        if(J3==0){
          IndD3[b]=0
          s3[b,]=s3[b-1,]
        }else{

          if(J3==1){
            Ind=2
          }else{

            Ind=sample(2:(J3+1),1)
          }


          s=s3[b,]
          s=s[-Ind]

          lam=lam3[b,]
          lambda=lam[!is.na(lam)]

          lam=lam[!is.na(lam)]
          lam=lam[-Ind]

          lam[Ind-1]=((s3[b,Ind]-s3[b,Ind-1])*lam3[b,Ind-1]+(s3[b,Ind+1]-s3[b,Ind])*lam3[b,Ind])/(s3[b,Ind+1]-s3[b,Ind-1])



          #############################################
          ####Sets up SigLam1 matrix for old density###
          #############################################


          W1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)
          Q1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)


          length1=diff(s3[b,])



          if(J3<2){
            if(J3==1){
              W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
              SigLam3=solve(diag(J3+1)-W1)%*%Q1

              do=log(dpois(J3,alpha3))+log(dmvnorm(lambda,rep(Mulam3[b],length(lambda)),SigLam3*Siglam3[b]))




            }else{


              do=log(dpois(J3,alpha3))+log(dnorm(lambda,Mulam3[b],Siglam3[b]))

            }
          }else{

            for(j in 2:J3){
              W1[j,j-1]=(clam3*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam3*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])


            SigLam3=solve(diag(J3+1)-W1)%*%Q1

            do=log(dpois(J3,alpha3))+log(dmvnorm(lambda,rep(Mulam3[b],length(lambda)),SigLam3*Siglam3[b]))


          }
          #############################################
          #############################################

          Lo=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s3[b,],lam1[b,],lam2[b,],lam3[b,],gam[b,])


          prior=((m3^2)*(s3[b,Ind+1]-s3[b,Ind-1]))/((2*J3+1)*(2*J3)*(s3[b,Ind]-s3[b,Ind-1])*(s3[b,Ind+1]-s3[b,Ind]))


          G3=G3-1
          J3=J3-1


          Ln=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s,lam1[b,],lam2[b,],lam,gam[b,])

          ###Make siglam matrix



          W1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)
          Q1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)



          length1=rep(0,J3+1)




          for(j in 1:length(length1)){
            length1[j]=s[j+1]-s[j]
          }


          if(J3<2){
            if(J3==1){
              W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
              SigLam3n=solve(diag(J3+1)-W1)%*%Q1


              dn=log(dpois(J3,alpha3))+log(dmvnorm(lam,rep(Mulam3[b],length(lam)),SigLam3n*Siglam3[b]))



            }else{

              dn=log(dpois(J3,alpha3))+log(dnorm(lam,Mulam3[b],Siglam3[b]))

            }
          }else{


            for(j in 2:J3){
              W1[j,j-1]=(clam3*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam3*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])


            SigLam3n=solve(diag(J3+1)-W1)%*%Q1

            dn=log(dpois(J3,alpha3))+log(dmvnorm(lam,rep(Mulam3[b],length(lam)),SigLam3n*Siglam3[b]))


          }
          ####




          alpha=Ln-Lo+dn-do+log(prior)+log(U3*(1-U3))

          if(is.nan(alpha)==TRUE){
            IndD3[b]=0
            J3=J3+1
            G3=G3+1
          }else{

            U=log(runif(1,0,1))

            iter[2]="AcceptRejDeath"

            if(U<alpha){
              s3[b,]=c(s,NA)
              IndD3[b]=1
              lam3[b,1:(J3+1)]=lam
              lam3[b,(J3+2):J3max]=rep(NA,J3max-J3-1)
            }else{
              IndD3[b]=0
              J3=J3+1
              G3=G3+1
            }
          }

          ####End else
        }
        split1[b]=J1
        split2[b]=J2
        split3[b]=J3
        ##
        sum1[b]=sum(eta1[b,])
        sum2[b]=sum(eta2[b,])
        sum3[b]=sum(eta3[b,])


      }





      ################End Samplers
      cat("




          ", "



          ", "


          ", "Posterior Inclusion Probabilities after half Burnin", "


          ", "Hazard 1", "

          ", colMeans(eta1[(B*burn+1):B,])*100, "

          ", "Hazard 2", "

          ", colMeans(eta2[(B*burn+1):B,])*100, "

          ", "Hazard 3", "

          ", colMeans(eta3[(B*burn+1):B,])*100,"

          ", "IndEta",mean(Indeta1[(B*burn+1):B])*100,mean(Indeta2[(B*burn+1):B])*100,mean(Indeta3[(B*burn+1):B])*100,"


          ","IndMix",mean(Indmix1[(B*burn+1):B])*100,mean(Indmix2[(B*burn+1):B])*100,mean(Indmix3[(B*burn+1):B])*100,"


          ", "Included Acceptance", "

          ", "Haz1", "
          ", colMeans(IncCond1[(B*burn+1):B,])*100, "

          ", colMeans(Indcond1[(B*burn+1):B,],na.rm=TRUE)*100,"


          ",  "Haz2", "
          ",colMeans(IncCond2[(B*burn+1):B,])*100, "

          ", colMeans(Indcond2[(B*burn+1):B,],na.rm=TRUE)*100,"



          ", "Haz3",colMeans(IncCond3[(B*burn+1):B,])*100,"



          ", colMeans(Indcond1[(B*burn+1):B,],na.rm=TRUE)*100,"


          ","Survival","

          ","IndDeath",mean(IndD1[(B*burn+1):B])*100,mean(IndD2[(B*burn+1):B])*100,mean(IndD3[(B*burn+1):B])*100,"

          ","IndBirth",mean(IndB1[(B*burn+1):B])*100,mean(IndB2[(B*burn+1):B])*100,mean(IndB3[(B*burn+1):B])*100,"

          ","Lambda","

          ", "Lam1",
          colMeans(Acceptlam1[(B*burn+1):B,],na.rm=TRUE)*100,"

          ","Lam2",
          colMeans(Acceptlam2[(B*burn+1):B,],na.rm=TRUE)*100,"

          ","Lam3",
          colMeans(Acceptlam3[(B*burn+1):B,],na.rm=TRUE)*100,"

          ","Indepsilon",mean(Indepsilon[(B*burn+1):B])*100)







      Path1= paste0(Path,"/IncCond1.txt")

      write.table(IncCond1[(burn*B+1):B,], Path1, sep="\t")

      Path1= paste0(Path,"/IncCond2.txt")

      write.table(IncCond2[(burn*B+1):B,], Path1, sep="\t")

      Path1= paste0(Path,"/IncCond3.txt")

      write.table(IncCond3[(burn*B+1):B,], Path1, sep="\t")






    }




    if(inc==1){
      cat("One Variable Included")

      Ind1s=rep(0,B)
      Ind2s=rep(0,B)
      Ind3s=rep(0,B)

      for(b in 2:B){



        if(b%%10000==0){cat(b, "iterations",date(), "  ")}else{
          if(b%%5000==0){cat(b, " iterations ")}}

        U=runif(1,0,1)


        iter[1]="etabeta1"

        ###eta1,beta1
        eta1[b,]=eta1[b-1,]
        beta1[b,]=beta1[b-1,]

        if(sum(eta1[b-1,])==0|sum(eta1[b-1,])==p1){
          if(sum(eta1[b-1,])==0){
            ###Add Automatically
            iter[2]="Add"
            Ind=sample(1:p1,1)
            eta1[b,Ind]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ####
            V1 = Sigmanew[spot2,spot2]
            V2 = as.matrix(Sigmanew[-spot2,-spot2])
            V12 = as.matrix(Sigmanew[spot2,-spot2])
            thetab=beta1[b-1,c(includednew, (p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot2])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            beta1[b,Ind]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta1[b,Ind],meannew,varnew))
            ######Accept reject###

            Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b-1,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo+dn + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
            U=log(runif(1,0,1))

            if(is.finite(alphab1)==FALSE){
              eta1[b,]=eta1[b-1,]
              beta1[b,]=beta1[b-1,]
              Indeta1[b]=0
            }else{

              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}
            }
          }
          if(sum(eta1[b-1,])==p1){
            ###Delete Automatically
            Ind=sample(1:p1,1)
            iter[2]="delete"
            eta1[b,Ind]=0
            beta1[b,Ind]=0
            includedold=rep(0,p1)
            for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(includedold[k]==Ind){spotold[k]=k}}
            spot1=max(spotold)


            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])

            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta1[b-1,c(includedold,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta1[b-1,Ind],meanold,varold))
            ######Accept reject###
            Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b-1,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo-do + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta1[b,]=eta1[b-1,]
              beta1[b,]=beta1[b-1,]
              Indeta1[b]=0
            }else{

              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}
            }}
        }else{

          U=runif(1,0,1)

          if(U<psi){
            ###Swapper
            includedold=rep(0,p1)
            iter[2]="swap"
            for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            ones=includedold
            zeros=rep(0,p1)
            for(k in 1:p1){if(eta1[b-1,k]==0){zeros[k]=k}}
            zeros=zeros[zeros != 0]
            ###Sample swap indices###
            if(length(ones)==1){
              Indone=ones}else{
                Indone=sample(ones,1)}
            if(length(zeros)==1){Indzero=zeros}else{
              Indzero=sample(zeros,1)}
            ####Change Beta/eta
            eta1[b,Indone]=0
            eta1[b,Indzero]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(Indone==includedold[k]){spotold[k]=k}}
            spot1=max(spotold)
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Indzero==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ###Generate new vector##
            beta1[b,Indone]=0

            ##meannew,varnew##
            V1 = Sigmanew[spot2,spot2]
            V2 = as.matrix(Sigmanew[-spot2,-spot2])
            V12 = as.matrix(Sigmanew[spot2,-spot2])
            thetab=beta1[b-1,c(includednew,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot2])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            beta1[b,Indzero]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta1[b,Indzero],meannew,varnew))
            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta1[b-1,c(includedold,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta1[b-1,Indone],meanold,varold))
            ######Accept reject###
            Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b-1,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo+dn-do
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta1[b,]=eta1[b-1,]
              beta1[b,]=beta1[b-1,]
              Indeta1[b]=0
            }else{
              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}

            }

          }else{
            ###Add/Delete
            Ind=sample(1:p1,1)
            if(eta1[b-1,Ind]==1){
              ##delete##
              iter[2]="delete"
              eta1[b,Ind]=0
              beta1[b,Ind]=0
              includedold=rep(0,p1)
              for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
              includedold=includedold[includedold != 0]
              spotold=rep(0,length(includedold))
              for(k in 1:length(includedold)){if(Ind==includedold[k]){spotold[k]=k}}
              spot1=max(spotold)


              ###Make sigma matrices##
              Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])

              ###Old density###
              V1 = Sigmaold[spot1,spot1]
              V2 = as.matrix(Sigmaold[-spot1,-spot1])
              V12 = as.matrix(Sigmaold[spot1,-spot1])
              thetab=beta1[b-1,c(includedold,(p1+1):(p1+inc))]
              thetano = as.matrix(thetab[-spot1])
              meanold = t(V12)%*%solve(V2)%*%thetano
              varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              do=log(dnorm(beta1[b-1,Ind],meanold,varold))
              ######Accept reject###
              Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b-1,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo-do + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
              U=log(runif(1,0,1))
              if(is.finite(alphab1)==FALSE){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{
                if(U>alphab1){
                  eta1[b,]=eta1[b-1,]
                  beta1[b,]=beta1[b-1,]
                  Indeta1[b]=0
                }else{Indeta1[b]=1}

              }

            }else{
              ###Add###
              eta1[b,Ind]=1

              iter[2]="add"
              includednew=rep(0,p1)
              for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
              includednew=includednew[includednew != 0]
              spotnew=rep(0,length(includednew))
              for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
              spot2=max(spotnew)
              ###Make sigma matrices##
              Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
              ####
              V1 = Sigmanew[spot2,spot2]
              V2 = as.matrix(Sigmanew[-spot2,-spot2])
              V12 = as.matrix(Sigmanew[spot2,-spot2])
              thetab=beta1[b-1,c(includednew,(p1+1):(p1+inc))]
              thetano = as.matrix(thetab[-spot2])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              beta1[b,Ind]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta1[b,Ind],meannew,varnew))
              ######Accept reject###
              Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b-1,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo+dn + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
              U=log(runif(1,0,1))
              if(is.finite(alphab1)==FALSE){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{
                if(U>alphab1){
                  eta1[b,]=eta1[b-1,]
                  beta1[b,]=beta1[b-1,]
                  Indeta1[b]=0
                }else{Indeta1[b]=1}
              }


            }

          }}



        ####ETABETA 2


        iter[1]="etabeta2"

        ###eta1,beta1
        eta2[b,]=eta2[b-1,]
        beta2[b,]=beta2[b-1,]

        if(sum(eta2[b-1,])==0|sum(eta2[b-1,])==p1){
          if(sum(eta2[b-1,])==0){
            ###Add Automatically
            iter[2]="Add"
            Ind=sample(1:p1,1)
            eta2[b,Ind]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta2[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ####
            V1 = Sigmanew[spot2,spot2]
            V2 = as.matrix(Sigmanew[-spot2,-spot2])
            V12 = as.matrix(Sigmanew[spot2,-spot2])
            thetab=beta2[b-1,c(includednew, (p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot2])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            beta2[b,Ind]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta2[b,Ind],meannew,varnew))
            ######Accept reject###

            Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo+dn + log(beta(sum(eta2[b,])+z1a,p1-sum(eta2[b,])+z1b)) - log(beta(sum(eta2[b-1,])+z1a,p1-sum(eta2[b-1,])+z1b))
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta2[b,]=eta2[b-1,]
              beta2[b,]=beta2[b-1,]
              Indeta2[b]=0
            }else{
              if(U>alphab1){
                eta2[b,]=eta2[b-1,]
                beta2[b,]=beta2[b-1,]
                Indeta2[b]=0
              }else{Indeta2[b]=1}
            }}
          if(sum(eta2[b-1,])==p1){
            ###Delete Automatically
            Ind=sample(1:p1,1)
            iter[2]="delete"
            eta2[b,Ind]=0
            beta2[b,Ind]=0
            includedold=rep(0,p1)
            for(k in 1:p1){if(eta2[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(includedold[k]==Ind){spotold[k]=k}}
            spot1=max(spotold)


            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])

            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta2[b-1,c(includedold,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta2[b-1,Ind],meanold,varold))
            ######Accept reject###
            Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo-do + log(beta(sum(eta2[b,])+z1a,p1-sum(eta2[b,])+z1b)) - log(beta(sum(eta2[b-1,])+z1a,p1-sum(eta2[b-1,])+z1b))
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta2[b,]=eta2[b-1,]
              beta2[b,]=beta2[b-1,]
              Indeta2[b]=0
            }else{
              if(U>alphab1){
                eta2[b,]=eta2[b-1,]
                beta2[b,]=beta2[b-1,]
                Indeta2[b]=0
              }else{Indeta2[b]=1}
            }}
        }else{

          U=runif(1,0,1)

          if(U<psi){
            ###Swapper
            includedold=rep(0,p1)
            iter[2]="swap"
            for(k in 1:p1){if(eta2[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            ones=includedold
            zeros=rep(0,p1)
            for(k in 1:p1){if(eta2[b-1,k]==0){zeros[k]=k}}
            zeros=zeros[zeros != 0]
            ###Sample swap indices###
            if(length(ones)==1){
              Indone=ones}else{
                Indone=sample(ones,1)}
            if(length(zeros)==1){Indzero=zeros}else{
              Indzero=sample(zeros,1)}
            ####Change Beta/eta
            eta2[b,Indone]=0
            eta2[b,Indzero]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta2[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(Indone==includedold[k]){spotold[k]=k}}
            spot1=max(spotold)
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Indzero==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ###Generate new vector##
            beta2[b,Indone]=0

            ##meannew,varnew##
            V1 = Sigmanew[spot2,spot2]
            V2 = as.matrix(Sigmanew[-spot2,-spot2])
            V12 = as.matrix(Sigmanew[spot2,-spot2])
            thetab=beta2[b-1,c(includednew,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot2])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            beta2[b,Indzero]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta2[b,Indzero],meannew,varnew))
            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta2[b-1,c(includedold,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta2[b-1,Indone],meanold,varold))
            ######Accept reject###
            Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo+dn-do
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta2[b,]=eta2[b-1,]
              beta2[b,]=beta2[b-1,]
              Indeta2[b]=0
            }else{
              if(U>alphab1){
                eta2[b,]=eta2[b-1,]
                beta2[b,]=beta2[b-1,]
                Indeta2[b]=0
              }else{Indeta2[b]=1}
            }


          }else{
            ###Add/Delete
            Ind=sample(1:p1,1)
            if(eta2[b-1,Ind]==1){
              ##delete##
              iter[2]="delete"
              eta2[b,Ind]=0
              beta2[b,Ind]=0
              includedold=rep(0,p1)
              for(k in 1:p1){if(eta2[b-1,k]==1){includedold[k]=k}}
              includedold=includedold[includedold != 0]
              spotold=rep(0,length(includedold))
              for(k in 1:length(includedold)){if(Ind==includedold[k]){spotold[k]=k}}
              spot1=max(spotold)


              ###Make sigma matrices##
              Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])

              ###Old density###
              V1 = Sigmaold[spot1,spot1]
              V2 = as.matrix(Sigmaold[-spot1,-spot1])
              V12 = as.matrix(Sigmaold[spot1,-spot1])
              thetab=beta2[b-1,c(includedold,(p1+1):(p1+inc))]
              thetano = as.matrix(thetab[-spot1])
              meanold = t(V12)%*%solve(V2)%*%thetano
              varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              do=log(dnorm(beta2[b-1,Ind],meanold,varold))
              ######Accept reject###
              Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo-do + log(beta(sum(eta2[b,])+z1a,p1-sum(eta2[b,])+z1b)) - log(beta(sum(eta2[b-1,])+z1a,p1-sum(eta2[b-1,])+z1b))
              U=log(runif(1,0,1))
              if(is.finite(alphab1)==FALSE){
                eta2[b,]=eta2[b-1,]
                beta2[b,]=beta2[b-1,]
                Indeta2[b]=0
              }else{
                if(U>alphab1){
                  eta2[b,]=eta2[b-1,]
                  beta2[b,]=beta2[b-1,]
                  Indeta2[b]=0
                }else{Indeta2[b]=1}
              }


            }else{
              ###Add###
              eta2[b,Ind]=1

              iter[2]="add"
              includednew=rep(0,p1)
              for(k in 1:p1){if(eta2[b,k]==1){includednew[k]=k}}
              includednew=includednew[includednew != 0]
              spotnew=rep(0,length(includednew))
              for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
              spot2=max(spotnew)
              ###Make sigma matrices##
              Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
              ####
              V1 = Sigmanew[spot2,spot2]
              V2 = as.matrix(Sigmanew[-spot2,-spot2])
              V12 = as.matrix(Sigmanew[spot2,-spot2])
              thetab=beta2[b-1,c(includednew,(p1+1):(p1+inc))]
              thetano = as.matrix(thetab[-spot2])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              beta2[b,Ind]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta2[b,Ind],meannew,varnew))
              ######Accept reject###
              Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo+dn + log(beta(sum(eta2[b,])+z1a,p1-sum(eta2[b,])+z1b)) - log(beta(sum(eta2[b-1,])+z1a,p1-sum(eta2[b-1,])+z1b))
              U=log(runif(1,0,1))
              if(is.finite(alphab1)==FALSE){
                eta2[b,]=eta2[b-1,]
                beta2[b,]=beta2[b-1,]
                Indeta2[b]=0
              }else{
                if(U>alphab1){
                  eta2[b,]=eta2[b-1,]
                  beta2[b,]=beta2[b-1,]
                  Indeta2[b]=0
                }else{Indeta2[b]=1}
              }


            }

          }}



        #####ETA3###


        ####

        iter[1]="etabeta3"

        ###eta1,beta1
        eta3[b,]=eta3[b-1,]
        beta3[b,]=beta3[b-1,]

        if(sum(eta3[b-1,])==0|sum(eta3[b-1,])==p1){
          if(sum(eta3[b-1,])==0){
            ###Add Automatically
            iter[2]="Add"
            Ind=sample(1:p1,1)
            eta3[b,Ind]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta3[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ####
            V1 = Sigmanew[spot2,spot2]
            V2 = as.matrix(Sigmanew[-spot2,-spot2])
            V12 = as.matrix(Sigmanew[spot2,-spot2])
            thetab=beta3[b-1,c(includednew, (p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot2])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            beta3[b,Ind]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta3[b,Ind],meannew,varnew))
            ######Accept reject###

            Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo+dn + log(beta(sum(eta3[b,])+z1a,p1-sum(eta3[b,])+z1b)) - log(beta(sum(eta3[b-1,])+z1a,p1-sum(eta3[b-1,])+z1b))
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta3[b,]=eta3[b-1,]
              beta3[b,]=beta3[b-1,]
              Indeta3[b]=0
            }else{
              if(U>alphab1){
                eta3[b,]=eta3[b-1,]
                beta3[b,]=beta3[b-1,]
                Indeta3[b]=0
              }else{Indeta3[b]=1}
            }}
          if(sum(eta3[b-1,])==p1){
            ###Delete Automatically
            Ind=sample(1:p1,1)
            iter[2]="delete"
            eta3[b,Ind]=0
            beta3[b,Ind]=0
            includedold=rep(0,p1)
            for(k in 1:p1){if(eta3[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(includedold[k]==Ind){spotold[k]=k}}
            spot1=max(spotold)


            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])

            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta3[b-1,c(includedold,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta3[b-1,Ind],meanold,varold))
            ######Accept reject###
            Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo-do + log(beta(sum(eta3[b,])+z1a,p1-sum(eta3[b,])+z1b)) - log(beta(sum(eta3[b-1,])+z1a,p1-sum(eta3[b-1,])+z1b))
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta3[b,]=eta3[b-1,]
              beta3[b,]=beta3[b-1,]
              Indeta3[b]=0
            }else{
              if(U>alphab1){
                eta3[b,]=eta3[b-1,]
                beta3[b,]=beta3[b-1,]
                Indeta3[b]=0
              }else{Indeta3[b]=1}
            }}
        }else{

          U=runif(1,0,1)

          if(U<psi){
            ###Swapper
            includedold=rep(0,p1)
            iter[2]="swap"
            for(k in 1:p1){if(eta3[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            ones=includedold
            zeros=rep(0,p1)
            for(k in 1:p1){if(eta3[b-1,k]==0){zeros[k]=k}}
            zeros=zeros[zeros != 0]
            ###Sample swap indices###
            if(length(ones)==1){
              Indone=ones}else{
                Indone=sample(ones,1)}
            if(length(zeros)==1){Indzero=zeros}else{
              Indzero=sample(zeros,1)}
            ####Change Beta/eta
            eta3[b,Indone]=0
            eta3[b,Indzero]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta3[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(Indone==includedold[k]){spotold[k]=k}}
            spot1=max(spotold)
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Indzero==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ###Generate new vector##
            beta3[b,Indone]=0

            ##meannew,varnew##
            V1 = Sigmanew[spot2,spot2]
            V2 = as.matrix(Sigmanew[-spot2,-spot2])
            V12 = as.matrix(Sigmanew[spot2,-spot2])
            thetab=beta3[b-1,c(includednew,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot2])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            beta3[b,Indzero]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta3[b,Indzero],meannew,varnew))
            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta3[b-1,c(includedold,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta3[b-1,Indone],meanold,varold))
            ######Accept reject###
            Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo+dn-do
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta3[b,]=eta3[b-1,]
              beta3[b,]=beta3[b-1,]
              Indeta3[b]=0
            }else{
              if(U>alphab1){
                eta3[b,]=eta3[b-1,]
                beta3[b,]=beta3[b-1,]
                Indeta3[b]=0
              }else{Indeta3[b]=1}}



          }else{
            ###Add/Delete
            Ind=sample(1:p1,1)
            if(eta3[b-1,Ind]==1){
              ##delete##
              iter[2]="delete"
              eta3[b,Ind]=0
              beta3[b,Ind]=0
              includedold=rep(0,p1)
              for(k in 1:p1){if(eta3[b-1,k]==1){includedold[k]=k}}
              includedold=includedold[includedold != 0]
              spotold=rep(0,length(includedold))
              for(k in 1:length(includedold)){if(Ind==includedold[k]){spotold[k]=k}}
              spot1=max(spotold)


              ###Make sigma matrices##
              Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])

              ###Old density###
              V1 = Sigmaold[spot1,spot1]
              V2 = as.matrix(Sigmaold[-spot1,-spot1])
              V12 = as.matrix(Sigmaold[spot1,-spot1])
              thetab=beta3[b-1,c(includedold,(p1+1):(p1+inc))]
              thetano = as.matrix(thetab[-spot1])
              meanold = t(V12)%*%solve(V2)%*%thetano
              varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              do=log(dnorm(beta3[b-1,Ind],meanold,varold))
              ######Accept reject###
              Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo-do + log(beta(sum(eta3[b,])+z1a,p1-sum(eta3[b,])+z1b)) - log(beta(sum(eta3[b-1,])+z1a,p1-sum(eta3[b-1,])+z1b))
              U=log(runif(1,0,1))
              if(is.finite(alphab1)==FALSE){
                eta3[b,]=eta3[b-1,]
                beta3[b,]=beta3[b-1,]
                Indeta3[b]=0
              }else{
                if(U>alphab1){
                  eta3[b,]=eta3[b-1,]
                  beta3[b,]=beta3[b-1,]
                  Indeta3[b]=0
                }else{Indeta3[b]=1}}



            }else{
              ###Add###
              eta3[b,Ind]=1

              iter[2]="add"
              includednew=rep(0,p1)
              for(k in 1:p1){if(eta3[b,k]==1){includednew[k]=k}}
              includednew=includednew[includednew != 0]
              spotnew=rep(0,length(includednew))
              for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
              spot2=max(spotnew)
              ###Make sigma matrices##
              Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
              ####
              V1 = Sigmanew[spot2,spot2]
              V2 = as.matrix(Sigmanew[-spot2,-spot2])
              V12 = as.matrix(Sigmanew[spot2,-spot2])
              thetab=beta3[b-1,c(includednew,(p1+1):(p1+inc))]
              thetano = as.matrix(thetab[-spot2])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              beta3[b,Ind]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta3[b,Ind],meannew,varnew))
              ######Accept reject###
              Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo+dn + log(beta(sum(eta3[b,])+z1a,p1-sum(eta3[b,])+z1b)) - log(beta(sum(eta3[b-1,])+z1a,p1-sum(eta3[b-1,])+z1b))
              U=log(runif(1,0,1))
              if(is.finite(alphab1)==FALSE){
                eta3[b,]=eta3[b-1,]
                beta3[b,]=beta3[b-1,]
                Indeta3[b]=0
              }else{
                if(U>alphab1){
                  eta3[b,]=eta3[b-1,]
                  beta3[b,]=beta3[b-1,]
                  Indeta3[b]=0
                }else{Indeta3[b]=1}}



            }

          }}





        ###End SVSS


        ###INCLUDED SAMPLERS

        iter[1]="Beta1"
        iter[2]="Included"

        if(sum(eta1[b,])==0){
          ##Sample Included
          Sigmanew= c*solve(t(X[,(p1+1):(p1+inc)])%*%X[,(p1+1):(p1+inc)])
          zeta1n=beta1[b,(p1+1):(p1+inc)]
          meannew=0
          varnew = sqrt(Sigmanew)
          zeta1=rnorm(1,meannew,varnew)
          dn=log(dnorm(zeta1,meannew,varnew))
          ###density old
          do=log(dnorm(zeta1n,meannew,varnew))
          beta=beta1[b,]
          beta[(p1+1):(p1+inc)]=zeta1
          Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
          Liken=LK1(Y1,Y2,I1,I2,X,beta,beta2[b-1,],beta3[b-1,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

          alphab1m=Liken-Likeo+dn -do
          U=log(runif(1,0,1))

          if(is.finite(alphab1m)==FALSE){
            Ind1s[b]=0
          }else{
            if(U>alphab1m){
              Ind1s[b]=0
            }else{Ind1s[b]=1
            beta1[b,]=beta
            zeta1n=zeta1
            }}
          ##End Inc Sampler
        }else{
          includednew=rep(0,p1)
          for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
          includednew=includednew[includednew != 0]
          zeta1n=beta1[b,c(includednew,(p1+1):(p1+inc))]

          ###Make sigma matrices##
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
          ####
          p=length(includednew)+inc
          ####Update All included variables

          for(k in (length(includednew)+1):(length(includednew)+inc)){
            zeta1=zeta1n
            V1 = Sigmanew[k,k]
            V2 = as.matrix(Sigmanew[-k,-k])
            V12 = as.matrix(Sigmanew[k,-k])
            thetano = as.matrix(zeta1[-k])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            zeta1[k]=rnorm(1,meannew,varnew)
            dn=log(dnorm(zeta1[k],meannew,varnew))
            ###density old
            do=log(dnorm(beta1[b,(p1+k-length(includednew))],meannew,varnew))


            ######Accept reject###
            Likeo=LK1(Y1,Y2,I1,I2,X,c(beta1[b,1:p1],zeta1n[(length(zeta1n)-inc+1):length(zeta1n)]),beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK1(Y1,Y2,I1,I2,X,c(beta1[b,1:p1],zeta1[(length(zeta1n)-inc+1):length(zeta1n)]),beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1s=Liken-Likeo+dn -do
            U=log(runif(1,0,1))

            if(is.finite(alphab1s)==FALSE){
              Ins1s[b]=0

            }else{

              if(U>alphab1s){



                Ind1s[b]=0

              }else{Ind1s[b]=1
              zeta1n=zeta1
              beta1[b,]=c(beta1[b,1:p1],zeta1[(length(zeta1)-inc+1):length(zeta1)])

              }
            }

          }

          ###End included sampler###
        }



        #####Conditional Sampler for Included!###


        if(sum(eta1[b,])>0){

          iter[2]="Conditional Inclusion"
          ##Jointly Update nonzero betas
          zeta1=beta1[b,]
          zeta1=zeta1[zeta1!=0]
          zeta1n=zeta1
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])



          ###############
          ####

          for(k in 1:length(includednew)){


            V1 = Sigmanew[k,k]
            V2 = as.matrix(Sigmanew[-k,-k])
            V12 = as.matrix(Sigmanew[k,-k])
            thetab=beta1[b,c(includednew,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-k])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            zeta1n[k]=rnorm(1,meannew,varnew)
            dn=log(dnorm(zeta1n[k],meannew,varnew))
            ###density old
            do=log(dnorm(zeta1[k],meannew,varnew))






            beta=beta1[b,]
            beta[c(includednew,(p1+1):(p1+inc))]=zeta1n


            Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK1(Y1,Y2,I1,I2,X,beta,beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1m=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphab1m)==FALSE){
              Indcond1[b,k]=0
            }else{
              if(U>alphab1m){
                Indcond1[b,includednew[k]]=0
                zeta1n[k]=zeta1[k]
              }else{Indcond1[b,includednew[k]]=1
              beta1[b,]=beta
              zeta1[k]=zeta1n[k]
              }}

          }



          ##Jointly Update nonzero betas
          iter[2]="mixing"
          zeta1n=beta1[b,]
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])


          zeta1n[c(includednew,(p1+1):(p1+inc))]=rmvnorm(1,rep(0,length(includednew)+inc),Sigmanew)

          beta=beta1[b,]
          beta=beta[beta!=0]

          dn=log(dmvnorm(zeta1n[c(includednew,(p1+1):(p1+inc))],rep(0,length(includednew)+inc),Sigmanew))
          ###density old
          do=log(dmvnorm(beta,rep(0,length(includednew)+inc),Sigmanew))

          ######Accept reject###
          Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
          Liken=LK1(Y1,Y2,I1,I2,X,zeta1n,beta2[b-1,],beta3[b-1,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

          alphamix1=Liken-Likeo+dn -do
          U=log(runif(1,0,1))
          if(is.finite(alphamix1)==FALSE){
            Indmix1[b]=0
          }else{
            if(U>alphamix1){

              Indmix1[b]=0
            }else{Indmix1[b]=1
            beta1[b,]=zeta1n
            }}

        }
        iter[1]="Beta2"
        iter[2]="Included"

        if(sum(eta2[b,])==0){
          ##Sample Included
          Sigmanew= c*solve(t(X[,(p1+1):(p1+inc)])%*%X[,(p1+1):(p1+inc)])
          zeta1n=beta2[b,(p1+1):(p1+inc)]

          meannew = 0
          varnew = sqrt(Sigmanew)
          zeta1=rnorm(1,meannew,varnew)
          dn=log(dnorm(zeta1,meannew,varnew))
          ###density old
          do=log(dnorm(zeta1n,meannew,varnew))
          beta=beta2[b,]
          beta[(p1+1):(p1+inc)]=zeta1
          Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
          Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta,beta3[b-1,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

          alphab1m=Liken-Likeo+dn -do
          U=log(runif(1,0,1))
          if(is.finite(alphab1m)==FALSE){
            Ind2s[b]=0
          }else{
            if(U>alphab1m){
              Ind2s[b]=0
            }else{Ind2s[b]=1
            beta2[b,]=beta
            zeta1n=zeta1
            }}
          ##End Inc Sampler
        }else{
          includednew=rep(0,p1)
          for(k in 1:p1){if(eta2[b,k]==1){includednew[k]=k}}
          includednew=includednew[includednew != 0]
          zeta1n=beta2[b,c(includednew,(p1+1):(p1+inc))]

          ###Make sigma matrices##
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
          ####
          p=length(includednew)+inc
          ####Update All included variables

          for(k in (length(includednew)+1):(length(includednew)+inc)){
            zeta1=zeta1n
            V1 = Sigmanew[k,k]
            V2 = as.matrix(Sigmanew[-k,-k])
            V12 = as.matrix(Sigmanew[k,-k])
            thetano = as.matrix(zeta1[-k])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            zeta1[k]=rnorm(1,meannew,varnew)
            dn=log(dnorm(zeta1[k],meannew,varnew))
            ###density old
            do=log(dnorm(beta2[b,(p1+k-length(includednew))],meannew,varnew))



            ######Accept reject###
            Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],c(beta2[b,1:p1],zeta1n[(length(zeta1n)-inc+1):length(zeta1n)]),beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],c(beta2[b,1:p1],zeta1[(length(zeta1n)-inc+1):length(zeta1n)]),beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1s=Liken-Likeo+dn -do
            U=log(runif(1,0,1))

            if(is.finite(alphab1s)==FALSE){
              Ind2s[b]=0

            }else{

              if(U>alphab1s){



                Ind2s[b]=0

              }else{Ind2s[b]=1
              zeta1n=zeta1
              beta2[b,]=c(beta2[b,1:p1],zeta1[(length(zeta1)-inc+1):length(zeta1)])

              }}

          }

          ###End included sampler###
        }



        #####Conditional Sampler for Included!###

        if(sum(eta2[b,])>0){


          iter[2]="Conditional Inclusion"
          ##Jointly Update nonzero betas
          zeta1=beta2[b,]
          zeta1=zeta1[zeta1!=0]
          zeta1n=zeta1
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])



          ###############
          ####

          for(k in 1:length(includednew)){


            V1 = Sigmanew[k,k]
            V2 = as.matrix(Sigmanew[-k,-k])
            V12 = as.matrix(Sigmanew[k,-k])
            thetab=beta1[b,c(includednew,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-k])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            zeta1n[k]=rnorm(1,meannew,varnew)
            dn=log(dnorm(zeta1n[k],meannew,varnew))
            ###density old
            do=log(dnorm(zeta1[k],meannew,varnew))

            beta=beta2[b,]
            beta[c(includednew,(p1+1):(p1+inc))]=zeta1n


            Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta,beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1m=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphab1m)==FALSE){
              Indcond2[b,k]=0
            }else{
              if(U>alphab1m){
                Indcond2[b,includednew[k]]=0
                zeta1n[k]=zeta1[k]
              }else{Indcond2[b,includednew[k]]=1
              beta2[b,]=beta
              zeta1[k]=zeta1n[k]
              }}

          }



          ##Jointly Update nonzero betas
          iter[2]="mixing"
          zeta1n=beta2[b,]
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])


          zeta1n[c(includednew,(p1+1):(p1+inc))]=rmvnorm(1,rep(0,length(includednew)+inc),Sigmanew)

          beta=beta2[b,]
          beta=beta[beta!=0]

          dn=log(dmvnorm(zeta1n[c(includednew,(p1+1):(p1+inc))],rep(0,length(includednew)+inc),Sigmanew))
          ###density old
          do=log(dmvnorm(beta,rep(0,length(includednew)+inc),Sigmanew))

          ######Accept reject###
          Likeo=LK2(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b-1,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
          Liken=LK2(Y1,Y2,I1,I2,X,beta1[b,],zeta1n,beta3[b-1,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

          alphamix1=Liken-Likeo+dn -do
          U=log(runif(1,0,1))
          if(is.finite(alphamix1)==FALSE){
            Indmix2[b]=0
          }else{

            if(U>alphamix1){

              Indmix2[b]=0
            }else{Indmix2[b]=1
            beta2[b,]=zeta1n
            }
          }


        }




        iter[1]="Beta3"
        iter[2]="Included"

        if(sum(eta3[b,])==0){
          ##Sample Included
          Sigmanew= c*solve(t(X[,(p1+1):(p1+inc)])%*%X[,(p1+1):(p1+inc)])
          zeta1n=beta3[b,(p1+1):(p1+inc)]
          zeta1=zeta1n
          meannew=0
          varnew = sqrt(Sigmanew)
          zeta1=rnorm(1,meannew,varnew)
          dn=log(dnorm(zeta1,meannew,varnew))
          ###density old
          do=log(dnorm(zeta1n,meannew,varnew))
          beta=beta3[b,]
          beta[(p1+1):(p1+inc)]=zeta1
          Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
          Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta,
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

          alphab1m=Liken-Likeo+dn -do
          U=log(runif(1,0,1))
          if(is.finite(alphab1m)==FALSE){
            Ind3s[b]=0
          }else{
            if(U>alphab1m){
              Ind3s[b]=0
            }else{Ind3s[b]=1
            beta3[b,]=beta
            zeta1n=zeta1
            }}
          ##End Inc Sampler
        }else{
          includednew=rep(0,p1)
          for(k in 1:p1){if(eta3[b,k]==1){includednew[k]=k}}
          includednew=includednew[includednew != 0]
          zeta1n=beta3[b-1,c(includednew,(p1+1):(p1+inc))]

          ###Make sigma matrices##
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
          ####
          p=length(includednew)+inc
          ####Update All included variables

          for(k in (length(includednew)+1):(length(includednew)+inc)){
            zeta1=zeta1n
            V1 = Sigmanew[k,k]
            V2 = as.matrix(Sigmanew[-k,-k])
            V12 = as.matrix(Sigmanew[k,-k])
            thetano = as.matrix(zeta1[-k])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            zeta1[k]=rnorm(1,meannew,varnew)
            dn=log(dnorm(zeta1[k],meannew,varnew))
            ###density old
            do=log(dnorm(beta3[b-1,(p1+k-length(includednew))],meannew,varnew))



            ######Accept reject###
            Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],c(beta3[b,1:p1],zeta1n[(length(zeta1n)-inc+1):length(zeta1n)]),
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],c(beta3[b,1:p1],zeta1[(length(zeta1n)-inc+1):length(zeta1n)]),
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1s=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphab1s)==FALSE){
              Ind3s[b]=0

            }else{
              if(U>alphab1s){



                Ind3s[b]=0

              }else{Ind3s[b]=1
              zeta1n=zeta1
              beta3[b,]=c(beta3[b,1:p1],zeta1[(length(zeta1)-inc+1):length(zeta1)])

              }}

          }

          ###End included sampler###
        }



        #####Conditional Sampler for Included!###
        if(sum(eta3[b,])>0){


          iter[2]="Conditional Inclusion"
          ##Jointly Update nonzero betas
          zeta1=beta3[b,]
          zeta1=zeta1[zeta1!=0]
          zeta1n=zeta1
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])



          ###############
          ####

          for(k in 1:length(includednew)){


            V1 = Sigmanew[k,k]
            V2 = as.matrix(Sigmanew[-k,-k])
            V12 = as.matrix(Sigmanew[k,-k])
            thetab=beta1[b,c(includednew,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-k])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            zeta1n[k]=rnorm(1,meannew,varnew)
            dn=log(dnorm(zeta1n[k],meannew,varnew))
            ###density old
            do=log(dnorm(zeta1[k],meannew,varnew))

            beta=beta3[b,]
            beta[c(includednew,(p1+1):(p1+inc))]=zeta1n


            Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta,
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1m=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphab1m)==FALSE){
              Indcond3[b,k]=0
            }else{
              if(U>alphab1m){
                Indcond3[b,includednew[k]]=0
                zeta1n[k]=zeta1[k]
              }else{
                Indcond3[b,includednew[k]]=1
                beta3[b,]=beta
                zeta1[k]=zeta1n[k]
              }}

          }



          ##Jointly Update nonzero betas
          iter[2]="mixing"
          zeta1n=beta3[b,]
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])


          zeta1n[c(includednew,(p1+1):(p1+inc))]=rmvnorm(1,rep(0,length(includednew)+inc),Sigmanew)

          beta=beta3[b,]
          beta=beta[beta!=0]

          dn=log(dmvnorm(zeta1n[c(includednew,(p1+1):(p1+inc))],rep(0,length(includednew)+inc),Sigmanew))
          ###density old
          do=log(dmvnorm(beta,rep(0,length(includednew)+inc),Sigmanew))

          ######Accept reject###
          Likeo=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],beta3[b,],
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
          Liken=LK3(Y1,Y2,I1,I2,X,beta1[b,],beta2[b,],zeta1n,
                    s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

          alphamix1=Liken-Likeo+dn -do


          U=log(runif(1,0,1))

          if(is.finite(alphamix1)==FALSE){
            Indmix3[b]=0

          }else{

            if(U>alphamix1){

              Indmix3[b]=0
            }else{Indmix3[b]=1
            beta3[b,]=zeta1n
            }
          }

        }





        ########################
        ###Frailty Samplers#####
        ########################
        ############Epsilon Sampler#####

        iter[1]="frailty"
        iter[2]="hier"
        Der1o=D1(epsilon[b-1],gam[b-1,])
        Der2o=D2(epsilon[b-1])

        epsilon[b]=rgamma(1,((epsilon[b-1]-min(0,Der1o/Der2o))^2)/(-(cep^2)/Der2o),rate=(epsilon[b-1]-min(0,Der1o/Der2o))/(-(cep^2)/Der2o))

        Der1n=D1(epsilon[b],gam[b-1,])
        Der2n=D2(epsilon[b])

        dn=dgamma(epsilon[b-1],((epsilon[b]-min(0,Der1n/Der2n))^2)/(-(cep^2)/Der2n),rate=(epsilon[b]-min(0,Der1n/Der2n))/(-(cep^2)/Der2n))
        do=dgamma(epsilon[b],((epsilon[b-1]-min(0,Der1o/Der2o))^2)/(-(cep^2)/Der2o),rate=(epsilon[b-1]-min(0,Der1o/Der2o))/(-(cep^2)/Der2o))
        pn=(n*epsilon[b]+psi1-1)*log(epsilon[b])-epsilon[b]*(sum(gam[b-1,])+w)+(epsilon[b]-1)*sum(log(gam[b-1,]))-n*log(gamma(epsilon[b]))
        po=(n*epsilon[b-1]+psi1-1)*log(epsilon[b-1])-epsilon[b-1]*(sum(gam[b-1,])+w)+(epsilon[b-1]-1)*sum(log(gam[b-1,]))-n*log(gamma(epsilon[b-1]))




        alphaep=log(dn)-log(do)+pn-po

        if(is.nan(alphaep)==TRUE){
          epsilon[b]=epsilon[b-1]
          Indepsilon[b]=0
        }else{
          U=log(runif(1,0,1))

          if(U>alphaep){
            epsilon[b]=epsilon[b-1]
            Indepsilon[b]=0
          }else{Indepsilon[b]=1}
        }

        ####Frailty Sampler here
        ####Gam here is not how it's done
        iter[2]="gamma"

        S1=s1[b-1,]
        S1=S1[!is.na(S1)]
        S2=s2[b-1,]
        S2=S2[!is.na(S2)]
        S3=s3[b-1,]
        S3=S3[!is.na(S3)]

        L1=lam1[b-1,]
        L1=as.matrix(L1[!is.na(L1)])
        L2=lam2[b-1,]
        L2=as.matrix(L2[!is.na(L2)])
        L3=lam3[b-1,]
        L3=as.matrix(L3[!is.na(L3)])

        phi1=phifun(Y1,Y1,I1,I2,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),S1,S2,S3,
                    L1,L2,L3,epsilon[b],X)
        ##Sample
        for(i in 1:n){
          gam[b,i]=rgamma(1,1/epsilon[b]+I1[i]+I2[i],rate=phi1[i])
        }


        ############################################
        #####Start LogBH Samplers###################
        ############################################
        ####Lam1####

        iter[1]="LogBH1"
        iter[2]="matrixsetup"

        W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
        Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


        length1=rep(0,J1+1)




        for(j in 1:length(length1)){
          length1[j]=s1[b-1,j+1]-s1[b-1,j]
        }


        if(J1<2){
          if(J1==1){
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            SigLam1=solve(diag(J1+1)-W1)%*%Q1





          }else{

            Q1=as.matrix(2/(m1))
            SigLam1=Q1
          }
        }else{


          for(j in 2:J1){
            W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
          }


          Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
          Q1[1,1]=2/(2*length1[1]+length1[2])
          W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
          W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


          SigLam1=solve(diag(J1+1)-W1)%*%Q1

        }



        iter[2]="Mu"
        ##Lambda1 Hierarchical Sampler
        ##Mulam

        if(J1>0){

          Mulam1[b]=rnorm(1,(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%L1)/(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%as.matrix(rep(1,J1+1))),sqrt(Siglam1[b-1]/(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%as.matrix(rep(1,J1+1)))))


          Siglam1[b]=1/rgamma(1,a1+(J1+1)/2,rate=b1+.5*(t(as.matrix(rep(Mulam1[b],J1+1))-L1)%*%solve(SigLam1)%*%(as.matrix(rep(Mulam1[b],J1+1))-L1)))


          ##Siglam

          iter[2]="Sigma"
        }else{



          Mulam1[b]=rnorm(1,lam1[b-1,1],sqrt(Siglam1[b-1]))


          Siglam1[b]=1/rgamma(1,a1+1/2,rate=b1+.5*(Mulam1[b]-lam1[b-1,1])^2)



        }

        #if(is.finite(Mulam1[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}
        #if(is.finite(Siglam1[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}



        #lambda1
        iter[2]="lam1"
        lam1[b,]=lam1[b-1,]
        #######

        for(m in 1:(J1+1)){



          lam=lam1[b,]
          lam=lam[is.na(lam)==FALSE]
          lambda=lam

          lam[m]=lambda[m]+runif(1,-cl1,cl1)


          if(J1==0){
            do=log(dnorm(lambda[m],Mulam1[b],sqrt(Siglam1[b])))
            dn=log(dnorm(lam[m],Mulam1[b],sqrt(Siglam1[b])))
          }else{


            #do=-(t(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])

            #dn=-(t(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])


            do=dmvnorm(lambda,rep(Mulam1[b],J1+1),Siglam1[b]*SigLam1)
            do=dmvnorm(lam,rep(Mulam1[b],J1+1),Siglam1[b]*SigLam1)
          }

          Likeo=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b-1,],lam3[b-1,],gam[b,])

          Liken=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam,lam2[b-1,],lam3[b-1,],gam[b,])



          U=log(runif(1,0,1))
          alphalam=Liken-Likeo+dn-do

          if(is.nan(alphalam)==TRUE){
            lam1[b,m]=lam1[b-1,m]
            Acceptlam1[b,m]=0
          }else{

            if(U<alphalam){
              Acceptlam1[b,m]=1
              lam1[b,m]=lam[m]
            }else{Acceptlam1[b,m]=0}
          }


        }




        ####Lam2####

        iter[1]="LogBH2"
        iter[2]="matrixsetup"

        W1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)
        Q1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)


        length1=diff(s2[b-1,])







        if(J2<2){
          if(J2==1){
            W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
            SigLam2=solve(diag(J2+1)-W1)%*%Q1





          }else{

            Q1=as.matrix(2/(m2))
            SigLam2=Q1
          }
        }else{


          for(j in 2:J2){
            W1[j,j-1]=(clam2*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            W1[j,j+1]=(clam2*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
          }


          Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
          Q1[1,1]=2/(2*length1[1]+length1[2])
          W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
          W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])


          SigLam2=solve(diag(J2+1)-W1)%*%Q1

        }

        iter[2]="Mu"
        ##Lambda1 Hierarchical Sampler
        ##Mulam
        if(J2>0){

          Mulam2[b]=rnorm(1,(t(as.matrix(rep(1,J2+1)))%*%solve(SigLam2)%*%L2)/(t(as.matrix(rep(1,J2+1)))%*%solve(SigLam2)%*%as.matrix(rep(1,J2+1))),sqrt(Siglam2[b-1]/(t(as.matrix(rep(1,J2+1)))%*%solve(SigLam2)%*%as.matrix(rep(1,J2+1)))))


          Siglam2[b]=1/rgamma(1,a2+(J2+1)/2,rate=b2+.5*(t(as.matrix(rep(Mulam2[b],J2+1))-L2)%*%solve(SigLam2)%*%(as.matrix(rep(Mulam2[b],J2+1))-L2)))


          ##Siglam
          iter[2]="Sigma"

        }else{



          Mulam2[b]=rnorm(1,lam2[b-1,1],sqrt(Siglam2[b-1]))

          Siglam2[b]=1/rgamma(1,a2+1/2,rate=b2+.5*(Mulam2[b]-lam2[b-1,1])^2)




        }

        #if(is.finite(Mulam2[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}
        #if(is.finite(Siglam2[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}
        #lambda1
        iter[2]="lam2"
        lam2[b,]=lam2[b-1,]
        #######


        for(m in 1:(J2+1)){



          lam=lam2[b,]
          lam=lam[is.na(lam)==FALSE]
          lambda=lam

          lam[m]=lambda[m]+runif(1,-cl2,cl2)






          if(J2==0){
            do=log(dnorm(lambda[m],Mulam2[b],sqrt(Siglam2[b])))
            dn=log(dnorm(lam[m],Mulam2[b],sqrt(Siglam2[b])))
          }else{


            #do=-(t(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])

            #dn=-(t(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])


            do=dmvnorm(lambda,rep(Mulam2[b],J2+1),Siglam2[b]*SigLam2)
            do=dmvnorm(lam,rep(Mulam2[b],J2+1),Siglam2[b]*SigLam2)
          }



          #do=-(t(as.matrix(lambda)-as.matrix(rep(Mulam2[b],J2+1)))%*%solve(SigLam2)%*%(as.matrix(lambda)-as.matrix(rep(Mulam2[b],J2+1))))/(2*Siglam2[b])

          #dn=-(t(as.matrix(lam)-as.matrix(rep(Mulam2[b],J2+1)))%*%solve(SigLam2)%*%(as.matrix(lam)-as.matrix(rep(Mulam2[b],J2+1))))/(2*Siglam2[b])


          Likeo=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam3[b-1,],gam[b,])

          Liken=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam,lam3[b-1,],gam[b,])



          U=log(runif(1,0,1))
          alphalam=Liken-Likeo+dn-do

          if(is.nan(alphalam)==TRUE){
            lam2[b,m]=lam2[b-1,m]
            Acceptlam2[b,m]=0
          }else{

            if(U<alphalam){
              Acceptlam2[b,m]=1
              lam2[b,m]=lam[m]
            }else{Acceptlam2[b,m]=0}
          }


        }


        ####Lam2####

        iter[1]="LogBH3"
        iter[2]="matrixsetup"

        W1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)
        Q1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)


        length1=diff(s3[b-1,])







        if(J3<2){
          if(J3==1){
            W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
            SigLam3=solve(diag(J3+1)-W1)%*%Q1





          }else{

            Q1=as.matrix(2/(m3))
            SigLam3=Q1
          }
        }else{


          for(j in 2:J3){
            W1[j,j-1]=(clam3*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            W1[j,j+1]=(clam3*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
          }


          Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
          Q1[1,1]=2/(2*length1[1]+length1[2])
          W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
          W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])


          SigLam3=solve(diag(J3+1)-W1)%*%Q1

        }

        iter[2]="Mu"
        ##Lambda1 Hierarchical Sampler
        ##Mulam

        if(J3>0){

          iter[2]="Sigma"


          Mulam3[b]=rnorm(1,(t(as.matrix(rep(1,J3+1)))%*%solve(SigLam3)%*%L3)/(t(as.matrix(rep(1,J3+1)))%*%solve(SigLam3)%*%as.matrix(rep(1,J3+1))),sqrt(Siglam3[b-1]/(t(as.matrix(rep(1,J3+1)))%*%solve(SigLam3)%*%as.matrix(rep(1,J3+1)))))
          ##Siglam

          Siglam3[b]=1/rgamma(1,a3+(J3+1)/2,rate=b3+.5*(t(as.matrix(rep(Mulam3[b],J3+1))-L3)%*%solve(SigLam3)%*%(as.matrix(rep(Mulam3[b],J3+1))-L3)))


        }else{



          Mulam3[b]=rnorm(1,lam3[b-1,1],sqrt(Siglam3[b-1]))

          Siglam3[b]=1/rgamma(1,a3+1/2,rate=b3+.5*(Mulam3[b]-lam3[b-1,1])^2)



        }

        #if(is.finite(Mulam3[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}
        #if(is.finite(Siglam3[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}

        #lambda3
        iter[2]="lam3"
        lam3[b,]=lam3[b-1,]
        #######

        for(m in 1:(J3+1)){


          lam=lam3[b,]
          lam=lam[is.na(lam)==FALSE]
          lambda=lam

          lam[m]=lambda[m]+runif(1,-cl3,cl3)



          if(J3==0){
            do=log(dnorm(lambda[m],Mulam3[b],sqrt(Siglam3[b])))
            dn=log(dnorm(lam[m],Mulam3[b],sqrt(Siglam3[b])))
          }else{


            #do=-(t(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])

            #dn=-(t(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])


            do=dmvnorm(lambda,rep(Mulam3[b],J3+1),Siglam3[b]*SigLam3)
            do=dmvnorm(lam,rep(Mulam3[b],J3+1),Siglam3[b]*SigLam3)
          }



          Likeo=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])

          Liken=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam,gam[b,])


          U=log(runif(1,0,1))
          alphalam=Liken-Likeo+dn-do

          if(is.nan(alphalam)==TRUE){
            lam3[b,m]=lam3[b-1,m]
            Acceptlam3[b,m]=0
          }else{

            if(U<alphalam){
              Acceptlam3[b,m]=1
              lam3[b,m]=lam[m]
            }else{Acceptlam3[b,m]=0}
          }


        }


        ##############################################
        ######## PUT BACK LAMBDA SAMPLERS HERE!!! ###

        ###Delete these later
        s2[b,]=s2[b-1,]
        s3[b,]=s3[b-1,]

        #####################################################
        ###################################################

        iter[1]="Haz1"
        iter[2]="Birth"

        ###Random Perturbation###
        U1=runif(1,0,1)
        #####

        s=s1[b-1,]
        s=s[!is.na(s)]

        if(length(s)<J1max){
          Birth=runif(1,0,m1)

          s1[b,1:(J1+3)]=sort(c(s,Birth))

          for(k in 2:(J1+2)){
            if(Birth>s1[b-1,k-1] & Birth<s1[b-1,k]){
              Ind=k-1
            }
          }

          lam=rep(0,J1+2)

          if(Ind==1 | Ind==J1+1){
            if(Ind==1){
              lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[(Ind+2):length(lam)]=lam1[b,(Ind+1):(J1+1)]
            }else{
              lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[1:(Ind-1)]=lam1[b,1:(Ind-1)]
            }
          }else{
            lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
            lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
            lam[1:(Ind-1)]=lam1[b,1:(Ind-1)]
            lam[(Ind+2):length(lam)]=lam1[b,(Ind+1):(J1+1)]
          }

          lam=lam[!is.na(lam)]

          lambda=lam1[b,]
          lambda=lambda[!is.na(lambda)]

          Lo=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])
          if(J1>0){
            do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))
          }else{
            do=log(dpois(J1,alpha1))+log(dnorm(lambda,Mulam1[b],Siglam1[b]))
          }

          prior=((2*J1+3)*(2*J1+2)*(Birth-s1[b-1,Ind])*(s1[b-1,Ind+1]-Birth))/((m1^2)*(s1[b-1,Ind+1]-s1[b-1,Ind]))

          G1=G1+1
          J1=J1+1

          Ln=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b-1,],s3[b-1,],lam,lam2[b,],lam3[b,],gam[b,])


          ##Make SigLam1



          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


          length1=diff(s1[b,])


          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1n=solve(diag(J1+1)-W1)%*%Q1





            }else{

              SigLam1n=2/m1
            }
          }else{


            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1n=solve(diag(J1+1)-W1)%*%Q1

          }



          dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),Siglam1[b]*SigLam1n))




          alpha=Ln-Lo+dn-do-log(U1*(1-U1)) + log(prior)

          if(is.nan(alpha)==TRUE){
            IndB1[b]=0
            s1[b,]=s1[b-1,]
            J1=J1-1
            G1=G1-1
          }else{

            U=log(runif(1,0,1))

            if(U<alpha){
              IndB1[b]=1
              lam1[b,1:(J1+1)]=lam
            }else{
              s1[b,]=s1[b-1,]
              IndB1[b]=0
              J1=J1-1
              G1=G1-1
            }

          }


        }else{
          s1[b,]=s1[b-1,]
          IndB1[b]=0
        }


        #########################################################
        ###################Death Sampler#########################
        ##########################################################
        iter[2]="Death"

        U1=runif(1,0,1)

        if(J1==0){
          IndD1[b]=0
          s1[b,]=s1[b-1,]
        }else{

          if(J1==1){
            Ind=2
          }else{

            Ind=sample(2:(J1+1),1)
          }


          s=s1[b,]
          s=s[-Ind]

          lam=lam1[b,]
          lambda=lam[!is.na(lam)]

          lam=lam[!is.na(lam)]
          lam=lam[-Ind]

          lam[Ind-1]=((s1[b,Ind]-s1[b,Ind-1])*lam1[b,Ind-1]+(s1[b,Ind+1]-s1[b,Ind])*lam1[b,Ind])/(s1[b,Ind+1]-s1[b,Ind-1])



          #############################################
          ####Sets up SigLam1 matrix for old density###
          #############################################


          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


          length1=diff(s1[b,])



          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1=solve(diag(J1+1)-W1)%*%Q1

              do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))




            }else{


              do=log(dpois(J1,alpha1))+log(dnorm(lambda,Mulam1[b],Siglam1[b]))

            }
          }else{

            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1=solve(diag(J1+1)-W1)%*%Q1

            do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))


          }
          #############################################
          #############################################

          Lo=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])


          prior=((m1^2)*(s1[b,Ind+1]-s1[b,Ind-1]))/((2*J1+1)*(2*J1)*(s1[b,Ind]-s1[b,Ind-1])*(s1[b,Ind+1]-s1[b,Ind]))


          G1=G1-1
          J1=J1-1


          Ln=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s,s2[b-1,],s3[b-1,],lam,lam2[b,],lam3[b,],gam[b,])

          ###Make siglam matrix



          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)



          length1=rep(0,J1+1)




          for(j in 1:length(length1)){
            length1[j]=s[j+1]-s[j]
          }


          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1n=solve(diag(J1+1)-W1)%*%Q1


              dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),SigLam1n*Siglam1[b]))



            }else{

              SigLam1n=2/m1
              dn=log(dpois(J1,alpha1))+log(dnorm(lam,Mulam1[b],Siglam1[b]))

            }
          }else{


            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1n=solve(diag(J1+1)-W1)%*%Q1

            dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),SigLam1n*Siglam1[b]))


          }
          ####




          alpha=Ln-Lo+dn-do+log(prior)+log(U1*(1-U1))

          if(is.nan(alpha)==TRUE){
            IndD1[b]=0
            J1=J1+1
            G1=G1+1
          }else{

            U=log(runif(1,0,1))

            iter[2]="AcceptRejDeath"

            if(U<alpha){
              s1[b,]=c(s,NA)
              IndD1[b]=1
              lam1[b,1:(J1+1)]=lam
              lam1[b,(J1+2):J1max]=rep(NA,J1max-J1-1)
            }else{
              IndD1[b]=0
              J1=J1+1
              G1=G1+1
            }
          }

          ####End else
        }
        ##






        #######################
        #####End of Death sampler
        ######################


        #####################################################
        ###################################################

        iter[1]="Haz2"
        iter[2]="Birth"

        ###Random Perturbation###
        U2=runif(1,0,1)
        #####

        s=s2[b-1,]
        s=s[!is.na(s)]

        if(length(s)<J2max){
          Birth=runif(1,0,m2)

          s2[b,1:(J2+3)]=sort(c(s,Birth))

          for(k in 2:(J2+2)){
            if(Birth>s2[b-1,k-1] & Birth<s2[b-1,k]){
              Ind=k-1
            }
          }

          lam=rep(0,J2+2)

          if(Ind==1 | Ind==J2+1){
            if(Ind==1){
              lam[Ind]=lam2[b,Ind] - ((s2[b-1,Ind+1]-Birth)/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
              lam[Ind+1]=lam2[b,Ind] + ((Birth-s2[b-1,Ind])/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
              lam[(Ind+2):length(lam)]=lam2[b,(Ind+1):(J2+1)]
            }else{
              lam[Ind]=lam2[b,Ind] - ((s2[b-1,Ind+1]-Birth)/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
              lam[Ind+1]=lam2[b,Ind] + ((Birth-s2[b-1,Ind])/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
              lam[1:(Ind-1)]=lam2[b,1:(Ind-1)]
            }
          }else{
            lam[Ind]=lam2[b,Ind] - ((s2[b-1,Ind+1]-Birth)/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
            lam[Ind+1]=lam2[b,Ind] + ((Birth-s2[b-1,Ind])/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
            lam[1:(Ind-1)]=lam2[b,1:(Ind-1)]
            lam[(Ind+2):length(lam)]=lam2[b,(Ind+1):(J2+1)]
          }

          lam=lam[!is.na(lam)]

          lambda=lam2[b,]
          lambda=lambda[!is.na(lambda)]

          Lo=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])
          if(J2>0){
            do=log(dpois(J2,alpha2))+log(dmvnorm(lambda,rep(Mulam2[b],length(lambda)),SigLam2*Siglam2[b]))
          }else{
            do=log(dpois(J2,alpha2))+log(dnorm(lambda,Mulam2[b],Siglam2[b]))
          }

          prior=((2*J2+3)*(2*J2+2)*(Birth-s2[b-1,Ind])*(s2[b-1,Ind+1]-Birth))/((m2^2)*(s2[b-1,Ind+1]-s2[b-1,Ind]))

          G2=G2+1
          J2=J2+1

          Ln=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s3[b-1,],lam1[b,],lam,lam3[b,],gam[b,])


          ##Make SigLam1



          W1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)
          Q1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)


          length1=diff(s2[b,])




          if(J2<2){
            if(J2==1){
              W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
              SigLam2n=solve(diag(J2+1)-W1)%*%Q1





            }else{

              SigLam2n=2/m2
            }
          }else{


            for(j in 2:J2){
              W1[j,j-1]=(clam2*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam2*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])


            SigLam2n=solve(diag(J2+1)-W1)%*%Q1

          }



          dn=log(dpois(J2,alpha2))+log(dmvnorm(lam,rep(Mulam2[b],length(lam)),Siglam2[b]*SigLam2n))




          alpha=Ln-Lo+dn-do-log(U2*(1-U2)) + log(prior)

          if(is.nan(alpha)==TRUE){
            IndB2[b]=0
            s2[b,]=s2[b-1,]
            J2=J2-1
            G2=G2-1
          }else{

            U=log(runif(1,0,1))

            if(U<alpha){
              IndB2[b]=1
              lam2[b,1:(J2+1)]=lam
            }else{
              s2[b,]=s2[b-1,]
              IndB2[b]=0
              J2=J2-1
              G2=G2-1
            }

          }


        }else{
          s2[b,]=s2[b-1,]
          IndB2[b]=0
        }


        #########################################################
        ###################Death Sampler#########################
        ##########################################################
        iter[2]="Death"

        U2=runif(1,0,1)
        if(J2==0){
          IndD2[b]=0
          s2[b,]=s2[b-1,]
        }else{

          if(J2==1){
            Ind=2
          }else{

            Ind=sample(2:(J2+1),1)
          }


          s=s2[b,]
          s=s[-Ind]

          lam=lam2[b,]
          lambda=lam[!is.na(lam)]

          lam=lam[!is.na(lam)]
          lam=lam[-Ind]

          lam[Ind-1]=((s2[b,Ind]-s2[b,Ind-1])*lam2[b,Ind-1]+(s2[b,Ind+1]-s2[b,Ind])*lam2[b,Ind])/(s2[b,Ind+1]-s2[b,Ind-1])



          #############################################
          ####Sets up SigLam1 matrix for old density###
          #############################################


          W1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)
          Q1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)


          length1=diff(s2[b,])







          if(J2<2){
            if(J2==1){
              W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
              SigLam2=solve(diag(J2+1)-W1)%*%Q1

              do=log(dpois(J2,alpha2))+log(dmvnorm(lambda,rep(Mulam2[b],length(lambda)),SigLam2*Siglam2[b]))




            }else{

              SigLam2=2/m2

              do=log(dpois(J2,alpha2))+log(dnorm(lambda,Mulam2[b],Siglam2[b]))

            }
          }else{

            for(j in 2:J2){
              W1[j,j-1]=(clam2*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam2*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])


            SigLam2=solve(diag(J2+1)-W1)%*%Q1

            do=log(dpois(J2,alpha2))+log(dmvnorm(lambda,rep(Mulam2[b],length(lambda)),SigLam2*Siglam2[b]))


          }
          #############################################
          #############################################

          Lo=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])


          prior=((m2^2)*(s2[b,Ind+1]-s2[b,Ind-1]))/((2*J2+1)*(2*J2)*(s2[b,Ind]-s2[b,Ind-1])*(s2[b,Ind+1]-s2[b,Ind]))


          G2=G2-1
          J2=J2-1


          Ln=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s,s3[b-1,],lam1[b,],lam,lam3[b,],gam[b,])

          ###Make siglam matrix



          W1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)
          Q1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)



          length1=rep(0,J2+1)




          for(j in 1:length(length1)){
            length1[j]=s[j+1]-s[j]
          }


          if(J2<2){
            if(J2==1){
              W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
              SigLam2n=solve(diag(J2+1)-W1)%*%Q1


              dn=log(dpois(J2,alpha2))+log(dmvnorm(lam,rep(Mulam2[b],length(lam)),SigLam2n*Siglam2[b]))



            }else{

              dn=log(dpois(J2,alpha2))+log(dnorm(lam,Mulam2[b],Siglam2[b]))

            }
          }else{


            for(j in 2:J2){
              W1[j,j-1]=(clam2*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam2*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])


            SigLam2n=solve(diag(J2+1)-W1)%*%Q1

            dn=log(dpois(J2,alpha2))+log(dmvnorm(lam,rep(Mulam2[b],length(lam)),SigLam2n*Siglam2[b]))


          }
          ####




          alpha=Ln-Lo+dn-do+log(prior)+log(U2*(1-U2))

          if(is.nan(alpha)==TRUE){
            IndD2[b]=0
            J2=J2+1
            G2=G2+1
          }else{

            U=log(runif(1,0,1))

            iter[2]="AcceptRejDeath"

            if(U<alpha){
              s2[b,]=c(s,NA)
              IndD2[b]=1
              lam2[b,1:(J2+1)]=lam
              lam2[b,(J2+2):J2max]=rep(NA,J2max-J2-1)
            }else{
              IndD2[b]=0
              J2=J2+1
              G2=G2+1
            }
          }

          ####End else
        }
        ##






        #######################
        #####End of Death sampler
        ######################


        #####################################################
        ###################################################

        iter[1]="Haz3"
        iter[2]="Birth"

        ###Random Perturbation###
        U3=runif(1,0,1)
        #####

        s=s3[b-1,]
        s=s[!is.na(s)]

        if(length(s)<J3max){
          Birth=runif(1,0,m3)

          s3[b,1:(J3+3)]=sort(c(s,Birth))

          for(k in 2:(J3+2)){
            if(Birth>s3[b-1,k-1] & Birth<s3[b-1,k]){
              Ind=k-1
            }
          }

          lam=rep(0,J3+2)

          if(Ind==1 | Ind==J3+1){
            if(Ind==1){
              lam[Ind]=lam3[b,Ind] - ((s3[b-1,Ind+1]-Birth)/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
              lam[Ind+1]=lam3[b,Ind] + ((Birth-s3[b-1,Ind])/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
              lam[(Ind+2):length(lam)]=lam3[b,(Ind+1):(J3+1)]
            }else{
              lam[Ind]=lam3[b,Ind] - ((s3[b-1,Ind+1]-Birth)/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
              lam[Ind+1]=lam3[b,Ind] + ((Birth-s3[b-1,Ind])/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
              lam[1:(Ind-1)]=lam3[b,1:(Ind-1)]
            }
          }else{
            lam[Ind]=lam3[b,Ind] - ((s3[b-1,Ind+1]-Birth)/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
            lam[Ind+1]=lam3[b,Ind] + ((Birth-s3[b-1,Ind])/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
            lam[1:(Ind-1)]=lam3[b,1:(Ind-1)]
            lam[(Ind+2):length(lam)]=lam3[b,(Ind+1):(J3+1)]
          }

          lam=lam[!is.na(lam)]

          lambda=lam3[b,]
          lambda=lambda[!is.na(lambda)]

          Lo=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])
          if(J3>0){
            do=log(dpois(J3,alpha3))+log(dmvnorm(lambda,rep(Mulam3[b],length(lambda)),SigLam3*Siglam3[b]))
          }else{
            do=log(dpois(J3,alpha3))+log(dnorm(lambda,Mulam3[b],Siglam3[b]))
          }

          prior=((2*J3+3)*(2*J3+2)*(Birth-s3[b-1,Ind])*(s3[b-1,Ind+1]-Birth))/((m3^2)*(s3[b-1,Ind+1]-s3[b-1,Ind]))

          G3=G3+1
          J3=J3+1

          Ln=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s3[b,],lam1[b,],lam2[b,],lam,gam[b,])


          ##Make SigLam1



          W1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)
          Q1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)


          length1=diff(s3[b,])



          if(J3<2){
            if(J3==1){
              W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
              SigLam3n=solve(diag(J3+1)-W1)%*%Q1





            }else{

              SigLam3n=2/m3
            }
          }else{


            for(j in 2:J3){
              W1[j,j-1]=(clam3*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam3*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])


            SigLam3n=solve(diag(J3+1)-W1)%*%Q1

          }



          dn=log(dpois(J3,alpha3))+log(dmvnorm(lam,rep(Mulam3[b],length(lam)),Siglam3[b]*SigLam3n))




          alpha=Ln-Lo+dn-do-log(U3*(1-U3)) + log(prior)

          if(is.nan(alpha)==TRUE){
            IndB3[b]=0
            s3[b,]=s3[b-1,]
            J3=J3-1
            G3=G3-1
          }else{

            U=log(runif(1,0,1))

            if(U<alpha){
              IndB3[b]=1
              lam3[b,1:(J3+1)]=lam
            }else{
              s3[b,]=s3[b-1,]
              IndB3[b]=0
              J3=J3-1
              G3=G3-1
            }

          }


        }else{
          s3[b,]=s3[b-1,]
          IndB3[b]=0
        }


        #########################################################
        ###################Death Sampler#########################
        ##########################################################
        iter[2]="Death"

        U3=runif(1,0,1)

        if(J3==0){
          IndD3[b]=0
          s3[b,]=s3[b-1,]
        }else{

          if(J3==1){
            Ind=2
          }else{

            Ind=sample(2:(J3+1),1)
          }


          s=s3[b,]
          s=s[-Ind]

          lam=lam3[b,]
          lambda=lam[!is.na(lam)]

          lam=lam[!is.na(lam)]
          lam=lam[-Ind]

          lam[Ind-1]=((s3[b,Ind]-s3[b,Ind-1])*lam3[b,Ind-1]+(s3[b,Ind+1]-s3[b,Ind])*lam3[b,Ind])/(s3[b,Ind+1]-s3[b,Ind-1])



          #############################################
          ####Sets up SigLam1 matrix for old density###
          #############################################


          W1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)
          Q1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)


          length1=diff(s3[b,])



          if(J3<2){
            if(J3==1){
              W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
              SigLam3=solve(diag(J3+1)-W1)%*%Q1

              do=log(dpois(J3,alpha3))+log(dmvnorm(lambda,rep(Mulam3[b],length(lambda)),SigLam3*Siglam3[b]))




            }else{


              do=log(dpois(J3,alpha3))+log(dnorm(lambda,Mulam3[b],Siglam3[b]))

            }
          }else{

            for(j in 2:J3){
              W1[j,j-1]=(clam3*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam3*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])


            SigLam3=solve(diag(J3+1)-W1)%*%Q1

            do=log(dpois(J3,alpha3))+log(dmvnorm(lambda,rep(Mulam3[b],length(lambda)),SigLam3*Siglam3[b]))


          }
          #############################################
          #############################################

          Lo=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s3[b,],lam1[b,],lam2[b,],lam3[b,],gam[b,])


          prior=((m3^2)*(s3[b,Ind+1]-s3[b,Ind-1]))/((2*J3+1)*(2*J3)*(s3[b,Ind]-s3[b,Ind-1])*(s3[b,Ind+1]-s3[b,Ind]))


          G3=G3-1
          J3=J3-1


          Ln=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s,lam1[b,],lam2[b,],lam,gam[b,])

          ###Make siglam matrix



          W1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)
          Q1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)



          length1=rep(0,J3+1)




          for(j in 1:length(length1)){
            length1[j]=s[j+1]-s[j]
          }


          if(J3<2){
            if(J3==1){
              W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
              SigLam3n=solve(diag(J3+1)-W1)%*%Q1


              dn=log(dpois(J3,alpha3))+log(dmvnorm(lam,rep(Mulam3[b],length(lam)),SigLam3n*Siglam3[b]))



            }else{

              dn=log(dpois(J3,alpha3))+log(dnorm(lam,Mulam3[b],Siglam3[b]))

            }
          }else{


            for(j in 2:J3){
              W1[j,j-1]=(clam3*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam3*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])


            SigLam3n=solve(diag(J3+1)-W1)%*%Q1

            dn=log(dpois(J3,alpha3))+log(dmvnorm(lam,rep(Mulam3[b],length(lam)),SigLam3n*Siglam3[b]))


          }
          ####




          alpha=Ln-Lo+dn-do+log(prior)+log(U3*(1-U3))

          if(is.nan(alpha)==TRUE){
            IndD3[b]=0
            J3=J3+1
            G3=G3+1
          }else{

            U=log(runif(1,0,1))

            iter[2]="AcceptRejDeath"

            if(U<alpha){
              s3[b,]=c(s,NA)
              IndD3[b]=1
              lam3[b,1:(J3+1)]=lam
              lam3[b,(J3+2):J3max]=rep(NA,J3max-J3-1)
            }else{
              IndD3[b]=0
              J3=J3+1
              G3=G3+1
            }
          }

          ####End else
        }
        split1[b]=J1
        split2[b]=J2
        split3[b]=J3
        ##
        sum1[b]=sum(eta1[b,])
        sum2[b]=sum(eta2[b,])
        sum3[b]=sum(eta3[b,])


      }





      ################End Samplers
      cat("




          ", "



          ", "


          ", "Posterior Inclusion Probabilities after half Burnin", "


          ", "Hazard 1", "

          ", colMeans(eta1[(B*burn+1):B,])*100, "

          ", "Hazard 2", "

          ", colMeans(eta2[(B*burn+1):B,])*100, "

          ", "Hazard 3", "

          ", colMeans(eta3[(B*burn+1):B,])*100,"

          ", "IndEta",mean(Indeta1[(B*burn+1):B])*100,mean(Indeta2[(B*burn+1):B])*100,mean(Indeta3[(B*burn+1):B])*100,"


          ","IndMix",mean(Indmix1[(B*burn+1):B])*100,mean(Indmix2[(B*burn+1):B])*100,mean(Indmix3[(B*burn+1):B])*100,"


          ", "Included Acceptance", "

          ", "Haz1", "
          ", mean(Ind1s[(B*burn+1):B])*100, "

          ", colMeans(Indcond1[(B*burn+1):B,],na.rm=TRUE)*100,"


          ",  "Haz2", "
          ",mean(Ind2s[(B*burn+1):B])*100, "

          ", colMeans(Indcond2[(B*burn+1):B,],na.rm=TRUE)*100,"



          ", "Haz3",mean(Ind3s[(B*burn+1):B])*100,"



          ", colMeans(Indcond1[(B*burn+1):B,],na.rm=TRUE)*100,"


          ","Survival","

          ","IndDeath",mean(IndD1[(B*burn+1):B])*100,mean(IndD2[(B*burn+1):B])*100,mean(IndD3[(B*burn+1):B])*100,"

          ","IndBirth",mean(IndB1[(B*burn+1):B])*100,mean(IndB2[(B*burn+1):B])*100,mean(IndB3[(B*burn+1):B])*100,"

          ","Lambda","

          ", "Lam1",
          colMeans(Acceptlam1[(B*burn+1):B,],na.rm=TRUE)*100,"

          ","Lam2",
          colMeans(Acceptlam2[(B*burn+1):B,],na.rm=TRUE)*100,"

          ","Lam3",
          colMeans(Acceptlam3[(B*burn+1):B,],na.rm=TRUE)*100,"

          ","Indepsilon",mean(Indepsilon[(B*burn+1):B])*100)




      Path1= paste0(Path,"/Ind1s.txt")

      write.table(Ind1s[(burn*B+1):B], Path1, sep="\t")






      Path1= paste0(Path,"/Ind2s.txt")

      write.table(Ind2s[(burn*B+1):B], Path1, sep="\t")





      Path1= paste0(Path,"/Ind3s.txt")

      write.table(Ind3s[(burn*B+1):B], Path1, sep="\t")

      par(mfrow=c(3,1))

      plot(1:B,sum1,type="l",xlab="",ylab="Haz1: # Included", main="Traceplot: # Included")
      plot(1:B,sum2,type="l",xlab="",ylab="Haz2: # Included")
      plot(1:B,sum3,type="l",xlab="",ylab="Haz3: # Included")


      plot(1:B,split1,type="l",xlab="",ylab="Haz1: Split #", main="Traceplot: # Split points")
      plot(1:B,split2,type="l",xlab="",ylab="Haz2: Split #")
      plot(1:B,split3,type="l",xlab="",ylab="Haz3: Split #")


    }



    ###If 0 inc



    if(inc==0){

      cat("No Variables Included")

      for(b in 2:B){


        if(b%%10000==0){cat(b, "iterations",date(), "  ")}else{
          if(b%%5000==0){cat(b, " iterations ")}}




        ###eta1,beta1
        eta1[b,]=eta1[b-1,]
        beta1[b,]=beta1[b-1,]

        if(sum(eta1[b-1,])==0|sum(eta1[b-1,])==p1){
          if(sum(eta1[b-1,])==0){
            ###Add Automatically
            Ind=sample(1:p1,1)
            eta1[b,Ind]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])
            ####

            meannew = 0
            varnew = sqrt(Sigmanew)
            ##################
            beta1[b,Ind]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta1[b,Ind],meannew,varnew))
            ######Accept reject###

            Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b-1,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo+dn + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
            U=log(runif(1,0,1))

            if(U>alphab1){
              eta1[b,]=eta1[b-1,]
              beta1[b,]=beta1[b-1,]
              Indeta1[b]=0
            }else{Indeta1[b]=1}
          }
          if(sum(eta1[b-1,])==p1){
            ###Delete Automatically
            Ind=sample(1:p1,1)
            eta1[b,Ind]=0
            beta1[b,Ind]=0
            includedold=rep(0,p1)
            for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(includedold[k]==Ind){spotold[k]=k}}
            spot1=max(spotold)


            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,includedold])%*%X[,includedold])

            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta1[b-1,includedold]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta1[b-1,Ind],meanold,varold))
            ######Accept reject###
            Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b-1,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo-do + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
            U=log(runif(1,0,1))

            if(U>alphab1){
              eta1[b,]=eta1[b-1,]
              beta1[b,]=beta1[b-1,]
              Indeta1[b]=0
            }else{Indeta1[b]=1}
          }
        }else{

          U=runif(1,0,1)

          if(U<psi){

            if(sum(eta1[b-1,])==1){

              includedold=rep(0,p1)
              for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
              includedold=includedold[includedold != 0]
              ones=includedold
              zeros=rep(0,p1)
              for(k in 1:p1){if(eta1[b-1,k]==0){zeros[k]=k}}
              zeros=zeros[zeros != 0]
              ###Sample swap indices###
              Indone=ones
              Indzero=sample(zeros,1)
              ####Change Beta/eta
              eta1[b,Indone]=0
              eta1[b,Indzero]=1
              beta1[b,Indone]=0
              ##

              Sigmaold=c*solve(t(X[,Indone])%*%X[,Indone])
              Sigmanew=c*solve(t(X[,Indzero])%*%X[,Indzero])

              meannew = 0
              varnew = sqrt(Sigmanew)
              ##################
              beta1[b,Indzero]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta1[b,Indzero],meannew,varnew))
              ###Old density###

              meanold = 0
              varold = sqrt(Sigmaold)
              do=log(dnorm(beta1[b-1,Indone],meanold,varold))
              ######Accept reject###
              Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b-1,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo+dn-do
              U=log(runif(1,0,1))

              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}
            }else{



              ###Swapper
              includedold=rep(0,p1)
              for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
              includedold=includedold[includedold != 0]
              ones=includedold
              zeros=rep(0,p1)
              for(k in 1:p1){if(eta1[b-1,k]==0){zeros[k]=k}}
              zeros=zeros[zeros != 0]
              ###Sample swap indices###
              Indone=sample(ones,1)
              Indzero=sample(zeros,1)
              ####Change Beta/eta
              eta1[b,Indone]=0
              eta1[b,Indzero]=1
              includednew=rep(0,p1)
              for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
              includednew=includednew[includednew != 0]
              spotold=rep(0,length(includedold))
              for(k in 1:length(includedold)){if(Indone==includedold[k]){spotold[k]=k}}
              spot1=max(spotold)
              spotnew=rep(0,length(includednew))
              for(k in 1:length(includednew)){if(Indzero==includednew[k]){spotnew[k]=k}}
              spot2=max(spotnew)
              ###Make sigma matrices##
              Sigmaold=c*solve(t(X[,includedold])%*%X[,includedold])
              Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])
              ###Generate new vector##
              beta1[b,Indone]=0

              ##meannew,varnew##
              V1 = Sigmanew[spot2,spot2]
              V2 = as.matrix(Sigmanew[-spot2,-spot2])
              V12 = as.matrix(Sigmanew[spot2,-spot2])
              thetab=beta1[b-1,includednew]
              thetano = as.matrix(thetab[-spot2])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              beta1[b,Indzero]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta1[b,Indzero],meannew,varnew))
              ###Old density###
              V1 = Sigmaold[spot1,spot1]
              V2 = as.matrix(Sigmaold[-spot1,-spot1])
              V12 = as.matrix(Sigmaold[spot1,-spot1])
              thetab=beta1[b-1,includedold]
              thetano = as.matrix(thetab[-spot1])
              meanold = t(V12)%*%solve(V2)%*%thetano
              varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              do=log(dnorm(beta1[b-1,Indone],meanold,varold))
              ######Accept reject###
              Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b-1,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo+dn-do
              U=log(runif(1,0,1))

              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}
            }

          }else{
            ###Add/Delete
            Ind=sample(1:p1,1)
            if(eta1[b-1,Ind]==1){
              ##delete##

              if(sum(eta1[b-1,])==1){


                eta1[b,Ind]=0
                beta1[b,Ind]=0
                includedold=rep(0,p1)
                for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
                includedold=includedold[includedold != 0]
                spotold=rep(0,length(includedold))
                for(k in 1:length(includedold)){if(Ind==includedold[k]){spotold[k]=k}}
                spot1=max(spotold)


                ###Make sigma matrices##
                Sigmaold=c*solve(t(X[,includedold])%*%X[,includedold])

                ###Old density###
                V1 = Sigmaold[spot1,spot1]
                V2 = as.matrix(Sigmaold[-spot1,-spot1])
                V12 = as.matrix(Sigmaold[spot1,-spot1])
                thetab=beta1[b-1,includedold]
                thetano = as.matrix(thetab[-spot1])
                meanold = 0
                varold = sqrt(Sigmaold)
                do=log(dnorm(beta1[b-1,Ind],meanold,varold))
                ######Accept reject###
                Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b-1,],beta2[b-1,],beta3[b-1,],
                          s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
                Liken=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                          s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

                alphab1=Liken-Likeo-do + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
                U=log(runif(1,0,1))

                if(U>alphab1){
                  eta1[b,]=eta1[b-1,]
                  beta1[b,]=beta1[b-1,]
                  Indeta1[b]=0
                }else{Indeta1[b]=1}



              }else{


                eta1[b,Ind]=0
                beta1[b,Ind]=0
                includedold=rep(0,p1)
                for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
                includedold=includedold[includedold != 0]
                spotold=rep(0,length(includedold))
                for(k in 1:length(includedold)){if(Ind==includedold[k]){spotold[k]=k}}
                spot1=max(spotold)


                ###Make sigma matrices##
                Sigmaold=c*solve(t(X[,includedold])%*%X[,includedold])

                ###Old density###
                V1 = Sigmaold[spot1,spot1]
                V2 = as.matrix(Sigmaold[-spot1,-spot1])
                V12 = as.matrix(Sigmaold[spot1,-spot1])
                thetab=beta1[b-1,includedold]
                thetano = as.matrix(thetab[-spot1])
                meanold = t(V12)%*%solve(V2)%*%thetano
                varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
                do=log(dnorm(beta1[b-1,Ind],meanold,varold))
                ######Accept reject###
                Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b-1,],beta2[b-1,],beta3[b-1,],
                          s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
                Liken=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                          s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

                alphab1=Liken-Likeo-do + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
                U=log(runif(1,0,1))

                if(U>alphab1){
                  eta1[b,]=eta1[b-1,]
                  beta1[b,]=beta1[b-1,]
                  Indeta1[b]=0
                }else{Indeta1[b]=1}

              }


            }else{
              ###Add###



              eta1[b,Ind]=1


              includednew=rep(0,p1)
              for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
              includednew=includednew[includednew != 0]
              spotnew=rep(0,length(includednew))
              for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
              spot2=max(spotnew)
              ###Make sigma matrices##
              Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])
              ####
              V1 = Sigmanew[spot2,spot2]
              V2 = as.matrix(Sigmanew[-spot2,-spot2])
              V12 = as.matrix(Sigmanew[spot2,-spot2])
              thetab=beta1[b-1,includednew]
              thetano = as.matrix(thetab[-spot2])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              beta1[b,Ind]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta1[b,Ind],meannew,varnew))
              ######Accept reject###
              Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b-1,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo+dn + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
              U=log(runif(1,0,1))

              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}



            }

          }

        }

        ##End Eta Beta

        includednew=rep(0,p1)
        for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
        includednew=includednew[includednew != 0]

        if(sum(eta1[b,])>0){


          if(sum(eta1[b,])==1){

            iter[2]="Conditional Inclusion"


            includednew=rep(0,p1)
            for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]

            Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])



            meannew = 0
            varnew = sqrt(Sigmanew)

            beta=beta1[b,]

            ##################
            beta[includednew]=rnorm(1,meannew,varnew)

            dn=log(dnorm(beta[includednew],meannew,varnew))
            ###density old
            do=log(dnorm(beta1[b,includednew],meannew,varnew))







            Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK1(Y1,Y2,I1,I2,X,beta,beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1m=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphab1m)==FALSE){
              Indcond1[b,includednew]=0
            }else{
              if(U>alphab1m){
                Indcond1[b,includednew]=0
              }else{Indcond1[b,includednew]=1
              beta1[b,]=beta
              }}

          }else{

            iter[2]="Conditional Inclusion"
            ##Jointly Update nonzero betas
            zeta1=beta1[b,]
            zeta1=zeta1[zeta1!=0]
            zeta1n=zeta1
            Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])



            ###############
            ####

            for(k in 1:length(includednew)){


              V1 = Sigmanew[k,k]
              V2 = as.matrix(Sigmanew[-k,-k])
              V12 = as.matrix(Sigmanew[k,-k])
              thetab=beta1[b,includednew]
              thetano = as.matrix(thetab[-k])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              zeta1n[k]=rnorm(1,meannew,varnew)
              dn=log(dnorm(zeta1n[k],meannew,varnew))
              ###density old
              do=log(dnorm(zeta1[k],meannew,varnew))






              beta=beta1[b,]
              beta[includednew]=zeta1n


              Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK1(Y1,Y2,I1,I2,X,beta,beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1m=Liken-Likeo+dn -do
              U=log(runif(1,0,1))
              if(is.finite(alphab1m)==FALSE){
                Indcond1[b,includednew[k]]=0
              }else{
                if(U>alphab1m){
                  Indcond1[b,includednew[k]]=0
                  zeta1n[k]=zeta1[k]
                }else{Indcond1[b,includednew[k]]=1
                beta1[b,]=beta
                zeta1[k]=zeta1n[k]
                }}

            }



            ##Jointly Update nonzero betas
            iter[2]="mixing"
            zeta1n=beta1[b,]
            Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])


            zeta1n[includednew]=rmvnorm(1,rep(0,length(includednew)),Sigmanew)

            beta=beta1[b,]
            beta=beta[beta!=0]

            dn=log(dmvnorm(zeta1n[includednew],rep(0,length(includednew)),Sigmanew))
            ###density old
            do=log(dmvnorm(beta,rep(0,length(includednew)),Sigmanew))

            ######Accept reject###
            Likeo=LK1(Y1,Y2,I1,I2,X,beta1[b,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK1(Y1,Y2,I1,I2,X,zeta1n,beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphamix1=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphamix1)==FALSE){
              Indmix1[b]=0
            }else{
              if(U>alphamix1){

                Indmix1[b]=0
              }else{Indmix1[b]=1
              beta1[b,]=zeta1n
              }}

          }



        }





        ###eta2,beta2
        eta2[b,]=eta2[b-1,]
        beta2[b,]=beta2[b-1,]

        if(sum(eta2[b-1,])==0|sum(eta2[b-1,])==p1){
          if(sum(eta2[b-1,])==0){
            ###Add Automatically
            Ind=sample(1:p1,1)
            eta2[b,Ind]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta2[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])
            ####

            meannew = 0
            varnew = sqrt(Sigmanew)
            ##################
            beta2[b,Ind]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta2[b,Ind],meannew,varnew))
            ######Accept reject###

            Likeo=LK2(Y1,Y2,I1,I2,X,beta2[b-1,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK2(Y1,Y2,I1,I2,X,beta2[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo+dn + log(beta(sum(eta2[b,])+z2a,p1-sum(eta2[b,])+z2b)) - log(beta(sum(eta2[b-1,])+z2a,p1-sum(eta2[b-1,])+z2b))
            U=log(runif(1,0,1))

            if(U>alphab1){
              eta2[b,]=eta2[b-1,]
              beta2[b,]=beta2[b-1,]
              Indeta2[b]=0
            }else{Indeta2[b]=1}
          }
          if(sum(eta2[b-1,])==p1){
            ###Delete Automatically
            Ind=sample(1:p1,1)
            eta2[b,Ind]=0
            beta2[b,Ind]=0
            includedold=rep(0,p1)
            for(k in 1:p1){if(eta2[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(includedold[k]==Ind){spotold[k]=k}}
            spot1=max(spotold)


            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,includedold])%*%X[,includedold])

            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta2[b-1,includedold]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta2[b-1,Ind],meanold,varold))
            ######Accept reject###
            Likeo=LK2(Y1,Y2,I1,I2,X,beta2[b-1,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK2(Y1,Y2,I1,I2,X,beta2[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo-do + log(beta(sum(eta2[b,])+z2a,p1-sum(eta2[b,])+z2b)) - log(beta(sum(eta2[b-1,])+z2a,p1-sum(eta2[b-1,])+z2b))
            U=log(runif(1,0,1))

            if(U>alphab1){
              eta2[b,]=eta2[b-1,]
              beta2[b,]=beta2[b-1,]
              Indeta2[b]=0
            }else{Indeta2[b]=1}
          }
        }else{

          U=runif(1,0,1)

          if(U<psi){

            if(sum(eta2[b-1,])==1){

              includedold=rep(0,p1)
              for(k in 1:p1){if(eta2[b-1,k]==1){includedold[k]=k}}
              includedold=includedold[includedold != 0]
              ones=includedold
              zeros=rep(0,p1)
              for(k in 1:p1){if(eta2[b-1,k]==0){zeros[k]=k}}
              zeros=zeros[zeros != 0]
              ###Sample swap indices###
              Indone=ones
              Indzero=sample(zeros,1)
              ####Change Beta/eta
              eta2[b,Indone]=0
              eta2[b,Indzero]=1
              beta2[b,Indone]=0
              ##

              Sigmaold=c*solve(t(X[,Indone])%*%X[,Indone])
              Sigmanew=c*solve(t(X[,Indzero])%*%X[,Indzero])

              meannew = 0
              varnew = sqrt(Sigmanew)
              ##################
              beta2[b,Indzero]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta2[b,Indzero],meannew,varnew))
              ###Old density###

              meanold = 0
              varold = sqrt(Sigmaold)
              do=log(dnorm(beta2[b-1,Indone],meanold,varold))
              ######Accept reject###
              Likeo=LK2(Y1,Y2,I1,I2,X,beta2[b-1,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK2(Y1,Y2,I1,I2,X,beta2[b,],beta2[b,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo+dn-do
              U=log(runif(1,0,1))

              if(U>alphab1){
                eta2[b,]=eta2[b-1,]
                beta2[b,]=beta2[b-1,]
                Indeta2[b]=0
              }else{Indeta2[b]=1}
            }else{



              ###Swapper
              includedold=rep(0,p1)
              for(k in 1:p1){if(eta2[b-1,k]==1){includedold[k]=k}}
              includedold=includedold[includedold != 0]
              ones=includedold
              zeros=rep(0,p1)
              for(k in 1:p1){if(eta2[b-1,k]==0){zeros[k]=k}}
              zeros=zeros[zeros != 0]
              ###Sample swap indices###
              Indone=sample(ones,1)
              Indzero=sample(zeros,1)
              ####Change Beta/eta
              eta2[b,Indone]=0
              eta2[b,Indzero]=1
              includednew=rep(0,p1)
              for(k in 1:p1){if(eta2[b,k]==1){includednew[k]=k}}
              includednew=includednew[includednew != 0]
              spotold=rep(0,length(includedold))
              for(k in 1:length(includedold)){if(Indone==includedold[k]){spotold[k]=k}}
              spot1=max(spotold)
              spotnew=rep(0,length(includednew))
              for(k in 1:length(includednew)){if(Indzero==includednew[k]){spotnew[k]=k}}
              spot2=max(spotnew)
              ###Make sigma matrices##
              Sigmaold=c*solve(t(X[,includedold])%*%X[,includedold])
              Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])
              ###Generate new vector##
              beta2[b,Indone]=0

              ##meannew,varnew##
              V1 = Sigmanew[spot2,spot2]
              V2 = as.matrix(Sigmanew[-spot2,-spot2])
              V12 = as.matrix(Sigmanew[spot2,-spot2])
              thetab=beta2[b-1,includednew]
              thetano = as.matrix(thetab[-spot2])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              beta2[b,Indzero]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta2[b,Indzero],meannew,varnew))
              ###Old density###
              V1 = Sigmaold[spot1,spot1]
              V2 = as.matrix(Sigmaold[-spot1,-spot1])
              V12 = as.matrix(Sigmaold[spot1,-spot1])
              thetab=beta2[b-1,includedold]
              thetano = as.matrix(thetab[-spot1])
              meanold = t(V12)%*%solve(V2)%*%thetano
              varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              do=log(dnorm(beta2[b-1,Indone],meanold,varold))
              ######Accept reject###
              Likeo=LK2(Y1,Y2,I1,I2,X,beta2[b-1,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK2(Y1,Y2,I1,I2,X,beta2[b,],beta2[b,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo+dn-do
              U=log(runif(1,0,1))

              if(U>alphab1){
                eta2[b,]=eta2[b-1,]
                beta2[b,]=beta2[b-1,]
                Indeta2[b]=0
              }else{Indeta2[b]=1}
            }

          }else{
            ###Add/Delete
            Ind=sample(1:p1,1)
            if(eta2[b-1,Ind]==1){
              ##delete##

              if(sum(eta2[b-1,])==1){


                eta2[b,Ind]=0
                beta2[b,Ind]=0
                includedold=rep(0,p1)
                for(k in 1:p1){if(eta2[b-1,k]==1){includedold[k]=k}}
                includedold=includedold[includedold != 0]
                spotold=rep(0,length(includedold))
                for(k in 1:length(includedold)){if(Ind==includedold[k]){spotold[k]=k}}
                spot1=max(spotold)


                ###Make sigma matrices##
                Sigmaold=c*solve(t(X[,includedold])%*%X[,includedold])

                ###Old density###
                V1 = Sigmaold[spot1,spot1]
                V2 = as.matrix(Sigmaold[-spot1,-spot1])
                V12 = as.matrix(Sigmaold[spot1,-spot1])
                thetab=beta2[b-1,includedold]
                thetano = as.matrix(thetab[-spot1])
                meanold = 0
                varold = sqrt(Sigmaold)
                do=log(dnorm(beta2[b-1,Ind],meanold,varold))
                ######Accept reject###
                Likeo=LK2(Y1,Y2,I1,I2,X,beta2[b-1,],beta2[b-1,],beta3[b-1,],
                          s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
                Liken=LK2(Y1,Y2,I1,I2,X,beta2[b,],beta2[b,],beta3[b-1,],
                          s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

                alphab1=Liken-Likeo-do + log(beta(sum(eta2[b,])+z2a,p1-sum(eta2[b,])+z2b)) - log(beta(sum(eta2[b-1,])+z2a,p1-sum(eta2[b-1,])+z2b))
                U=log(runif(1,0,1))

                if(U>alphab1){
                  eta2[b,]=eta2[b-1,]
                  beta2[b,]=beta2[b-1,]
                  Indeta2[b]=0
                }else{Indeta2[b]=1}



              }else{


                eta2[b,Ind]=0
                beta2[b,Ind]=0
                includedold=rep(0,p1)
                for(k in 1:p1){if(eta2[b-1,k]==1){includedold[k]=k}}
                includedold=includedold[includedold != 0]
                spotold=rep(0,length(includedold))
                for(k in 1:length(includedold)){if(Ind==includedold[k]){spotold[k]=k}}
                spot1=max(spotold)


                ###Make sigma matrices##
                Sigmaold=c*solve(t(X[,includedold])%*%X[,includedold])

                ###Old density###
                V1 = Sigmaold[spot1,spot1]
                V2 = as.matrix(Sigmaold[-spot1,-spot1])
                V12 = as.matrix(Sigmaold[spot1,-spot1])
                thetab=beta2[b-1,includedold]
                thetano = as.matrix(thetab[-spot1])
                meanold = t(V12)%*%solve(V2)%*%thetano
                varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
                do=log(dnorm(beta2[b-1,Ind],meanold,varold))
                ######Accept reject###
                Likeo=LK2(Y1,Y2,I1,I2,X,beta2[b-1,],beta2[b-1,],beta3[b-1,],
                          s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
                Liken=LK2(Y1,Y2,I1,I2,X,beta2[b,],beta2[b,],beta3[b-1,],
                          s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

                alphab1=Liken-Likeo-do + log(beta(sum(eta2[b,])+z2a,p1-sum(eta2[b,])+z2b)) - log(beta(sum(eta2[b-1,])+z2a,p1-sum(eta2[b-1,])+z2b))
                U=log(runif(1,0,1))

                if(U>alphab1){
                  eta2[b,]=eta2[b-1,]
                  beta2[b,]=beta2[b-1,]
                  Indeta2[b]=0
                }else{Indeta2[b]=1}

              }


            }else{
              ###Add###



              eta2[b,Ind]=1


              includednew=rep(0,p1)
              for(k in 1:p1){if(eta2[b,k]==1){includednew[k]=k}}
              includednew=includednew[includednew != 0]
              spotnew=rep(0,length(includednew))
              for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
              spot2=max(spotnew)
              ###Make sigma matrices##
              Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])
              ####
              V1 = Sigmanew[spot2,spot2]
              V2 = as.matrix(Sigmanew[-spot2,-spot2])
              V12 = as.matrix(Sigmanew[spot2,-spot2])
              thetab=beta2[b-1,includednew]
              thetano = as.matrix(thetab[-spot2])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              beta2[b,Ind]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta2[b,Ind],meannew,varnew))
              ######Accept reject###
              Likeo=LK2(Y1,Y2,I1,I2,X,beta2[b-1,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK2(Y1,Y2,I1,I2,X,beta2[b,],beta2[b,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo+dn + log(beta(sum(eta2[b,])+z2a,p1-sum(eta2[b,])+z2b)) - log(beta(sum(eta2[b-1,])+z2a,p1-sum(eta2[b-1,])+z2b))
              U=log(runif(1,0,1))

              if(U>alphab1){
                eta2[b,]=eta2[b-1,]
                beta2[b,]=beta2[b-1,]
                Indeta2[b]=0
              }else{Indeta2[b]=1}



            }

          }

        }

        ##End Eta Beta

        includednew=rep(0,p1)
        for(k in 1:p1){if(eta2[b,k]==1){includednew[k]=k}}
        includednew=includednew[includednew != 0]

        if(sum(eta2[b,])>0){


          if(sum(eta2[b,])==1){

            iter[2]="Conditional Inclusion"


            includednew=rep(0,p1)
            for(k in 1:p1){if(eta2[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]

            Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])



            meannew = 0
            varnew = sqrt(Sigmanew)

            beta=beta2[b,]

            ##################
            beta[includednew]=rnorm(1,meannew,varnew)

            dn=log(dnorm(beta[includednew],meannew,varnew))
            ###density old
            do=log(dnorm(beta2[b,includednew],meannew,varnew))







            Likeo=LK2(Y1,Y2,I1,I2,X,beta2[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK2(Y1,Y2,I1,I2,X,beta,beta,beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1m=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphab1m)==FALSE){
              Indcond2[b,includednew]=0
            }else{
              if(U>alphab1m){
                Indcond2[b,includednew]=0
              }else{Indcond2[b,includednew]=1
              beta2[b,]=beta
              }}

          }else{

            iter[2]="Conditional Inclusion"
            ##Jointly Update nonzero betas
            zeta2=beta2[b,]
            zeta2=zeta2[zeta2!=0]
            zeta2n=zeta2
            Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])



            ###############
            ####

            for(k in 1:length(includednew)){


              V1 = Sigmanew[k,k]
              V2 = as.matrix(Sigmanew[-k,-k])
              V12 = as.matrix(Sigmanew[k,-k])
              thetab=beta2[b,includednew]
              thetano = as.matrix(thetab[-k])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              zeta2n[k]=rnorm(1,meannew,varnew)
              dn=log(dnorm(zeta2n[k],meannew,varnew))
              ###density old
              do=log(dnorm(zeta2[k],meannew,varnew))






              beta=beta2[b,]
              beta[includednew]=zeta2n


              Likeo=LK2(Y1,Y2,I1,I2,X,beta2[b,],beta2[b,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK2(Y1,Y2,I1,I2,X,beta,beta,beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1m=Liken-Likeo+dn -do
              U=log(runif(1,0,1))
              if(is.finite(alphab1m)==FALSE){
                Indcond2[b,includednew[k]]=0
              }else{
                if(U>alphab1m){
                  Indcond2[b,includednew[k]]=0
                  zeta2n[k]=zeta2[k]
                }else{Indcond2[b,includednew[k]]=1
                beta2[b,]=beta
                zeta2[k]=zeta2n[k]
                }}

            }



            ##Jointly Update nonzero betas
            iter[2]="mixing"
            zeta2n=beta2[b,]
            Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])


            zeta2n[includednew]=rmvnorm(1,rep(0,length(includednew)),Sigmanew)

            beta=beta2[b,]
            beta=beta[beta!=0]

            dn=log(dmvnorm(zeta2n[includednew],rep(0,length(includednew)),Sigmanew))
            ###density old
            do=log(dmvnorm(beta,rep(0,length(includednew)),Sigmanew))

            ######Accept reject###
            Likeo=LK2(Y1,Y2,I1,I2,X,beta2[b,],beta2[b,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK2(Y1,Y2,I1,I2,X,zeta2n,zeta2n,beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphamix1=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphamix1)==FALSE){
              Indmix2[b]=0
            }else{
              if(U>alphamix1){

                Indmix2[b]=0
              }else{Indmix2[b]=1
              beta2[b,]=zeta2n
              }}

          }



        }



        ###eta3,beta3
        eta3[b,]=eta3[b-1,]
        beta3[b,]=beta3[b-1,]

        if(sum(eta3[b-1,])==0|sum(eta3[b-1,])==p1){
          if(sum(eta3[b-1,])==0){
            ###Add Automatically
            Ind=sample(1:p1,1)
            eta3[b,Ind]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta3[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])
            ####

            meannew = 0
            varnew = sqrt(Sigmanew)
            ##################
            beta3[b,Ind]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta3[b,Ind],meannew,varnew))
            ######Accept reject###

            Likeo=LK3(Y1,Y2,I1,I2,X,beta3[b-1,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK3(Y1,Y2,I1,I2,X,beta3[b,],beta2[b-1,],beta3[b,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo+dn + log(beta(sum(eta3[b,])+z3a,p1-sum(eta3[b,])+z3b)) - log(beta(sum(eta3[b-1,])+z3a,p1-sum(eta3[b-1,])+z3b))
            U=log(runif(1,0,1))

            if(U>alphab1){
              eta3[b,]=eta3[b-1,]
              beta3[b,]=beta3[b-1,]
              Indeta3[b]=0
            }else{Indeta3[b]=1}
          }
          if(sum(eta3[b-1,])==p1){
            ###Delete Automatically
            Ind=sample(1:p1,1)
            eta3[b,Ind]=0
            beta3[b,Ind]=0
            includedold=rep(0,p1)
            for(k in 1:p1){if(eta3[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(includedold[k]==Ind){spotold[k]=k}}
            spot1=max(spotold)


            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,includedold])%*%X[,includedold])

            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta3[b-1,includedold]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta3[b-1,Ind],meanold,varold))
            ######Accept reject###
            Likeo=LK3(Y1,Y2,I1,I2,X,beta3[b-1,],beta2[b-1,],beta3[b-1,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK3(Y1,Y2,I1,I2,X,beta3[b,],beta2[b-1,],beta3[b,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1=Liken-Likeo-do + log(beta(sum(eta3[b,])+z3a,p1-sum(eta3[b,])+z3b)) - log(beta(sum(eta3[b-1,])+z3a,p1-sum(eta3[b-1,])+z3b))
            U=log(runif(1,0,1))

            if(U>alphab1){
              eta3[b,]=eta3[b-1,]
              beta3[b,]=beta3[b-1,]
              Indeta3[b]=0
            }else{Indeta3[b]=1}
          }
        }else{

          U=runif(1,0,1)

          if(U<psi){

            if(sum(eta3[b-1,])==1){

              includedold=rep(0,p1)
              for(k in 1:p1){if(eta3[b-1,k]==1){includedold[k]=k}}
              includedold=includedold[includedold != 0]
              ones=includedold
              zeros=rep(0,p1)
              for(k in 1:p1){if(eta3[b-1,k]==0){zeros[k]=k}}
              zeros=zeros[zeros != 0]
              ###Sample swap indices###
              Indone=ones
              Indzero=sample(zeros,1)
              ####Change Beta/eta
              eta3[b,Indone]=0
              eta3[b,Indzero]=1
              beta3[b,Indone]=0
              ##

              Sigmaold=c*solve(t(X[,Indone])%*%X[,Indone])
              Sigmanew=c*solve(t(X[,Indzero])%*%X[,Indzero])

              meannew = 0
              varnew = sqrt(Sigmanew)
              ##################
              beta3[b,Indzero]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta3[b,Indzero],meannew,varnew))
              ###Old density###

              meanold = 0
              varold = sqrt(Sigmaold)
              do=log(dnorm(beta3[b-1,Indone],meanold,varold))
              ######Accept reject###
              Likeo=LK3(Y1,Y2,I1,I2,X,beta3[b-1,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK3(Y1,Y2,I1,I2,X,beta3[b,],beta2[b-1,],beta3[b,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo+dn-do
              U=log(runif(1,0,1))

              if(U>alphab1){
                eta3[b,]=eta3[b-1,]
                beta3[b,]=beta3[b-1,]
                Indeta3[b]=0
              }else{Indeta3[b]=1}
            }else{



              ###Swapper
              includedold=rep(0,p1)
              for(k in 1:p1){if(eta3[b-1,k]==1){includedold[k]=k}}
              includedold=includedold[includedold != 0]
              ones=includedold
              zeros=rep(0,p1)
              for(k in 1:p1){if(eta3[b-1,k]==0){zeros[k]=k}}
              zeros=zeros[zeros != 0]
              ###Sample swap indices###
              Indone=sample(ones,1)
              Indzero=sample(zeros,1)
              ####Change Beta/eta
              eta3[b,Indone]=0
              eta3[b,Indzero]=1
              includednew=rep(0,p1)
              for(k in 1:p1){if(eta3[b,k]==1){includednew[k]=k}}
              includednew=includednew[includednew != 0]
              spotold=rep(0,length(includedold))
              for(k in 1:length(includedold)){if(Indone==includedold[k]){spotold[k]=k}}
              spot1=max(spotold)
              spotnew=rep(0,length(includednew))
              for(k in 1:length(includednew)){if(Indzero==includednew[k]){spotnew[k]=k}}
              spot2=max(spotnew)
              ###Make sigma matrices##
              Sigmaold=c*solve(t(X[,includedold])%*%X[,includedold])
              Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])
              ###Generate new vector##
              beta3[b,Indone]=0

              ##meannew,varnew##
              V1 = Sigmanew[spot2,spot2]
              V2 = as.matrix(Sigmanew[-spot2,-spot2])
              V12 = as.matrix(Sigmanew[spot2,-spot2])
              thetab=beta3[b-1,includednew]
              thetano = as.matrix(thetab[-spot2])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              beta3[b,Indzero]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta3[b,Indzero],meannew,varnew))
              ###Old density###
              V1 = Sigmaold[spot1,spot1]
              V2 = as.matrix(Sigmaold[-spot1,-spot1])
              V12 = as.matrix(Sigmaold[spot1,-spot1])
              thetab=beta3[b-1,includedold]
              thetano = as.matrix(thetab[-spot1])
              meanold = t(V12)%*%solve(V2)%*%thetano
              varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              do=log(dnorm(beta3[b-1,Indone],meanold,varold))
              ######Accept reject###
              Likeo=LK3(Y1,Y2,I1,I2,X,beta3[b-1,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK3(Y1,Y2,I1,I2,X,beta3[b,],beta2[b-1,],beta3[b,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo+dn-do
              U=log(runif(1,0,1))

              if(U>alphab1){
                eta3[b,]=eta3[b-1,]
                beta3[b,]=beta3[b-1,]
                Indeta3[b]=0
              }else{Indeta3[b]=1}
            }

          }else{
            ###Add/Delete
            Ind=sample(1:p1,1)
            if(eta3[b-1,Ind]==1){
              ##delete##

              if(sum(eta3[b-1,])==1){


                eta3[b,Ind]=0
                beta3[b,Ind]=0
                includedold=rep(0,p1)
                for(k in 1:p1){if(eta3[b-1,k]==1){includedold[k]=k}}
                includedold=includedold[includedold != 0]
                spotold=rep(0,length(includedold))
                for(k in 1:length(includedold)){if(Ind==includedold[k]){spotold[k]=k}}
                spot1=max(spotold)


                ###Make sigma matrices##
                Sigmaold=c*solve(t(X[,includedold])%*%X[,includedold])

                ###Old density###
                V1 = Sigmaold[spot1,spot1]
                V2 = as.matrix(Sigmaold[-spot1,-spot1])
                V12 = as.matrix(Sigmaold[spot1,-spot1])
                thetab=beta3[b-1,includedold]
                thetano = as.matrix(thetab[-spot1])
                meanold = 0
                varold = sqrt(Sigmaold)
                do=log(dnorm(beta3[b-1,Ind],meanold,varold))
                ######Accept reject###
                Likeo=LK3(Y1,Y2,I1,I2,X,beta3[b-1,],beta2[b-1,],beta3[b-1,],
                          s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
                Liken=LK3(Y1,Y2,I1,I2,X,beta3[b,],beta2[b-1,],beta3[b,],
                          s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

                alphab1=Liken-Likeo-do + log(beta(sum(eta3[b,])+z3a,p1-sum(eta3[b,])+z3b)) - log(beta(sum(eta3[b-1,])+z3a,p1-sum(eta3[b-1,])+z3b))
                U=log(runif(1,0,1))

                if(U>alphab1){
                  eta3[b,]=eta3[b-1,]
                  beta3[b,]=beta3[b-1,]
                  Indeta3[b]=0
                }else{Indeta3[b]=1}



              }else{


                eta3[b,Ind]=0
                beta3[b,Ind]=0
                includedold=rep(0,p1)
                for(k in 1:p1){if(eta3[b-1,k]==1){includedold[k]=k}}
                includedold=includedold[includedold != 0]
                spotold=rep(0,length(includedold))
                for(k in 1:length(includedold)){if(Ind==includedold[k]){spotold[k]=k}}
                spot1=max(spotold)


                ###Make sigma matrices##
                Sigmaold=c*solve(t(X[,includedold])%*%X[,includedold])

                ###Old density###
                V1 = Sigmaold[spot1,spot1]
                V2 = as.matrix(Sigmaold[-spot1,-spot1])
                V12 = as.matrix(Sigmaold[spot1,-spot1])
                thetab=beta3[b-1,includedold]
                thetano = as.matrix(thetab[-spot1])
                meanold = t(V12)%*%solve(V2)%*%thetano
                varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
                do=log(dnorm(beta3[b-1,Ind],meanold,varold))
                ######Accept reject###
                Likeo=LK3(Y1,Y2,I1,I2,X,beta3[b-1,],beta2[b-1,],beta3[b-1,],
                          s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
                Liken=LK3(Y1,Y2,I1,I2,X,beta3[b,],beta2[b-1,],beta3[b,],
                          s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

                alphab1=Liken-Likeo-do + log(beta(sum(eta3[b,])+z3a,p1-sum(eta3[b,])+z3b)) - log(beta(sum(eta3[b-1,])+z3a,p1-sum(eta3[b-1,])+z3b))
                U=log(runif(1,0,1))

                if(U>alphab1){
                  eta3[b,]=eta3[b-1,]
                  beta3[b,]=beta3[b-1,]
                  Indeta3[b]=0
                }else{Indeta3[b]=1}

              }


            }else{
              ###Add###



              eta3[b,Ind]=1


              includednew=rep(0,p1)
              for(k in 1:p1){if(eta3[b,k]==1){includednew[k]=k}}
              includednew=includednew[includednew != 0]
              spotnew=rep(0,length(includednew))
              for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
              spot2=max(spotnew)
              ###Make sigma matrices##
              Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])
              ####
              V1 = Sigmanew[spot2,spot2]
              V2 = as.matrix(Sigmanew[-spot2,-spot2])
              V12 = as.matrix(Sigmanew[spot2,-spot2])
              thetab=beta3[b-1,includednew]
              thetano = as.matrix(thetab[-spot2])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              beta3[b,Ind]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta3[b,Ind],meannew,varnew))
              ######Accept reject###
              Likeo=LK3(Y1,Y2,I1,I2,X,beta3[b-1,],beta2[b-1,],beta3[b-1,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK3(Y1,Y2,I1,I2,X,beta3[b,],beta2[b-1,],beta3[b,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1=Liken-Likeo+dn + log(beta(sum(eta3[b,])+z3a,p1-sum(eta3[b,])+z3b)) - log(beta(sum(eta3[b-1,])+z3a,p1-sum(eta3[b-1,])+z3b))
              U=log(runif(1,0,1))

              if(U>alphab1){
                eta3[b,]=eta3[b-1,]
                beta3[b,]=beta3[b-1,]
                Indeta3[b]=0
              }else{Indeta3[b]=1}



            }

          }

        }

        ##End Eta Beta

        includednew=rep(0,p1)
        for(k in 1:p1){if(eta3[b,k]==1){includednew[k]=k}}
        includednew=includednew[includednew != 0]

        if(sum(eta3[b,])>0){


          if(sum(eta3[b,])==1){

            iter[2]="Conditional Inclusion"


            includednew=rep(0,p1)
            for(k in 1:p1){if(eta3[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]

            Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])



            meannew = 0
            varnew = sqrt(Sigmanew)

            beta=beta3[b,]

            ##################
            beta[includednew]=rnorm(1,meannew,varnew)

            dn=log(dnorm(beta[includednew],meannew,varnew))
            ###density old
            do=log(dnorm(beta3[b,includednew],meannew,varnew))







            Likeo=LK3(Y1,Y2,I1,I2,X,beta3[b,],beta2[b-1,],beta3[b,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK3(Y1,Y2,I1,I2,X,beta,beta2[b-1,],beta,
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphab1m=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphab1m)==FALSE){
              Indcond3[b,includednew]=0
            }else{
              if(U>alphab1m){
                Indcond3[b,includednew]=0
              }else{Indcond3[b,includednew]=1
              beta3[b,]=beta
              }}

          }else{

            iter[2]="Conditional Inclusion"
            ##Jointly Update nonzero betas
            zeta3=beta3[b,]
            zeta3=zeta3[zeta3!=0]
            zeta3n=zeta3
            Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])



            ###############
            ####

            for(k in 1:length(includednew)){


              V1 = Sigmanew[k,k]
              V2 = as.matrix(Sigmanew[-k,-k])
              V12 = as.matrix(Sigmanew[k,-k])
              thetab=beta3[b,includednew]
              thetano = as.matrix(thetab[-k])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              zeta3n[k]=rnorm(1,meannew,varnew)
              dn=log(dnorm(zeta3n[k],meannew,varnew))
              ###density old
              do=log(dnorm(zeta3[k],meannew,varnew))






              beta=beta3[b,]
              beta[includednew]=zeta3n


              Likeo=LK3(Y1,Y2,I1,I2,X,beta3[b,],beta2[b-1,],beta3[b,],
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
              Liken=LK3(Y1,Y2,I1,I2,X,beta,beta2[b-1,],beta,
                        s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

              alphab1m=Liken-Likeo+dn -do
              U=log(runif(1,0,1))
              if(is.finite(alphab1m)==FALSE){
                Indcond3[b,includednew[k]]=0
              }else{
                if(U>alphab1m){
                  Indcond3[b,includednew[k]]=0
                  zeta3n[k]=zeta3[k]
                }else{Indcond3[b,includednew[k]]=1
                beta3[b,]=beta
                zeta3[k]=zeta3n[k]
                }}

            }



            ##Jointly Update nonzero betas
            iter[2]="mixing"
            zeta3n=beta3[b,]
            Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])


            zeta3n[includednew]=rmvnorm(1,rep(0,length(includednew)),Sigmanew)

            beta=beta3[b,]
            beta=beta[beta!=0]

            dn=log(dmvnorm(zeta3n[includednew],rep(0,length(includednew)),Sigmanew))
            ###density old
            do=log(dmvnorm(beta,rep(0,length(includednew)),Sigmanew))

            ######Accept reject###
            Likeo=LK3(Y1,Y2,I1,I2,X,beta3[b,],beta2[b-1,],beta3[b,],
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])
            Liken=LK3(Y1,Y2,I1,I2,X,zeta3n,beta2[b-1,],zeta3n,
                      s1[b-1,],s2[b-1,],s3[b-1,],lam1[b-1,],lam2[b-1,],lam3[b-1,],gam[b-1,])

            alphamix1=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphamix1)==FALSE){
              Indmix3[b]=0
            }else{
              if(U>alphamix1){

                Indmix3[b]=0
              }else{Indmix3[b]=1
              beta3[b,]=zeta3n
              }}

          }



        }

















        ########################
        ###Frailty Samplers#####
        ########################
        ############Epsilon Sampler#####

        iter[1]="frailty"
        iter[2]="hier"
        Der1o=D1(epsilon[b-1],gam[b-1,])
        Der2o=D2(epsilon[b-1])

        epsilon[b]=rgamma(1,((epsilon[b-1]-min(0,Der1o/Der2o))^2)/(-(cep^2)/Der2o),rate=(epsilon[b-1]-min(0,Der1o/Der2o))/(-(cep^2)/Der2o))

        Der1n=D1(epsilon[b],gam[b-1,])
        Der2n=D2(epsilon[b])

        dn=dgamma(epsilon[b-1],((epsilon[b]-min(0,Der1n/Der2n))^2)/(-(cep^2)/Der2n),rate=(epsilon[b]-min(0,Der1n/Der2n))/(-(cep^2)/Der2n))
        do=dgamma(epsilon[b],((epsilon[b-1]-min(0,Der1o/Der2o))^2)/(-(cep^2)/Der2o),rate=(epsilon[b-1]-min(0,Der1o/Der2o))/(-(cep^2)/Der2o))
        pn=(n*epsilon[b]+psi1-1)*log(epsilon[b])-epsilon[b]*(sum(gam[b-1,])+w)+(epsilon[b]-1)*sum(log(gam[b-1,]))-n*log(gamma(epsilon[b]))
        po=(n*epsilon[b-1]+psi1-1)*log(epsilon[b-1])-epsilon[b-1]*(sum(gam[b-1,])+w)+(epsilon[b-1]-1)*sum(log(gam[b-1,]))-n*log(gamma(epsilon[b-1]))




        alphaep=log(dn)-log(do)+pn-po

        if(is.nan(alphaep)==TRUE){
          epsilon[b]=epsilon[b-1]
          Indepsilon[b]=0
        }else{
          U=log(runif(1,0,1))

          if(U>alphaep){
            epsilon[b]=epsilon[b-1]
            Indepsilon[b]=0
          }else{Indepsilon[b]=1}
        }

        ####Frailty Sampler here
        ####Gam here is not how it's done
        iter[2]="gamma"

        S1=s1[b-1,]
        S1=S1[!is.na(S1)]
        S2=s2[b-1,]
        S2=S2[!is.na(S2)]
        S3=s3[b-1,]
        S3=S3[!is.na(S3)]

        L1=lam1[b-1,]
        L1=as.matrix(L1[!is.na(L1)])
        L2=lam2[b-1,]
        L2=as.matrix(L2[!is.na(L2)])
        L3=lam3[b-1,]
        L3=as.matrix(L3[!is.na(L3)])

        phi1=phifun(Y1,Y1,I1,I2,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),S1,S2,S3,
                    L1,L2,L3,epsilon[b],X)
        ##Sample
        for(i in 1:n){
          gam[b,i]=rgamma(1,1/epsilon[b]+I1[i]+I2[i],rate=phi1[i])
        }


        ############################################
        #####Start LogBH Samplers###################
        ############################################
        ####Lam1####

        iter[1]="LogBH1"
        iter[2]="matrixsetup"

        W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
        Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


        length1=rep(0,J1+1)




        for(j in 1:length(length1)){
          length1[j]=s1[b-1,j+1]-s1[b-1,j]
        }


        if(J1<2){
          if(J1==1){
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            SigLam1=solve(diag(J1+1)-W1)%*%Q1





          }else{

            Q1=as.matrix(2/(m1))
            SigLam1=Q1
          }
        }else{


          for(j in 2:J1){
            W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
          }


          Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
          Q1[1,1]=2/(2*length1[1]+length1[2])
          W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
          W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


          SigLam1=solve(diag(J1+1)-W1)%*%Q1

        }



        iter[2]="Mu"
        ##Lambda1 Hierarchical Sampler
        ##Mulam

        if(J1>0){

          Mulam1[b]=rnorm(1,(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%L1)/(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%as.matrix(rep(1,J1+1))),sqrt(Siglam1[b-1]/(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%as.matrix(rep(1,J1+1)))))


          Siglam1[b]=1/rgamma(1,a1+(J1+1)/2,rate=b1+.5*(t(as.matrix(rep(Mulam1[b],J1+1))-L1)%*%solve(SigLam1)%*%(as.matrix(rep(Mulam1[b],J1+1))-L1)))


          ##Siglam

          iter[2]="Sigma"
        }else{



          Mulam1[b]=rnorm(1,lam1[b-1,1],sqrt(Siglam1[b-1]))


          Siglam1[b]=1/rgamma(1,a1+1/2,rate=b1+.5*(Mulam1[b]-lam1[b-1,1])^2)



        }

        #if(is.finite(Mulam1[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}
        #if(is.finite(Siglam1[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}



        #lambda1
        iter[2]="lam1"
        lam1[b,]=lam1[b-1,]
        #######

        for(m in 1:(J1+1)){



          lam=lam1[b,]
          lam=lam[is.na(lam)==FALSE]
          lambda=lam

          lam[m]=lambda[m]+runif(1,-cl1,cl1)


          if(J1==0){
            do=log(dnorm(lambda[m],Mulam1[b],sqrt(Siglam1[b])))
            dn=log(dnorm(lam[m],Mulam1[b],sqrt(Siglam1[b])))
          }else{


            #do=-(t(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])

            #dn=-(t(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])


            do=dmvnorm(lambda,rep(Mulam1[b],J1+1),Siglam1[b]*SigLam1)
            do=dmvnorm(lam,rep(Mulam1[b],J1+1),Siglam1[b]*SigLam1)
          }

          Likeo=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b-1,],lam3[b-1,],gam[b,])

          Liken=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam,lam2[b-1,],lam3[b-1,],gam[b,])



          U=log(runif(1,0,1))
          alphalam=Liken-Likeo+dn-do

          if(is.nan(alphalam)==TRUE){
            lam1[b,m]=lam1[b-1,m]
            Acceptlam1[b,m]=0
          }else{

            if(U<alphalam){
              Acceptlam1[b,m]=1
              lam1[b,m]=lam[m]
            }else{Acceptlam1[b,m]=0}
          }


        }




        ####Lam2####

        iter[1]="LogBH2"
        iter[2]="matrixsetup"

        W1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)
        Q1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)


        length1=diff(s2[b-1,])







        if(J2<2){
          if(J2==1){
            W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
            SigLam2=solve(diag(J2+1)-W1)%*%Q1





          }else{

            Q1=as.matrix(2/(m2))
            SigLam2=Q1
          }
        }else{


          for(j in 2:J2){
            W1[j,j-1]=(clam2*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            W1[j,j+1]=(clam2*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
          }


          Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
          Q1[1,1]=2/(2*length1[1]+length1[2])
          W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
          W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])


          SigLam2=solve(diag(J2+1)-W1)%*%Q1

        }

        iter[2]="Mu"
        ##Lambda1 Hierarchical Sampler
        ##Mulam
        if(J2>0){

          Mulam2[b]=rnorm(1,(t(as.matrix(rep(1,J2+1)))%*%solve(SigLam2)%*%L2)/(t(as.matrix(rep(1,J2+1)))%*%solve(SigLam2)%*%as.matrix(rep(1,J2+1))),sqrt(Siglam2[b-1]/(t(as.matrix(rep(1,J2+1)))%*%solve(SigLam2)%*%as.matrix(rep(1,J2+1)))))


          Siglam2[b]=1/rgamma(1,a2+(J2+1)/2,rate=b2+.5*(t(as.matrix(rep(Mulam2[b],J2+1))-L2)%*%solve(SigLam2)%*%(as.matrix(rep(Mulam2[b],J2+1))-L2)))


          ##Siglam
          iter[2]="Sigma"

        }else{



          Mulam2[b]=rnorm(1,lam2[b-1,1],sqrt(Siglam2[b-1]))

          Siglam2[b]=1/rgamma(1,a2+1/2,rate=b2+.5*(Mulam2[b]-lam2[b-1,1])^2)




        }

        #if(is.finite(Mulam2[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}
        #if(is.finite(Siglam2[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}
        #lambda1
        iter[2]="lam2"
        lam2[b,]=lam2[b-1,]
        #######


        for(m in 1:(J2+1)){



          lam=lam2[b,]
          lam=lam[is.na(lam)==FALSE]
          lambda=lam

          lam[m]=lambda[m]+runif(1,-cl2,cl2)






          if(J2==0){
            do=log(dnorm(lambda[m],Mulam2[b],sqrt(Siglam2[b])))
            dn=log(dnorm(lam[m],Mulam2[b],sqrt(Siglam2[b])))
          }else{


            #do=-(t(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])

            #dn=-(t(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])


            do=dmvnorm(lambda,rep(Mulam2[b],J2+1),Siglam2[b]*SigLam2)
            do=dmvnorm(lam,rep(Mulam2[b],J2+1),Siglam2[b]*SigLam2)
          }



          #do=-(t(as.matrix(lambda)-as.matrix(rep(Mulam2[b],J2+1)))%*%solve(SigLam2)%*%(as.matrix(lambda)-as.matrix(rep(Mulam2[b],J2+1))))/(2*Siglam2[b])

          #dn=-(t(as.matrix(lam)-as.matrix(rep(Mulam2[b],J2+1)))%*%solve(SigLam2)%*%(as.matrix(lam)-as.matrix(rep(Mulam2[b],J2+1))))/(2*Siglam2[b])


          Likeo=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam3[b-1,],gam[b,])

          Liken=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam,lam3[b-1,],gam[b,])



          U=log(runif(1,0,1))
          alphalam=Liken-Likeo+dn-do

          if(is.nan(alphalam)==TRUE){
            lam2[b,m]=lam2[b-1,m]
            Acceptlam2[b,m]=0
          }else{

            if(U<alphalam){
              Acceptlam2[b,m]=1
              lam2[b,m]=lam[m]
            }else{Acceptlam2[b,m]=0}
          }


        }


        ####Lam2####

        iter[1]="LogBH3"
        iter[2]="matrixsetup"

        W1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)
        Q1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)


        length1=diff(s3[b-1,])







        if(J3<2){
          if(J3==1){
            W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
            SigLam3=solve(diag(J3+1)-W1)%*%Q1





          }else{

            Q1=as.matrix(2/(m3))
            SigLam3=Q1
          }
        }else{


          for(j in 2:J3){
            W1[j,j-1]=(clam3*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            W1[j,j+1]=(clam3*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
          }


          Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
          Q1[1,1]=2/(2*length1[1]+length1[2])
          W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
          W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])


          SigLam3=solve(diag(J3+1)-W1)%*%Q1

        }

        iter[2]="Mu"
        ##Lambda1 Hierarchical Sampler
        ##Mulam

        if(J3>0){

          iter[2]="Sigma"


          Mulam3[b]=rnorm(1,(t(as.matrix(rep(1,J3+1)))%*%solve(SigLam3)%*%L3)/(t(as.matrix(rep(1,J3+1)))%*%solve(SigLam3)%*%as.matrix(rep(1,J3+1))),sqrt(Siglam3[b-1]/(t(as.matrix(rep(1,J3+1)))%*%solve(SigLam3)%*%as.matrix(rep(1,J3+1)))))
          ##Siglam

          Siglam3[b]=1/rgamma(1,a3+(J3+1)/2,rate=b3+.5*(t(as.matrix(rep(Mulam3[b],J3+1))-L3)%*%solve(SigLam3)%*%(as.matrix(rep(Mulam3[b],J3+1))-L3)))


        }else{



          Mulam3[b]=rnorm(1,lam3[b-1,1],sqrt(Siglam3[b-1]))

          Siglam3[b]=1/rgamma(1,a3+1/2,rate=b3+.5*(Mulam3[b]-lam3[b-1,1])^2)



        }

        #if(is.finite(Mulam3[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}
        #if(is.finite(Siglam3[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}

        #lambda3
        iter[2]="lam3"
        lam3[b,]=lam3[b-1,]
        #######

        for(m in 1:(J3+1)){


          lam=lam3[b,]
          lam=lam[is.na(lam)==FALSE]
          lambda=lam

          lam[m]=lambda[m]+runif(1,-cl3,cl3)



          if(J3==0){
            do=log(dnorm(lambda[m],Mulam3[b],sqrt(Siglam3[b])))
            dn=log(dnorm(lam[m],Mulam3[b],sqrt(Siglam3[b])))
          }else{


            #do=-(t(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])

            #dn=-(t(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])


            do=dmvnorm(lambda,rep(Mulam3[b],J3+1),Siglam3[b]*SigLam3)
            do=dmvnorm(lam,rep(Mulam3[b],J3+1),Siglam3[b]*SigLam3)
          }



          Likeo=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])

          Liken=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                     s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam,gam[b,])


          U=log(runif(1,0,1))
          alphalam=Liken-Likeo+dn-do

          if(is.nan(alphalam)==TRUE){
            lam3[b,m]=lam3[b-1,m]
            Acceptlam3[b,m]=0
          }else{

            if(U<alphalam){
              Acceptlam3[b,m]=1
              lam3[b,m]=lam[m]
            }else{Acceptlam3[b,m]=0}
          }


        }


        ##############################################
        ######## PUT BACK LAMBDA SAMPLERS HERE!!! ###

        ###Delete these later
        s2[b,]=s2[b-1,]
        s3[b,]=s3[b-1,]

        #####################################################
        ###################################################

        iter[1]="Haz1"
        iter[2]="Birth"

        ###Random Perturbation###
        U1=runif(1,0,1)
        #####

        s=s1[b-1,]
        s=s[!is.na(s)]

        if(length(s)<J1max){
          Birth=runif(1,0,m1)

          s1[b,1:(J1+3)]=sort(c(s,Birth))

          for(k in 2:(J1+2)){
            if(Birth>s1[b-1,k-1] & Birth<s1[b-1,k]){
              Ind=k-1
            }
          }

          lam=rep(0,J1+2)

          if(Ind==1 | Ind==J1+1){
            if(Ind==1){
              lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[(Ind+2):length(lam)]=lam1[b,(Ind+1):(J1+1)]
            }else{
              lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[1:(Ind-1)]=lam1[b,1:(Ind-1)]
            }
          }else{
            lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
            lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
            lam[1:(Ind-1)]=lam1[b,1:(Ind-1)]
            lam[(Ind+2):length(lam)]=lam1[b,(Ind+1):(J1+1)]
          }

          lam=lam[!is.na(lam)]

          lambda=lam1[b,]
          lambda=lambda[!is.na(lambda)]

          Lo=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b-1,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])
          if(J1>0){
            do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))
          }else{
            do=log(dpois(J1,alpha1))+log(dnorm(lambda,Mulam1[b],Siglam1[b]))
          }

          prior=((2*J1+3)*(2*J1+2)*(Birth-s1[b-1,Ind])*(s1[b-1,Ind+1]-Birth))/((m1^2)*(s1[b-1,Ind+1]-s1[b-1,Ind]))

          G1=G1+1
          J1=J1+1

          Ln=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b-1,],s3[b-1,],lam,lam2[b,],lam3[b,],gam[b,])


          ##Make SigLam1



          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


          length1=diff(s1[b,])


          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1n=solve(diag(J1+1)-W1)%*%Q1





            }else{

              SigLam1n=2/m1
            }
          }else{


            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1n=solve(diag(J1+1)-W1)%*%Q1

          }



          dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),Siglam1[b]*SigLam1n))




          alpha=Ln-Lo+dn-do-log(U1*(1-U1)) + log(prior)

          if(is.nan(alpha)==TRUE){
            IndB1[b]=0
            s1[b,]=s1[b-1,]
            J1=J1-1
            G1=G1-1
          }else{

            U=log(runif(1,0,1))

            if(U<alpha){
              IndB1[b]=1
              lam1[b,1:(J1+1)]=lam
            }else{
              s1[b,]=s1[b-1,]
              IndB1[b]=0
              J1=J1-1
              G1=G1-1
            }

          }


        }else{
          s1[b,]=s1[b-1,]
          IndB1[b]=0
        }


        #########################################################
        ###################Death Sampler#########################
        ##########################################################
        iter[2]="Death"

        U1=runif(1,0,1)

        if(J1==0){
          IndD1[b]=0
          s1[b,]=s1[b-1,]
        }else{

          if(J1==1){
            Ind=2
          }else{

            Ind=sample(2:(J1+1),1)
          }


          s=s1[b,]
          s=s[-Ind]

          lam=lam1[b,]
          lambda=lam[!is.na(lam)]

          lam=lam[!is.na(lam)]
          lam=lam[-Ind]

          lam[Ind-1]=((s1[b,Ind]-s1[b,Ind-1])*lam1[b,Ind-1]+(s1[b,Ind+1]-s1[b,Ind])*lam1[b,Ind])/(s1[b,Ind+1]-s1[b,Ind-1])



          #############################################
          ####Sets up SigLam1 matrix for old density###
          #############################################


          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


          length1=diff(s1[b,])



          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1=solve(diag(J1+1)-W1)%*%Q1

              do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))




            }else{


              do=log(dpois(J1,alpha1))+log(dnorm(lambda,Mulam1[b],Siglam1[b]))

            }
          }else{

            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1=solve(diag(J1+1)-W1)%*%Q1

            do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))


          }
          #############################################
          #############################################

          Lo=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])


          prior=((m1^2)*(s1[b,Ind+1]-s1[b,Ind-1]))/((2*J1+1)*(2*J1)*(s1[b,Ind]-s1[b,Ind-1])*(s1[b,Ind+1]-s1[b,Ind]))


          G1=G1-1
          J1=J1-1


          Ln=LK1L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s,s2[b-1,],s3[b-1,],lam,lam2[b,],lam3[b,],gam[b,])

          ###Make siglam matrix



          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)



          length1=rep(0,J1+1)




          for(j in 1:length(length1)){
            length1[j]=s[j+1]-s[j]
          }


          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1n=solve(diag(J1+1)-W1)%*%Q1


              dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),SigLam1n*Siglam1[b]))



            }else{

              SigLam1n=2/m1
              dn=log(dpois(J1,alpha1))+log(dnorm(lam,Mulam1[b],Siglam1[b]))

            }
          }else{


            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1n=solve(diag(J1+1)-W1)%*%Q1

            dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),SigLam1n*Siglam1[b]))


          }
          ####




          alpha=Ln-Lo+dn-do+log(prior)+log(U1*(1-U1))

          if(is.nan(alpha)==TRUE){
            IndD1[b]=0
            J1=J1+1
            G1=G1+1
          }else{

            U=log(runif(1,0,1))

            iter[2]="AcceptRejDeath"

            if(U<alpha){
              s1[b,]=c(s,NA)
              IndD1[b]=1
              lam1[b,1:(J1+1)]=lam
              lam1[b,(J1+2):J1max]=rep(NA,J1max-J1-1)
            }else{
              IndD1[b]=0
              J1=J1+1
              G1=G1+1
            }
          }

          ####End else
        }
        ##






        #######################
        #####End of Death sampler
        ######################


        #####################################################
        ###################################################

        iter[1]="Haz2"
        iter[2]="Birth"

        ###Random Perturbation###
        U2=runif(1,0,1)
        #####

        s=s2[b-1,]
        s=s[!is.na(s)]

        if(length(s)<J2max){
          Birth=runif(1,0,m2)

          s2[b,1:(J2+3)]=sort(c(s,Birth))

          for(k in 2:(J2+2)){
            if(Birth>s2[b-1,k-1] & Birth<s2[b-1,k]){
              Ind=k-1
            }
          }

          lam=rep(0,J2+2)

          if(Ind==1 | Ind==J2+1){
            if(Ind==1){
              lam[Ind]=lam2[b,Ind] - ((s2[b-1,Ind+1]-Birth)/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
              lam[Ind+1]=lam2[b,Ind] + ((Birth-s2[b-1,Ind])/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
              lam[(Ind+2):length(lam)]=lam2[b,(Ind+1):(J2+1)]
            }else{
              lam[Ind]=lam2[b,Ind] - ((s2[b-1,Ind+1]-Birth)/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
              lam[Ind+1]=lam2[b,Ind] + ((Birth-s2[b-1,Ind])/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
              lam[1:(Ind-1)]=lam2[b,1:(Ind-1)]
            }
          }else{
            lam[Ind]=lam2[b,Ind] - ((s2[b-1,Ind+1]-Birth)/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
            lam[Ind+1]=lam2[b,Ind] + ((Birth-s2[b-1,Ind])/(s2[b-1,Ind+1]-s2[b-1,Ind]))*log((1-U2)/U2)
            lam[1:(Ind-1)]=lam2[b,1:(Ind-1)]
            lam[(Ind+2):length(lam)]=lam2[b,(Ind+1):(J2+1)]
          }

          lam=lam[!is.na(lam)]

          lambda=lam2[b,]
          lambda=lambda[!is.na(lambda)]

          Lo=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b-1,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])
          if(J2>0){
            do=log(dpois(J2,alpha2))+log(dmvnorm(lambda,rep(Mulam2[b],length(lambda)),SigLam2*Siglam2[b]))
          }else{
            do=log(dpois(J2,alpha2))+log(dnorm(lambda,Mulam2[b],Siglam2[b]))
          }

          prior=((2*J2+3)*(2*J2+2)*(Birth-s2[b-1,Ind])*(s2[b-1,Ind+1]-Birth))/((m2^2)*(s2[b-1,Ind+1]-s2[b-1,Ind]))

          G2=G2+1
          J2=J2+1

          Ln=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s3[b-1,],lam1[b,],lam,lam3[b,],gam[b,])


          ##Make SigLam1



          W1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)
          Q1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)


          length1=diff(s2[b,])




          if(J2<2){
            if(J2==1){
              W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
              SigLam2n=solve(diag(J2+1)-W1)%*%Q1





            }else{

              SigLam2n=2/m2
            }
          }else{


            for(j in 2:J2){
              W1[j,j-1]=(clam2*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam2*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])


            SigLam2n=solve(diag(J2+1)-W1)%*%Q1

          }



          dn=log(dpois(J2,alpha2))+log(dmvnorm(lam,rep(Mulam2[b],length(lam)),Siglam2[b]*SigLam2n))




          alpha=Ln-Lo+dn-do-log(U2*(1-U2)) + log(prior)

          if(is.nan(alpha)==TRUE){
            IndB2[b]=0
            s2[b,]=s2[b-1,]
            J2=J2-1
            G2=G2-1
          }else{

            U=log(runif(1,0,1))

            if(U<alpha){
              IndB2[b]=1
              lam2[b,1:(J2+1)]=lam
            }else{
              s2[b,]=s2[b-1,]
              IndB2[b]=0
              J2=J2-1
              G2=G2-1
            }

          }


        }else{
          s2[b,]=s2[b-1,]
          IndB2[b]=0
        }


        #########################################################
        ###################Death Sampler#########################
        ##########################################################
        iter[2]="Death"

        U2=runif(1,0,1)
        if(J2==0){
          IndD2[b]=0
          s2[b,]=s2[b-1,]
        }else{

          if(J2==1){
            Ind=2
          }else{

            Ind=sample(2:(J2+1),1)
          }


          s=s2[b,]
          s=s[-Ind]

          lam=lam2[b,]
          lambda=lam[!is.na(lam)]

          lam=lam[!is.na(lam)]
          lam=lam[-Ind]

          lam[Ind-1]=((s2[b,Ind]-s2[b,Ind-1])*lam2[b,Ind-1]+(s2[b,Ind+1]-s2[b,Ind])*lam2[b,Ind])/(s2[b,Ind+1]-s2[b,Ind-1])



          #############################################
          ####Sets up SigLam1 matrix for old density###
          #############################################


          W1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)
          Q1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)


          length1=diff(s2[b,])







          if(J2<2){
            if(J2==1){
              W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
              SigLam2=solve(diag(J2+1)-W1)%*%Q1

              do=log(dpois(J2,alpha2))+log(dmvnorm(lambda,rep(Mulam2[b],length(lambda)),SigLam2*Siglam2[b]))




            }else{

              SigLam2=2/m2

              do=log(dpois(J2,alpha2))+log(dnorm(lambda,Mulam2[b],Siglam2[b]))

            }
          }else{

            for(j in 2:J2){
              W1[j,j-1]=(clam2*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam2*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])


            SigLam2=solve(diag(J2+1)-W1)%*%Q1

            do=log(dpois(J2,alpha2))+log(dmvnorm(lambda,rep(Mulam2[b],length(lambda)),SigLam2*Siglam2[b]))


          }
          #############################################
          #############################################

          Lo=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])


          prior=((m2^2)*(s2[b,Ind+1]-s2[b,Ind-1]))/((2*J2+1)*(2*J2)*(s2[b,Ind]-s2[b,Ind-1])*(s2[b,Ind+1]-s2[b,Ind]))


          G2=G2-1
          J2=J2-1


          Ln=LK2L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s,s3[b-1,],lam1[b,],lam,lam3[b,],gam[b,])

          ###Make siglam matrix



          W1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)
          Q1=matrix(rep(0,(J2+1)*(J2+1)),nrow=J2+1)



          length1=rep(0,J2+1)




          for(j in 1:length(length1)){
            length1[j]=s[j+1]-s[j]
          }


          if(J2<2){
            if(J2==1){
              W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
              SigLam2n=solve(diag(J2+1)-W1)%*%Q1


              dn=log(dpois(J2,alpha2))+log(dmvnorm(lam,rep(Mulam2[b],length(lam)),SigLam2n*Siglam2[b]))



            }else{

              dn=log(dpois(J2,alpha2))+log(dnorm(lam,Mulam2[b],Siglam2[b]))

            }
          }else{


            for(j in 2:J2){
              W1[j,j-1]=(clam2*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam2*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J2+1,J2+1]=2/(length1[J2]+2*length1[J2+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam2*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J2+1,J2]=(clam2*(length1[J2+1]+length1[J2]))/(length1[J2]+2*length1[J2+1])


            SigLam2n=solve(diag(J2+1)-W1)%*%Q1

            dn=log(dpois(J2,alpha2))+log(dmvnorm(lam,rep(Mulam2[b],length(lam)),SigLam2n*Siglam2[b]))


          }
          ####




          alpha=Ln-Lo+dn-do+log(prior)+log(U2*(1-U2))

          if(is.nan(alpha)==TRUE){
            IndD2[b]=0
            J2=J2+1
            G2=G2+1
          }else{

            U=log(runif(1,0,1))

            iter[2]="AcceptRejDeath"

            if(U<alpha){
              s2[b,]=c(s,NA)
              IndD2[b]=1
              lam2[b,1:(J2+1)]=lam
              lam2[b,(J2+2):J2max]=rep(NA,J2max-J2-1)
            }else{
              IndD2[b]=0
              J2=J2+1
              G2=G2+1
            }
          }

          ####End else
        }
        ##






        #######################
        #####End of Death sampler
        ######################


        #####################################################
        ###################################################

        iter[1]="Haz3"
        iter[2]="Birth"

        ###Random Perturbation###
        U3=runif(1,0,1)
        #####

        s=s3[b-1,]
        s=s[!is.na(s)]

        if(length(s)<J3max){
          Birth=runif(1,0,m3)

          s3[b,1:(J3+3)]=sort(c(s,Birth))

          for(k in 2:(J3+2)){
            if(Birth>s3[b-1,k-1] & Birth<s3[b-1,k]){
              Ind=k-1
            }
          }

          lam=rep(0,J3+2)

          if(Ind==1 | Ind==J3+1){
            if(Ind==1){
              lam[Ind]=lam3[b,Ind] - ((s3[b-1,Ind+1]-Birth)/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
              lam[Ind+1]=lam3[b,Ind] + ((Birth-s3[b-1,Ind])/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
              lam[(Ind+2):length(lam)]=lam3[b,(Ind+1):(J3+1)]
            }else{
              lam[Ind]=lam3[b,Ind] - ((s3[b-1,Ind+1]-Birth)/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
              lam[Ind+1]=lam3[b,Ind] + ((Birth-s3[b-1,Ind])/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
              lam[1:(Ind-1)]=lam3[b,1:(Ind-1)]
            }
          }else{
            lam[Ind]=lam3[b,Ind] - ((s3[b-1,Ind+1]-Birth)/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
            lam[Ind+1]=lam3[b,Ind] + ((Birth-s3[b-1,Ind])/(s3[b-1,Ind+1]-s3[b-1,Ind]))*log((1-U3)/U3)
            lam[1:(Ind-1)]=lam3[b,1:(Ind-1)]
            lam[(Ind+2):length(lam)]=lam3[b,(Ind+1):(J3+1)]
          }

          lam=lam[!is.na(lam)]

          lambda=lam3[b,]
          lambda=lambda[!is.na(lambda)]

          Lo=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s3[b-1,],lam1[b,],lam2[b,],lam3[b,],gam[b,])
          if(J3>0){
            do=log(dpois(J3,alpha3))+log(dmvnorm(lambda,rep(Mulam3[b],length(lambda)),SigLam3*Siglam3[b]))
          }else{
            do=log(dpois(J3,alpha3))+log(dnorm(lambda,Mulam3[b],Siglam3[b]))
          }

          prior=((2*J3+3)*(2*J3+2)*(Birth-s3[b-1,Ind])*(s3[b-1,Ind+1]-Birth))/((m3^2)*(s3[b-1,Ind+1]-s3[b-1,Ind]))

          G3=G3+1
          J3=J3+1

          Ln=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s3[b,],lam1[b,],lam2[b,],lam,gam[b,])


          ##Make SigLam1



          W1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)
          Q1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)


          length1=diff(s3[b,])



          if(J3<2){
            if(J3==1){
              W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
              SigLam3n=solve(diag(J3+1)-W1)%*%Q1





            }else{

              SigLam3n=2/m3
            }
          }else{


            for(j in 2:J3){
              W1[j,j-1]=(clam3*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam3*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])


            SigLam3n=solve(diag(J3+1)-W1)%*%Q1

          }



          dn=log(dpois(J3,alpha3))+log(dmvnorm(lam,rep(Mulam3[b],length(lam)),Siglam3[b]*SigLam3n))




          alpha=Ln-Lo+dn-do-log(U3*(1-U3)) + log(prior)

          if(is.nan(alpha)==TRUE){
            IndB3[b]=0
            s3[b,]=s3[b-1,]
            J3=J3-1
            G3=G3-1
          }else{

            U=log(runif(1,0,1))

            if(U<alpha){
              IndB3[b]=1
              lam3[b,1:(J3+1)]=lam
            }else{
              s3[b,]=s3[b-1,]
              IndB3[b]=0
              J3=J3-1
              G3=G3-1
            }

          }


        }else{
          s3[b,]=s3[b-1,]
          IndB3[b]=0
        }


        #########################################################
        ###################Death Sampler#########################
        ##########################################################
        iter[2]="Death"

        U3=runif(1,0,1)

        if(J3==0){
          IndD3[b]=0
          s3[b,]=s3[b-1,]
        }else{

          if(J3==1){
            Ind=2
          }else{

            Ind=sample(2:(J3+1),1)
          }


          s=s3[b,]
          s=s[-Ind]

          lam=lam3[b,]
          lambda=lam[!is.na(lam)]

          lam=lam[!is.na(lam)]
          lam=lam[-Ind]

          lam[Ind-1]=((s3[b,Ind]-s3[b,Ind-1])*lam3[b,Ind-1]+(s3[b,Ind+1]-s3[b,Ind])*lam3[b,Ind])/(s3[b,Ind+1]-s3[b,Ind-1])



          #############################################
          ####Sets up SigLam1 matrix for old density###
          #############################################


          W1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)
          Q1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)


          length1=diff(s3[b,])



          if(J3<2){
            if(J3==1){
              W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
              SigLam3=solve(diag(J3+1)-W1)%*%Q1

              do=log(dpois(J3,alpha3))+log(dmvnorm(lambda,rep(Mulam3[b],length(lambda)),SigLam3*Siglam3[b]))




            }else{


              do=log(dpois(J3,alpha3))+log(dnorm(lambda,Mulam3[b],Siglam3[b]))

            }
          }else{

            for(j in 2:J3){
              W1[j,j-1]=(clam3*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam3*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])


            SigLam3=solve(diag(J3+1)-W1)%*%Q1

            do=log(dpois(J3,alpha3))+log(dmvnorm(lambda,rep(Mulam3[b],length(lambda)),SigLam3*Siglam3[b]))


          }
          #############################################
          #############################################

          Lo=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s3[b,],lam1[b,],lam2[b,],lam3[b,],gam[b,])


          prior=((m3^2)*(s3[b,Ind+1]-s3[b,Ind-1]))/((2*J3+1)*(2*J3)*(s3[b,Ind]-s3[b,Ind-1])*(s3[b,Ind+1]-s3[b,Ind]))


          G3=G3-1
          J3=J3-1


          Ln=LK3L(Y1,Y2,I1,I2,X,as.matrix(beta1[b,]),as.matrix(beta2[b,]),as.matrix(beta3[b,]),
                  s1[b,],s2[b,],s,lam1[b,],lam2[b,],lam,gam[b,])

          ###Make siglam matrix



          W1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)
          Q1=matrix(rep(0,(J3+1)*(J3+1)),nrow=J3+1)



          length1=rep(0,J3+1)




          for(j in 1:length(length1)){
            length1[j]=s[j+1]-s[j]
          }


          if(J3<2){
            if(J3==1){
              W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
              SigLam3n=solve(diag(J3+1)-W1)%*%Q1


              dn=log(dpois(J3,alpha3))+log(dmvnorm(lam,rep(Mulam3[b],length(lam)),SigLam3n*Siglam3[b]))



            }else{

              dn=log(dpois(J3,alpha3))+log(dnorm(lam,Mulam3[b],Siglam3[b]))

            }
          }else{


            for(j in 2:J3){
              W1[j,j-1]=(clam3*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam3*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J3+1,J3+1]=2/(length1[J3]+2*length1[J3+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam3*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J3+1,J3]=(clam3*(length1[J3+1]+length1[J3]))/(length1[J3]+2*length1[J3+1])


            SigLam3n=solve(diag(J3+1)-W1)%*%Q1

            dn=log(dpois(J3,alpha3))+log(dmvnorm(lam,rep(Mulam3[b],length(lam)),SigLam3n*Siglam3[b]))


          }
          ####




          alpha=Ln-Lo+dn-do+log(prior)+log(U3*(1-U3))

          if(is.nan(alpha)==TRUE){
            IndD3[b]=0
            J3=J3+1
            G3=G3+1
          }else{

            U=log(runif(1,0,1))

            iter[2]="AcceptRejDeath"

            if(U<alpha){
              s3[b,]=c(s,NA)
              IndD3[b]=1
              lam3[b,1:(J3+1)]=lam
              lam3[b,(J3+2):J3max]=rep(NA,J3max-J3-1)
            }else{
              IndD3[b]=0
              J3=J3+1
              G3=G3+1
            }
          }

          ####End else
        }
        split1[b]=J1
        split2[b]=J2
        split3[b]=J3
        ##
        sum1[b]=sum(eta1[b,])
        sum2[b]=sum(eta2[b,])
        sum3[b]=sum(eta3[b,])


        ##End Sampler
      }
      ###End of Sampler




      ################End Samplers
      cat("

          ", "Posterior Inclusion Probabilities after half Burnin", "


          ", "Hazard 1", "

          ", colMeans(eta1[(B*burn+1):B,])*100, "

          ", "Hazard 2", "

          ", colMeans(eta2[(B*burn+1):B,])*100, "

          ", "Hazard 3", "

          ", colMeans(eta3[(B*burn+1):B,])*100,"

          ", "IndEta",mean(Indeta1[(B*burn+1):B])*100,mean(Indeta2[(B*burn+1):B])*100,mean(Indeta3[(B*burn+1):B])*100,"


          ","IndMix",mean(Indmix1[(B*burn+1):B])*100,mean(Indmix2[(B*burn+1):B])*100,mean(Indmix3[(B*burn+1):B])*100,"


          ", "Included Acceptance", "

          ", "Haz1", "
          ", "

          ", colMeans(Indcond1[(B*burn+1):B,],na.rm=TRUE)*100,"


          ",  "Haz2", "
          ", "

          ", colMeans(Indcond2[(B*burn+1):B,],na.rm=TRUE)*100,"



          ", "Haz3","



          ", colMeans(Indcond1[(B*burn+1):B,],na.rm=TRUE)*100,"


          ","Survival","

          ","IndDeath",mean(IndD1[(B*burn+1):B])*100,mean(IndD2[(B*burn+1):B])*100,mean(IndD3[(B*burn+1):B])*100,"

          ","IndBirth",mean(IndB1[(B*burn+1):B])*100,mean(IndB2[(B*burn+1):B])*100,mean(IndB3[(B*burn+1):B])*100,"

          ","Lambda","

          ", "Lam1",
          colMeans(Acceptlam1[(B*burn+1):B,],na.rm=TRUE)*100,"

          ","Lam2",
          colMeans(Acceptlam2[(B*burn+1):B,],na.rm=TRUE)*100,"

          ","Lam3",
          colMeans(Acceptlam3[(B*burn+1):B,],na.rm=TRUE)*100,"

          ","Indepsilon",mean(Indepsilon[(B*burn+1):B])*100)




      ##End


    }





    ###Return Values



    Path1= paste0(Path,"/beta1.txt")

    write.table(beta1[(burn*B+1):B,], Path1, sep="\t")

    Path1= paste0(Path,"/beta2.txt")

    write.table(beta2[(burn*B+1):B,], Path1, sep="\t")

    Path1= paste0(Path,"/beta3.txt")

    write.table(beta3[(burn*B+1):B,], Path1, sep="\t")

    Path1= paste0(Path,"/eta1.txt")

    write.table(eta1[(burn*B+1):B,], Path1, sep="\t")

    Path1= paste0(Path,"/eta2.txt")

    write.table(eta2[(burn*B+1):B,], Path1, sep="\t")

    Path1= paste0(Path,"/eta3.txt")

    write.table(eta3[(burn*B+1):B,], Path1, sep="\t")

    Path1= paste0(Path,"/gam.txt")

    write.table(gam[(burn*B+1):B,], Path1, sep="\t")

    Path1= paste0(Path,"/lam1.txt")

    write.table(lam1[(burn*B+1):B,], Path1, sep="\t")

    Path1= paste0(Path,"/lam2.txt")

    write.table(lam2[(burn*B+1):B,], Path1, sep="\t")

    Path1= paste0(Path,"/lam3.txt")

    write.table(lam3[(burn*B+1):B,], Path1, sep="\t")

    Path1= paste0(Path,"/s1.txt")

    write.table(s1[(burn*B+1):B,], Path1, sep="\t")

    Path1= paste0(Path,"/s2.txt")

    write.table(s2[(burn*B+1):B,], Path1, sep="\t")


    Path1= paste0(Path,"/s3.txt")

    write.table(s3[(burn*B+1):B,], Path1, sep="\t")

    Path1= paste0(Path,"/sum1.txt")

    write.table(sum1[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/sum2.txt")

    write.table(sum2[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/sum3.txt")

    write.table(sum3[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/split1.txt")

    write.table(split1[(burn*B+1):B], Path1, sep="\t")


    Path1= paste0(Path,"/split2.txt")

    write.table(split2[(burn*B+1):B], Path1, sep="\t")


    Path1= paste0(Path,"/split3.txt")

    write.table(split3[(burn*B+1):B], Path1, sep="\t")
    ###########

    Path1= paste0(Path,"/siglam1.txt")

    write.table(Siglam1[(burn*B+1):B], Path1, sep="\t")


    Path1= paste0(Path,"/siglam2.txt")

    write.table(Siglam2[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/siglam3.txt")

    write.table(Siglam3[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/mulam1.txt")

    write.table(Mulam1[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/mulam2.txt")

    write.table(Mulam2[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/mulam3.txt")

    write.table(Mulam3[(burn*B+1):B], Path1, sep="\t")


    Path1= paste0(Path,"/epsilon.txt")

    write.table(epsilon[(burn*B+1):B], Path1, sep="\t")

    ##Acceptances

    Path1= paste0(Path,"/Indeta1.txt")

    write.table(Indeta1[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/Indeta2.txt")

    write.table(Indeta2[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/Indeta3.txt")

    write.table(Indeta3[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/IndD1.txt")

    write.table(IndD1[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/IndB1.txt")

    write.table(IndB1[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/IndD2.txt")

    write.table(IndD2[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/IndD3.txt")

    write.table(IndD3[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/IndB2.txt")

    write.table(IndB2[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/IndB3.txt")

    write.table(IndB3[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/Acceptlam1.txt")

    write.table(Acceptlam1[(burn*B+1):B,], Path1, sep="\t")

    Path1= paste0(Path,"/Acceptlam2.txt")

    write.table(Acceptlam2[(burn*B+1):B,], Path1, sep="\t")

    Path1= paste0(Path,"/Acceptlam3.txt")

    write.table(Acceptlam3[(burn*B+1):B,], Path1, sep="\t")

    Path1= paste0(Path,"/Indepsilon.txt")


    write.table(Indepsilon[(burn*B+1):B], Path1, sep="\t")






    Path1= paste0(Path,"/Indmix1.txt")


    write.table(Indmix1[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/Indmix2.txt")


    write.table(Indmix2[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/Indmix3.txt")


    write.table(Indmix3[(burn*B+1):B], Path1, sep="\t")


    Path1= paste0(Path,"/meangam.txt")


    write.table(colMeans(gam[(burn*B+1):B,]), Path1, sep="\t")




    par(mfrow=c(3,1))

    plot(1:B,sum1,type="l",xlab="",ylab="Haz1: # Included", main="Traceplot: # Included")
    plot(1:B,sum2,type="l",xlab="",ylab="Haz2: # Included")
    plot(1:B,sum3,type="l",xlab="",ylab="Haz3: # Included")


    plot(1:B,split1,type="l",xlab="",ylab="Haz1: Split #", main="Traceplot: # Split points")
    plot(1:B,split2,type="l",xlab="",ylab="Haz2: Split #")
    plot(1:B,split3,type="l",xlab="",ylab="Haz3: Split #")





  }}
