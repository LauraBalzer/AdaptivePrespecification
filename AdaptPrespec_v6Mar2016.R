
#####################################
# Simulations to illustrate estimation and inference for the PATE
#  and SATE with the unadjusted estimator, the MLE with a priori specified
#	adjustment set (W.9), TMLE with adaptive pre-specification for initial 
# estimation of Qbar(A,W) and C-TMLE including collaborative estimation of g(A|W)
#
# Eq.# + Sec refer to the paper (Adaptive Pre-specification in Randomized
#	Trials with and without Pair-Matching)
#
# Programmer: Laura Balzer (lbbalzer@gmail.com)
# 
# Last update: 03.06.2016 
#	for a more general data structure (under pair-matching)
# 	for more general CV-scheme (previous was only LOOCV)

# For previous versions - please email
######################################


#------------------------------------
# simulate.Data.and.Run: function to generate the simulated data
#	and run the estimators 
#
# input: 
# output: true values of PATE and SATE, point estimates and inference from the 4
# estimators, results from CV-selection
#------------------------------------

simulate.Data.and.Run<- function(){
  
  X<- get.Data(n)
  
  # SATE is the average of the difference in counterfactuals for the sample
  SATE<- mean(X$Y1- X$Y0) 
  
  # True value of the target parameter depending on the scientific question
  truth<- ifelse(POP.EFFECT, PATE, SATE)
  
  # If Matching, use the nbpMatching package function to create pairs
  if(PAIRED){
    
    # see nbpMatching package for further details
    dist<- distancematrix(gendistance(data.frame(X[, matchOn])))
    matches<- nonbimatch(dist)
    # matches contains ids for the pair match & the distance measure
    grpA<- as.numeric(matches$halves[,'Group1.Row'])
    grpB<- as.numeric(matches$halves[,'Group2.Row'])
    
    # re-organize the data into paired data & randomize the txt within pairs
    X1<- data.frame(X[grpA, ], Pair=1:n.pairs, A= rbinom(n.pairs, 1, .5))
    X2<- data.frame(X[grpB, ], Pair=1:n.pairs, A= ifelse(X1$A==1, 0, 1 ))
    
    X.all<- NULL
    for(i in 1:n.pairs){
      X.all<- rbind(X.all, X1[i,], X2[i,]) 
    }
    
  } else{ # no matching
    
    # randomly assign treatment so that it's balanced overall
    A<- rbinom(n.pairs, 1, 0.5)
    A.2<- ifelse(A==1, 0, 1)
    # also include a dummy variable for pair
    X.all <- data.frame(X, Pair=1:n, A=c(A, A.2))
  }
  #---------------------
  # Generate observed data
  #---------------------
  
  # we observe the counterfactual outcome corresponding to 
  #	the observed exposure 
  Y<- ifelse(X.all$A, X.all$Y1, X.all$Y0)	
  
  # observed data are (W, A, Y)
  # keep the covariates in the adjustment set.
  data<- data.frame(X.all[,ADJ.SET], Pair=X.all$Pair, A=X.all$A, Y)
  
  #--------------------------
  # ESTIMATION AND INFERENCE
  #--------------------------
  
  # This is programmed so that the "Pair" represents the indpt unit
  # 	for a non-matched trial, Pair is a dummy variable for each observation
  
  # We are doing LOO-CV with each Fold="Pair"
  # For other CV-schemes, change the specification of folds here
  data<- data.frame(data, folds=data$Pair)
  
  # unadjusted estimator 
  unadj<- do.Estimation(Do.CV.Inference=F, truth=truth, data=data, QAdj='U', 
                        family=FAMILY)
  
  # a priori spec. working model for Qbar(A,W), adjusting for W9
  mle.W9<- do.Estimation(Do.CV.Inference=F, truth=truth, data=data, QAdj='W.9',
                         family=FAMILY)
  
  # Adaptive Pre-Specified Approach for Step 1. Initial Estimation
  select.Q<- suppressWarnings( CV.selector(data, family=FAMILY, forQ=T, 
                                           gAdj=NULL))
  
  # TMLE based on this selection (need to do CV-inference) 
  tmle<- do.Estimation(Do.CV.Inference=T, truth=truth, data=data, QAdj=select.Q, 
                       family=FAMILY)
  
  # Adaptive Pre-Specified Approach for Step 2. Targeting
  select.G<- CV.selector(data,  family=FAMILY, forQ=F, QAdj=select.Q)
  
  # C-TMLE using the adaptive selection for Step1 and Step2
  ctmle<-  do.Estimation(Do.CV.Inference=T, truth=truth, data=data,
                         QAdj=select.Q, family=FAMILY, gAdj=select.G)
  
  CV.select<- data.frame(select.Q, select.G)
  
  RETURN<- list(PATE=PATE, SATE=SATE, unadj=unadj, mle.W9=mle.W9, tmle=tmle, 
                ctmle=ctmle, CV.select=CV.select)	
  
  RETURN
}


#------------------
# do.Estimation: function to do TMLE + get inference
#
# input: Do.CV.Inference (get cross-validated estimate of IC), 
#   truth (true value of the target parameter), obs. data, 
#   QAdj (adjustment variable for Qbar(A,W)), family for Qbar(A,W), 
#		gAdj (adjustment variable for g(A|W)), significance level (alpha)
#
# output: point estimate, inference, CV-inference (if appropriate)
#-----------------

do.Estimation<- function(Do.CV.Inference=F, truth, data, QAdj, family, 
                         gAdj=NULL, alpha=0.05){
  
  # get a point estimate with TMLE based on the full data (i.e. without CV)
  est<- suppressWarnings(do.TMLE.with.CV(Do.CV=F, train=data, QAdj=QAdj, 
                                         family=family, gAdj=gAdj))
  
  # estimate the influence curve (and residuals)
  resid<- (data$Y - est$QbarAW)
  DY<- est$H.AW*(data$Y - est$QbarAW)
  DW<- est$Qbar1W - est$Qbar0W - est$psi.hat
  
  # get the variance of the influence curve (appropriate for the study design)
  var.IC<- get.Variance.IC(resid=resid, DY=DY, DW=DW, indpt=data$Pair)
  
  # degrees of freedom for t-distribution
  df<- ifelse(PAIRED, (n.pairs-1),  (n-2) )	
  
  # get confidence interval coverage and rejection
  if(POP.EFFECT){
    inference<-get.Inference.tdist(truth=truth, psi.hat=est$psi.hat, alpha=alpha, 
                                df=df, var= var.IC$var.PATE)
  } else{
    inference<-get.Inference.tdist(truth=truth, psi.hat=est$psi.hat, alpha=alpha, 
                                df=df, var= var.IC$var.SATE)
  }
  
  #------------------------------------
  # IF GETTING A CROSS-VALIDATED ESTIMATE OF THE INFLUENCE CURVE
  #		(Section 6 of manuscript)
  #-------------------------------------
  if(Do.CV.Inference){  
    
    var.IC.CV<- get.CV.inference(data=data, QAdj=QAdj, family=family, gAdj=gAdj, get.variance=T)
    
    if(POP.EFFECT){
      inference.CV<- get.Inference.tdist(truth=truth, psi.hat=est$psi.hat, 
                                         alpha=alpha, df=df, var= var.IC.CV$var.PATE)
    }else{
      inference.CV<- get.Inference.tdist(truth=truth, psi.hat=est$psi.hat,
                                         alpha=alpha, df=df, var= var.IC.CV$var.SATE)
    }	
    
  } else{ # otherwise, fill in the CV-inference with NA
    
    inference.CV<- data.frame(var=NA, tstat=NA, cover=NA, reject=NA)
  }
  #------------------------------------
  
  RETURN<- data.frame(est$psi.hat, 	inference, inference.CV)
  
  colnames(RETURN)<- c('psi.hat','var','tstat','cover','reject', 'var.cv', 
                       'tstat.cv', 'cover.cv', 'reject.cv')		
  RETURN
}

#---------------------------
# do.TMLE.with.CV: runs the standard TMLE algorithm (i.e. training set=data); 
#   will also fit the TMLE algorithm for Qbar*(A,W) on training set 
#   and get estimates for validation set
#
# input: Do.CV (whether not doing cross-validation), train (training set), 
#   valid (validation set), Qadj (adjustment variable for Qbar(A,W)), 
#   family for fitting Qbar(A,W)), gAdj (adjustment variable for g(A|W))
#
# output: estimates of the propensity score (exposure mechanism), 
#   clever covariate (H.AW), targeted predictions on obs exposure, txt & control,
#   and point estimate 
#----------------------

do.TMLE.with.CV<- function(Do.CV=F, train, valid, QAdj, family, gAdj=NULL){	
  
  # do entire TMLE algorithm on the training set	
  train.temp<- train[, c(QAdj, 'A', 'Y')]
  
  # STEP 1: INITIAL ESTIMATION
  # fit the working model for Qbar(A,W) on the training set
  glm.out<-  glm( Y~., family=family,  data=train.temp ) 
  
  # get the predicted outcomes under obs exp, txt and control
  X1<- X0<- train
  X1$A<-1; X0$A<- 0	
  
  QbarAW.train <- predict(glm.out, newdata=train, type='response')
  Qbar1W.train<-  predict(glm.out, newdata=X1, type='response')
  Qbar0W.train<-  predict(glm.out, newdata=X0, type='response')
  # ignore any warning messages
  
  #-------------------
  # STEP 2: TARGETING 
  # estimate the propensity score (i.e. the exposure mechanism)
  # code assumes 2 armed randomized trial with equal allocation g(A|W)=0.5
  # O/W change the following to be true prob of being randomized to the txt
  
  if( is.null(gAdj) ){
    # if adj var for g(A|W) is null, then g(A|W)=0.5
    pscore.train <-  rep(0.5, nrow(train) )
  } else if (gAdj=='U'){	
    # if adj var for g(A|W) is 'U' (unadj), then g(A|W)=0.5
    pscore.train <-  rep(0.5, nrow(train) )
  } else{  
    # otherwise fit the working model for the pscore on training set 	
    train.temp<- train[, c(gAdj, 'A')]
    p.out<-   glm( A~., family='binomial', data=train.temp) 
    pscore.train <- predict(p.out, newdata=train,  type="response")
  }
  
  # calculate the clever covariate
  H.AW.train<- train$A/pscore.train - (1-train$A)/(1-pscore.train)
  
  # update the initial estimator of the outcome regression
  
  if(family!='binomial'){  
    # if the outcome Y is continuous & unbounded outcome, do a linear update
    linUpdate<-  glm(train$Y ~ -1 +offset( QbarAW.train) + H.AW.train, 
                     family="gaussian")
    eps<-linUpdate$coef
    
    # updated QbarAW estimates for training set. 
    QbarAW.train<- QbarAW.train + eps*H.AW.train	
    Qbar1W.train<- Qbar1W.train + eps/pscore.train
    Qbar0W.train<- Qbar0W.train - eps/(1-pscore.train)	
    
  }else { 
    # if binary or bounded outcome in [0,1], logistic update
    # See Gruber and van der Laan 2010 for further details
    
    logitUpdate<- suppressWarnings( glm(train$Y ~ -1 +offset(qlogis(QbarAW.train))
                                        + H.AW.train, family="binomial"))
    eps<-logitUpdate$coef
    
    # updated QbarAW estimates for training set. 
    QbarAW.train<- plogis( qlogis(QbarAW.train) + eps*H.AW.train)	
    Qbar1W.train<- plogis( qlogis(Qbar1W.train) + eps/pscore.train)
    Qbar0W.train<- plogis( qlogis(Qbar0W.train) - eps/(1-pscore.train))
  } 		
  
  # STEP 3: PARAMETER ESTIMATION
  psi.hat<- mean(Qbar1W.train- Qbar0W.train)
  
  #-----------------------
  # if not doing cross-validation (i.e. if running the std TMLE algorithm)	
  if(!Do.CV){
    
    
    RETURN<- list(pscore=pscore.train, H.AW=H.AW.train, QbarAW=QbarAW.train, 
                  Qbar1W=Qbar1W.train, Qbar0W=Qbar0W.train, psi.hat=psi.hat)
    
  } else{  
    #----------------------------------------------------------
    # If want estimates for the validation set 
    # (for data-adaptive selection or for a CV-variance estimate)
    
    # get initial estimates based on the fit of Qbar(A,W) from the training set
    V1<- V0<- valid
    V1$A= 1; V0$A=0
    
    QbarAW.valid<- predict(glm.out, newdata=valid, type='response')
    Qbar1W.valid<- predict(glm.out, newdata=V1, type='response')
    Qbar0W.valid<- predict(glm.out, newdata=V0, type='response')
    # again ignore any warnings
    
    # get estimates of the propensity score based on the fit of g(A|W)
    #	from the training set
    # (change these probabilities if not 2-arm trial with balanced allocation)
    
    if( is.null(gAdj)){  #know g(A|W)=0.5
      pscore.valid<- rep(0.5, nrow(valid))
    } else if (gAdj=='U'){	
      pscore.valid<- rep(0.5, nrow(valid))
    } else{
      pscore.valid<- predict(p.out, newdata=valid,  type='response')
    }
    
    # calculate the clever covariate
    H.AW.valid<- valid$A/pscore.valid - (1-valid$A)/(1-pscore.valid)
    
    # update
    if(family!='binomial'){
      QbarAW.valid<- QbarAW.valid + eps*H.AW.valid	
      Qbar1W.valid<- Qbar1W.valid + eps/pscore.valid
      Qbar0W.valid<- Qbar0W.valid - eps/(1-pscore.valid)
    } else {	
      QbarAW.valid<- plogis( qlogis(QbarAW.valid) + eps*H.AW.valid)	
      Qbar1W.valid<- plogis( qlogis(Qbar1W.valid) + eps/pscore.valid)
      Qbar0W.valid<- plogis( qlogis(Qbar0W.valid) - eps/(1-pscore.valid))
    } 
    
    RETURN<- list(pscore=pscore.valid, H.AW=H.AW.valid, QbarAW=QbarAW.valid, 
                  Qbar1W=Qbar1W.valid, Qbar0W=Qbar0W.valid, psi.hat=psi.hat)
  }
  RETURN
}



#------------
# get.Variance.IC: function to get the variance based on the estimated influence
# curve (IC)
# 	input: residuals + estimates of the DY and DW components of the IC 
#		also need independent unit
#	output: variance for PATE or SATE depending on the study design
#-------------

get.Variance.IC<- function(resid, DY, DW, indpt){
  
  # without pair-matching
  if(!PAIRED){
    var.PATE<- var(DY + DW) /n	 # (Eq. 3)
    var.SATE<- var(DY)/n  # (Eq. 4)
    
  } else	{ 
    
    # if pair-matched trial, need to know the independent unit 	
    temp<- unique(indpt)
    n.pairs<- length(temp)
    
    # For PATE, need to acct for covariance of residuals (Section 4)
    # For SATE, can approx IC with DbarY= 1/2 sum_{i in pairs} HAW_i*(Y_i -Qbar_i)
    resid.paired<- DY.paired<- rep(NA, n.pairs)
    
    for(i in 1:n.pairs){		
      resid.paired[i]<- (resid[indpt==temp[i]][1])*(resid[indpt==temp[i]][2]) #eq7
      DY.paired[i]<- 0.5*sum(DY[ indpt== temp[i]])  
    }
    
    # For PATE, (variance as if iid - 2*covariance)/n
    var.PATE<- (var(DY + DW) - 2*mean(resid.paired))/n
    
    var.SATE<- var(DY.paired)/n.pairs
    
  }
  
  data.frame(var.PATE, var.SATE)
}

#----------------------------
# get.Inference.tdist: get inference assuming a Student's t-distribution
#
#	input: truth (true value of psi), psi.hat (point estimate), alpha 
# (signifance level), df (degrees of freedom) var (variance estimate)
# 
# output: variance, tstat, coverage, reject null.
#--------

get.Inference.tdist<- function(truth, psi.hat, alpha, df, var){
  
  # cutoff based on t-dist for testing and CI	
  cutoff <- qt(alpha/2, df=df, lower.tail=F)
  
  # standard error (square root of the variance)
  se<- sqrt(var)
  
  # confidence interval coverage
  cover<- ( (psi.hat - cutoff*se) <= truth & truth <= (psi.hat + cutoff*se) )
  
  # test statistic and pvalue
  tstat <- psi.hat/se
  reject<- 2*pt(abs(tstat), df=df, lower.tail=F) < alpha
  
  data.frame(var, tstat, cover, reject)
}

#--------------------------------------------
# CV.selector: function to select the candidate working model (TMLE)
# that minimizes the estimated variance for Step 1 Initial estimation 
#	and in Step 2 Targeting.
#
# input: data, family for Qbar(A,W), forQ (CV selection for Qbar(A,W) or g(A|W))
#	  QAdj (a priori spec adj variable for Qbar(A,W)), 	
#   gAdj (a priori spec adj variable for g(A|W))
#
#	output: candidate adjustment variable with smallest CV-risk
#-------------------------------------------

CV.selector<- function(data, family, forQ, QAdj, gAdj){
  
  # number of potential TMLEs correspond to num of possible adj variables
  num.tmles <- length(ADJ.SET)
  CV.risk <- rep(NA, num.tmles)
  
  # since the folds are initialized before running any algorithms,
  # we are going to loop through the candidate estimators
  for(k in 1:num.tmles){	
    
    if(forQ){ 
      # adaptive pre-specified approach for Step 1: Initial Estimation
      # (Sec. 3.1 for non-matched trial and Sec. 4.1 for matched trial)
      CV.risk[k]<- get.CV.inference(data=data, QAdj=ADJ.SET[k], family=family, 
                                    gAdj=NULL, get.variance=F)
      
    } else{
      # adaptive pre-specified approach for Step 2: Targeting (based on selected
      # adjustment variable for Qbar(A,W)). (Sec 5.1)
      CV.risk[k]<- get.CV.inference(data=data, QAdj=QAdj, family=family, 
                                    gAdj=ADJ.SET[k], get.variance=F)     
    }
    
  }	
  
  # SELECT THE ADJUSTMENT VARIABLE RESULTING IN THE SMALLEST
  # CV-RISK ESTIMATE 
  ADJ.SET[ which.min(CV.risk)]
}

#-------------------------------
# get.CV.inference: function to get an CV-estimate of the variance
#	also used by the CV-selector to choose the candidate 
# with the smallest CV-risk
#
# input: data, QAdj (adjustment variable for Qbar(A,W) ), family for Qbar(A,W), 
#		gAdj (adjustment variable for g(A|W) )
#		get.variance=T if calc CV-variance; O/W calc risk
#
# output: cross-validated estimates of variance or risk
#--------------

get.CV.inference<- function(data, QAdj, family, gAdj, get.variance){
  
  # this is programmed so that the "Pair" represents the indpt unit
  # 	for a non-matched trial, Pair is a dummy variable for each observation
  
  # We are doing LOO-CV with "Pair" as the indpt unit
  # For other CV-schemes, change the specification of folds
  folds<- data$folds
  index.folds<- unique(folds)
  nFolds<- length(index.folds) 
  
  # Need estimates of DY+DW components of the IC (and possibly residuals)
  #	for obs in the validation set
  DY<- DW<- resid<- rep(NA,n )
  
  # if goal is estimating the risk, record average loss for each fold
  risk<- rep(NA,  nFolds)
  
  # Now run the TMLE algorithm using data in the training set
  # and obtain estimates of the clever covariate H.AW and targeted predictions
  # for observations Qbar*(A,W), Qbar*(1,W), Qbar*(0,W) in the validation set
  
  for(i in 1: nFolds) {
    
    valid <- data[folds==index.folds[i], ]
    train <- data[folds!=index.folds[i],]
    
    # run full TMLE on training set, but evaluate for validation set
    valid.out<- do.TMLE.with.CV(Do.CV=T, train=train, valid=valid, QAdj=QAdj, 
                                family=family,  gAdj=gAdj)	
    
    resid[folds==index.folds[i]] <- valid$Y - valid.out$QbarAW                
    DY[folds==index.folds[i]]<- valid.out$H.AW*(valid$Y - valid.out$QbarAW)   
    DW[folds==index.folds[i]]<- valid.out$Qbar1W - valid.out$Qbar0W- valid.out$psi.hat
    
    if(!get.variance){ # estimate the risk for each fold
      risk[i]<- calc.risk(resid= resid[folds==index.folds[i]], DY=DY[folds==index.folds[i]], 
                          DW=DW[folds==index.folds[i]], indpt=valid$Pair)
    }    
    
  }
  
  # if goal is to get a CV variance estimate
  if(get.variance){
    RETURN<- get.Variance.IC(resid=resid, DY=DY, DW=DW, indpt=data$Pair)
  }else{
    # otherwise average the CV-risk estimates over the folds
    RETURN<- mean(risk)	
  }
  RETURN
}



#----------------------------------
# calc.risk: function to estimate the risk for each fold
# input: estimates of resid, DY, DQ + outcome Y for obs in valid set, 
#	also need indpt unit
#	also relies on POP.EFECT + PAIRED (global variables)
# output: estimated IC for PATE (rho) and SATE
#------------------------------------
calc.risk<- function(resid, DY, DW, indpt){
  
  index.indpt<- unique(indpt)
  
  n.indpt <- length(index.indpt)
  loss<- rep(NA,  n.indpt )
  # for PATE in pair-matching, we can ignore the design
  loss2<-matrix(NA, nrow=n.indpt, ncol=2)
  
  # for each independent unit (obs or pair)	
  for(k in 1: n.indpt) {
    
    these<- indpt==index.indpt[k]
    
    # evaluate the loss functions for obs in valid set
    if(POP.EFFECT){
      
      if(!PAIRED){
        loss[k]<- (DY[these] + DW[these])^2 		#Eq5
      } else{
        if(!DO.RHO){ 
          # calc loss ignoring the pair-matching scheme
          loss2[k,] <- (DY[these] + DW[these])^2  
        }else{
          # alternative loss function for the less conservative variance estimator #Eq8
          loss[k]<- 0.5*sum( (DY[these]+DW[these])^2) - 2*(resid[these][1])*(resid[these][2])
        }
      }
      
    } else{
      # if SATE
      if(PAIRED){
        loss[k]<- ( 0.5*(DY[these][1] + DY[these][2]) )^2  #E9
      } else{
        loss[k]<- (DY[these] )^2 		#Eq6
      }
    }	
    
  }
  if(POP.EFFECT & PAIRED & !DO.RHO){
    loss<- loss2
  }
  
  # estimate the risk for fold v with average loss 
  mean(loss)
}



######################################
# functions to generate the data for either Simulation Study
get.Data<- function(n){  
  
  # generate the binary {-1,1} variable R (only used in Sim2)
  R<- rbinom(n, 1, .5)
  R[R==0]<- -1
  
  # 9 baseline covariates with specified correlation structure	
  s<- 1
  Sigma<- matrix(0.5*s*s, nrow=3, ncol=3)
  diag(Sigma)<- s^2
  
  W1<- mvrnorm(n, rep(0,3), Sigma)
  W2<- mvrnorm(n, rep(0,3), Sigma)
  W3<- cbind(rnorm(n,0,s), rnorm(n,0,s), rnorm(n,0,s))
  
  # generate baseline Z (only used in Sim2)
  # correlated with the sum of the baseline covariates 	
  Z<- R*plogis( W1[,1] + W2[,1] + W3[,1] + .5*rnorm(n,0, s) ) 
  
  # unmeasured factor impacting the outcome	
  UY<- rnorm(n,0,1)
  
  # calculate the counterfactuals
  if(SIM==1){
    Y0<- get.Y.Sim1(W1=W1, W2=W2, W3=W3, UY=UY, A=0)
    Y1<- get.Y.Sim1(W1=W1, W2=W2, W3=W3, UY=UY, A=1)
    
  }else{
    Y0<- get.Y.Sim2(W1=W1, W2=W2, W3=W3, R=R, Z=Z, UY=UY, A=0)
    Y1<- get.Y.Sim2(W1=W1, W2=W2, W3=W3, R=R, Z=Z, UY=UY, A=1)
  }
  
  # add-on a dummy variable and return 
  data.frame(U=rep(1,n), R=R, W=cbind(W1,W2,W3), Z=Z, Y1=Y1, Y0=Y0)
}


#-------------------
# generate the outcome Y as a function of the baseline covariates, exposure 
#	and unmeasured UY

get.Y.Sim1<- function(W1, W2, W3, UY, A){
  0.4*A + 0.25*W1[,1]+ 0.25*W1[,2]+ 0.25*W2[,1]+ 0.25*W2[,2]+ 0.25*UY + 0.25*A*W1[,1] +0.25*A*UY
}

get.Y.Sim2<- function(W1, W2, W3, R, Z, UY, A){
  plogis(0.75*A + 0.5*W1[,2]+ 0.5*W2[,2]+ 0.5*W3[,2] + 1.5*Z+ 0.25*UY +0.75*A*W1[,2] - 0.75*A*W2[,2] + 0.5*A*Z )/7.5 
}


##########################################################
# calculate the population average treatment effect
#  as average diff in counterfactuals over a pop of 900,000
get.PATE<- function(N=900000){
  
  X<- get.Data(n=N)
  
  mean(X$Y1 - X$Y0)
}


###########################################
#----------------------------------------------------
###########################################


library('nbpMatching')
library('MASS')

set.seed(123)

# Specify the Simulation Study 1 or 2
SIM<<- 2

# Specify the target of inference as Population ATE or Sample ATE
POP.EFFECT<<- F

# Specify the study design 
PAIRED <<- T

# If Simulation2, specify the matching set 
SET<- 1

# The code includes two ways to calculate the loss for the PATE
# in a pair-trial: with and without the correction term rho for the
# covariance of the residuals within matched pairs
# (The correction rho is always used when estimating the variance 
# of the TMLE for the PATE in a pair-matched trial)

DO.RHO<<- T



#--------
if(SIM==1){
	n<<- 40
	# specify the matching set (covariates W1...W6)
 	matchOn<<- c('W.1', 'W.2', 'W.3', 'W.4', 'W.5', 'W.6')
 	
	# Specify the potential adjustment set (i.e. the library)
	ADJ.SET<<- c('U', paste('W', 1:9, sep='.') )
	
	# for the conditional mean outcome Qbar(A,W) 
	FAMILY<<- 'gaussian' 
 	
}else{ #SIM2
	n<<-30
		
	# also specify the matching set
	if(SET==1){
		matchOn<<- 'R'
	}else if(SET==2) {
		matchOn<<-  c('R', 'W.2', 'W.5', 'W.8')
	} else{
		matchOn<- NULL	
	}
	
	# Specify the potential adjustment set (i.e. the library)
		ADJ.SET<<- c('U', 'R', paste('W', 1:9, sep='.'),  'Z')
	
	#  for the conditional mean outcome Qbar(A,W) 
	FAMILY<<- 'binomial' 
	
}
n.pairs<<- n/2  



PATE<<- get.PATE()

simulate.Data.and.Run()
