## R package functions ##
## Author: Brenda Price  ##
## 09/23/2019 ##

## IPCW-TMLE ##

expit <- function(x){exp(x)/(1 + exp(x))}

MSEfromProb <- function(Y,pred){
  MSE <- mean((Y-pred)^2)
  return(MSE)
}
RMSEfromProb <-function(Y,pred){
  MSE <- mean((Y-pred)^2)
  RMSE <-sqrt(MSE)
  return(RMSE)
}
AUCfromProb <- function(Y,pred){
  if (!require(pROC)) install.packages('pROC')
  library(pROC)
  auc <- auc(Y, pred)[1]
  return(auc)
}


ProbFrom_tps <- function(tps_object,dset,model.dset,x){
  # need counts from Phase 2 data used for model fitting
  n <- table(model.dset$Y,model.dset$S)
  # X is model matrix for formula used in in tps.
  # dset needs an "S" variable supplied.
  expit <-function(l) {exp(l)/(1+exp(l))}
  XBeta <- x %*% (as.matrix(tps_object$coef))
  fitted.values <- expit(XBeta) ## output predicted values
  return(fitted.values)
}

#' @title IPCWTMLE2Phase
#'
#' IPCW-TMLE to output the TMLEs of predicted values
#'
#' @description IPCW-TMLE outputting TMLEs of predicted values. This function only performs the TMLE step for targeting the predicted values. And does not output any additional estimates.
#'
#' @param w The covariates
#' @param a The treatment or exposure binary variable.
#' @param y The outcome
#' @param Q.SL.library SuperLearner library of algorithms for Q estimation
#' @param g.SL.library SuperLearner library of algorithms for g estimation
#' @param wgts observation level weights; should be 0 for delta=0
#' @return fitted.values The predicted values for each observation from the TMLE
#'
#' @export
#' @examples
#' weighted.tmle.predicted.values

## TMLE of predicted values only
weighted.tmle.predicted.values <- function(w,a,y,delta,Q.SL.library,g.SL.library,wgts,max.wgt=Inf){
  if (!require(SuperLearner)) install.packages('SuperLearner')
  library(SuperLearner)
  n <- nrow(w)
  WA <- data.frame(w,A=a)
  wgts = n*wgts/sum(wgts)
  temp1 <- SuperLearner(y,data.frame(w,A=a),SL.library=Q.SL.library,newX=rbind(WA),
                        obsWeights=wgts,family=binomial,cvControl=list(V=10))
  Qbar.ests <- temp1$SL.predict[,1]
  # obtain an estimate of g_n
  temp <- SuperLearner(a,w,SL.library=g.SL.library,newX=w,obsWeights=wgts,
                       family=binomial,cvControl=list(V=10))
  g.ests <- temp$SL.predict[,1]
  rm(temp)
  g.ests[a==0] <- 1-g.ests[a==0]
  a.ind <- 2*a-1
  offset <- qlogis(pmin(pmax(Qbar.ests[1:n],0.0005),0.9995))
  g.and.wgts = pmin(wgts/g.ests,max.wgt)
  eps = coef(glm(y ~ -1 + offset(offset) + a.ind,weights=g.and.wgts,family=binomial))
  Q1 <- plogis(qlogis(Qbar.ests[1:n]) + eps)
  Q0 <- plogis(qlogis(Qbar.ests[1:n]) - eps)
  QA <- ifelse(a==1,Q1,Q0)
  QA <- ifelse(delta==1,QA,NA)
  fitted.values <- QA
  return(list(fitted.values = fitted.values))}

#' @title IPCWTMLE2Phase
#' Main IPCW-TMLE function for two-phase studies
#'
#' @description IPCW-TMLE function to estimate E[E[Y|A=1]] and E[E[Y|A=0]] separately, and to also estimate the marginalized risk difference, relative risk, and odds ratio.  This function assumes that the full phase one data set is input with weights equal to 0 for those observations not selected for the second phase.  An indicator variable, delta, should equal 1 for all phase-two data and 0 for all other observations (non-phase-two).
#'
#' @param w The covariates
#' @param a The treatment or exposure binary variable.
#' @param y The outcome
#' @param delta Indicator for inclusion into the second phase
#' @param Q.SL.library SuperLearner library of algorithms for Q estimation
#' @param g.SL.library SuperLearner library of algorithms for g estimation
#' @param wgts Observation level weights; should be 0 for delta=0
#' @return Estimates of marginalized risk difference, RR, and OR
#'
#' @export
#' @examples
#' weighted.tmle.separate

weighted.tmle.separate <- function(w,a,y,delta,Q.SL.library,g.SL.library,wgts,max.wgt=Inf){
  ## assumes full Phase 1 dataset used, but weights>0 for only those in Phase II
  if (!require(SuperLearner)) install.packages('SuperLearner')
  library(SuperLearner)
  n <- nrow(w) # total counts of Phase I data
  WA <- WA0 <- WA1 <-data.frame(w,A=a)
  WA0$A <- 0 # counterfactual for control
  WA1$A <- 1 # counterfacutal for case
  wgts = n*wgts/sum(wgts) # standardize weights for stability within these calculations
  temp1 <- SuperLearner(y,data.frame(w,A=a),SL.library=Q.SL.library,newX=rbind(WA,WA0,WA1),
                        obsWeights=wgts,family=binomial,cvControl=list(V=10))
  Qbar.ests <- temp1$SL.predict[,1]
  rm(temp1)
  temp <- SuperLearner(a,w,SL.library=g.SL.library,newX=w,obsWeights=wgts,
                       family=binomial,cvControl=list(V=10))
  g.ests <- temp$SL.predict[,1]
  rm(temp)
  g.ests[a==0] <- 1-g.ests[a==0]
  a.ind <- 2*a-1
  a1 <- ifelse(a.ind==1,1,0)
  a0 <- ifelse(a.ind==-1,1,0)
  offset <- qlogis(pmin(pmax(Qbar.ests[1:n],0.0005),0.9995))
  g.and.wgts = pmin(wgts/g.ests,max.wgt)
  # eps = coef(glm(y ~ -1 + offset(offset) + a.ind,weights=g.and.wgts,family=binomial))
  eps1 = coef(glm(y[a1==1] ~ -1 + offset(offset[a1==1]) + a1[a1==1],weights=g.and.wgts[a1==1],family=binomial))
  eps0 = coef(glm(y[a0==1] ~ -1 + offset(offset[a0==1]) + a0[a0==1],weights=g.and.wgts[a0==1],family=binomial))

  Q0 <- plogis(qlogis(Qbar.ests[(n+1):(2*n)]) + eps0) # Qbar.ests[(n+1):(2*n)]=WA0 predicted values
  Q1 <- plogis(qlogis(Qbar.ests[(2*n+1):(3*n)]) + eps1) # Qbar.ests[(2*n+1):(3*n)]=WA1 predicted values
  # estimates
  est.RR <- mean((Q1[a1==1])*wgts[a1==1])/mean((Q0[a0==1])*wgts[a0==1])
  est.RiskDiff <- mean((Q1[a1==1])*wgts[a1==1]) - mean((Q0[a0==1])*wgts[a0==1])
   est.logOR <- mean(wgts[a1==1]*log(Q1[a1==1]/(1-Q1[a1==1]))) - mean(wgts[a1==0]*log(Q0[a0==1]/(1-Q0[a0==1])))# weight on log scale
   est.OR <- exp(est.logOR)
  EY1 <- mean(Q1[a1==1]*wgts[a1==1])
  EY0 <- mean(Q0[a0==1]*wgts[a0==1])

  # confidence intervals
  ic1 = (a1*g.and.wgts*(y- Q1) + a1*wgts*(Q1-EY1)) #[a1==1] # equals only wgts*(Q1-EY1) for placebo subjects
  ic0 = (a0*g.and.wgts*(y- Q0) + a0*wgts*(Q0-EY0))#[a0==1] # equals only  wgts*(Q0-EY0) for vaccine subjects
  se1 <- sqrt(sum(ic1[wgts>0 & a1==1]^2)/sum(wgts[a1==1])^2) # sd(ic1)/sqrt(length(ic1)) 05022019
  se0 <- sqrt(sum(ic0[wgts>0 & a0==1]^2)/sum(wgts[a0==1])^2) # sum of squared IC / n
  ci1 <- c(EY1-qnorm(0.975)*se1,EY1+qnorm(0.975)*se1)
  ci0 <- c(EY0-qnorm(0.975)*se0,EY0+qnorm(0.975)*se0)
  ic.log.rr <- (1/EY1)*(a1*g.and.wgts*(y- Q1) + wgts*(Q1)) -
               (1/EY0)*((a0)*g.and.wgts*(y- Q0) + wgts*(Q0)) # all non-phaseII should be weighted to 0
  var.log.rr.sum <- sum( (ic.log.rr-mean(ic.log.rr))^2)#var(ic.log.rr)
  se.log.rr <- sqrt(var.log.rr.sum)/sqrt(n^2)
  se.rr <- sqrt(exp(2*est.RR)*var.log.rr.sum)/sqrt(n^2) #5022019
  ci.rr <- exp(c(log(est.RR)-qnorm(0.975)*se.log.rr,
              log(est.RR)+qnorm(0.975)*se.log.rr))

  ic.RiskDiff<- a.ind*g.and.wgts*(y-(a*Q1 + (1-a)*Q0)) + wgts*((Q1-Q0) - est.RiskDiff)
  ci.RiskDiff <- c(est.RiskDiff-qnorm(0.975)*sd(ic.RiskDiff)/sqrt(n),
                   est.RiskDiff+qnorm(0.975)*sd(ic.RiskDiff)/sqrt(n))
  QA <- ifelse(a==1,Q1,Q0)
  ic.logOR <- 1/(EY1*(1-EY1)) *(a*g.and.wgts*(y-QA) +wgts*(Q1 -EY1)) -
              1/(EY0*(1-EY0)) *((1-a)*g.and.wgts*(y-QA) +wgts*(Q0 -EY0))
  ci.logOR <- (c((log(est.OR))-qnorm(0.975)*sqrt(var(ic.logOR)/n),
                 (log(est.OR))+qnorm(0.975)*sqrt(var(ic.logOR)/n)))
  ci.OR <- exp(ci.logOR)
  return(list(est.RR=est.RR,ci.rr=ci.rr,EY1=EY1,EY0=EY0,Q0 = Q0, Q1=Q1,
              ci1=ci1, ci0=ci0,se1=se1,se0=se0,se.rr=se.rr,
              se.log.rr=se.log.rr,
              est.RiskDiff=est.RiskDiff,
              ci.RiskDiff=ci.RiskDiff,
              est.OR=est.OR,
              ci.OR = ci.OR))}

logit.var.p.all <- function(X,beta,sigmaB){
  gradient <- matrix(NA, nrow=nrow(X), ncol= length(beta))
  for(j in 1:nrow(X)){
    XB <- X[j,] %*% beta
    eXB <- exp(XB)
    num.list <-list()
    for(k in 1:length(beta)){
      num.list[[k]] <- (X[j,k])*eXB
    }
    num <- unlist(num.list)
    denom <- (1+eXB)^2
    gradient[j,] <- num/denom
  }
  # end j
  var <- (gradient) %*% sigmaB %*% t(gradient)
  return(var)
}

var.p_a <- function(var,weights=PhII$weights[PhII$A==1]){
  grad <- weights # weights are already standardized /sum(weights)
  var.p_a <-(grad %*% var %*% grad)
  return(var.p_a)
}

rr.var <- function(var1,var0,EY1,EY0,cov=0){
  coef <-(EY1/EY0)^2
  rr.var <- coef * (var1/(EY1^2) + var0/(EY0^2) - 2*cov/(EY1*EY0))
  return(rr.var)
}

MCSE<- function(est){
  M <- length(est)
  I.1 <- ((1/M)*sum(est))^2
  I.2 <- (1/M)* sum(est^2)
  MCSE <- sqrt(I.2 - I.1)
  return(MCSE)
}

# @title IPCWTMLE2Phase

#' Estimates for Breslow-Holubkov approach from the "osDesign" package
#'
#' Provides estimates using the BH approach for two stage designs for the marginal means, marginal risk difference, and marginal relative risk.  It requires as input the output object from the "tps" function from the "osDesign" package
#'
#' @param fmla The formula used in the "tps" call from "osDesign"
#' @param dat The phase 1 dataset coded with delta=1 for phase 2 sampled observations
#' @param tps.fit The fit output by the tps function in "osDesign"
#' @return Estimates of marginalized risk difference, relative risk and odds ratios
#'
#' @export
#' @examples
#' See package vignette
#' BH.Estimates

BH.Estimates <- function(fmla=fmla,dat=dat,tps.fit=tps.fit){
  x2 <- model.matrix(fmla[-2], dat)
  dat$Y <- dat[,"VCD"]
  BHprob.train <- ProbFrom_tps(tps_object=tps.fit, dset = dat, model.dset=dat,x2)
  ## calculate the EY1 and EY0 estimates; use standardized weights
  EY1 <- sum(BHprob.train[dat$A==1 & dat$delta==1]*dat$weight[dat$A==1 & dat$delta==1])/sum(dat$weight[dat$A==1 & dat$delta==1])
  EY0 <- sum(BHprob.train[dat$A==0 & dat$delta==1]*dat$weight[dat$A==0 & dat$delta==1])/sum(dat$weight[dat$A==0 & dat$delta==1])

  ## calculate the standard errors
  var1 <- logit.var.p.all(X=x2[dat$A==1 & dat$delta==1,],beta=tps.fit$coef,sigmaB=tps.fit$covm) # variance of A==1 only
  var0 <- logit.var.p.all(X=x2[dat$A==0 & dat$delta==1,],beta=tps.fit$coef,sigmaB=tps.fit$covm) # variance of A==0 only
  var.p_1 <- var.p_a(var1,weights=dat$weight[dat$A==1 & dat$delta==1])
  var.p_0 <- var.p_a(var0,weights=dat$weight[dat$A==0 & dat$delta==1])
  se1 <-sqrt(var.p_1/sum(dat$weight[dat$A==1  & dat$delta==1])^2 )
  se0 <-sqrt(var.p_0/sum(dat$weight[dat$A==0  & dat$delta==1])^2 )

  # calculate the confidence intervals
  ci1 <- c(EY1-qnorm(0.975)*se1,
           EY1+qnorm(0.975)*se1)
  ci0 <- c(EY0-qnorm(0.975)*se0,
           EY0+qnorm(0.975)*se0)

  ## Relative Risk
  RR <- rr <- EY1/EY0

  ## Relative risk CI shoulld be done on the log scale
 # var.log.p_1 <- var.p_1/((EY1)^2)
 # var.log.p_0 <- var.p_0/((EY0)^2 )
  var.log.p_1 <- se1^2/((EY1)^2)
  var.log.p_0 <- se0^2/((EY0)^2 )
 # log.rr.variance <- (var.log.p_1 + var.log.p_0)/nrow(dat)^2
  log.rr.variance <- (var.log.p_1) + (var.log.p_0)
  #log.rr.se <- sqrt(log.rr.variance)
  log.rr.se <- sqrt(log.rr.variance)
  ci.RR <- exp(c(log(rr)-qnorm(0.975)*log.rr.se,
                 log(rr)+qnorm(0.975)*log.rr.se))

  ## Risk Difference
  RiskDiff <- EY1-EY0
  var.RD <- se1^2 + se0^2
  ci.RiskDiff <- c(RiskDiff-qnorm(0.975)*sqrt(var.RD),
                   RiskDiff+qnorm(0.975)*sqrt(var.RD))


  ## Odds Ratio
  OR <- (EY1/(1-EY1))/(EY0/(1-EY0))
  LOR <- log(OR)
  ## OR on log scale
  n1 <- sum(dat$weight[dat$A==1 & dat$delta==1])
  n0 <- sum(dat$weight[dat$A==0 & dat$delta==1])
  # standard delta method for variance of log odds
  log.or.variance <- 1/(n1*EY1) + 1/(n1*(1-EY1)) + 1/(n0*EY0) +  1/(n0*(1-EY0))
  log.or.se <- sqrt(log.or.variance)
  ci.OR <- exp(c(log(rr)-qnorm(0.975)*log.or.se,
                 log(rr)+qnorm(0.975)*log.or.se))


  return(list(EY1=EY1,EY0=EY0,ci1=ci1,ci0=ci0,
              RiskDiff=RiskDiff,
              ci.RiskDiff=ci.RiskDiff,
              RR=RR,ci.RR=ci.RR,
              OR=OR,
              ci.OR=ci.OR))


}

## checking the push

