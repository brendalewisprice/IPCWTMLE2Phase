## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----data_read_in, echo=FALSE--------------------------------------------
# Load the data sets available in the package IPCWTMLE
load("/Users/brendaprice/Documents/Dissertation/IPCWTMLE/data/simCYD14_twostage.RData")
load("/Users/brendaprice/Documents/Dissertation/IPCWTMLE/data/simCYD15_twostage.RData")
source("/Users/brendaprice/Documents/Dissertation/IPCWTMLE/R/FunctionsForPackage09242019.R")

## ----data_read_in2, echo=TRUE, eval=FALSE--------------------------------
#  # Load the data sets and needed packages
#  # call libraries
#  library(SuperLearner)
#  library(IPCWTMLE)
#  data(simCYD14_twostage)
#  data(simCYD15_twostage)

## ----data_summary--------------------------------------------------------
dim(simCYD14_twostage)
dim(simCYD15_twostage)

## ----data_summary2-------------------------------------------------------
dat <- rbind(simCYD14_twostage,simCYD15_twostage)
dim(dat)

## ----data_summary3, eval=TRUE--------------------------------------------
dat$A <- ifelse(dat$M13_PRNT_SeroAverage >= median(dat$M13_PRNT_SeroAverage),1,0)

## ----data_summary4, eval=TRUE--------------------------------------------
SL.W <- c("VACC","AGEYRS","MALE","MYS","IDN","THA","VNM","COL","HND","PRI", "MEX","M13_MNv2_SeroAverage" )
fmla = as.formula(paste("VCD ~ A + factor(VACC) + AGEYRS + MALE + MYS + IDN + THA + VNM + COL + HND + PRI + MEX"))


## ----data_summary5, eval=TRUE,comment=NA, message=FALSE, warning=FALSE----
out <- weighted.tmle.separate(w=dat[,SL.W], 
                       a=dat$A, 
                       y=dat$VCD, 
                       delta=dat$delta, 
                       Q.SL.library=c("SL.mean","SL.glm" ), #,"SL.step","SL.glmnet"), 
                       g.SL.library=c("SL.mean","SL.glm"), 
                       wgts=dat$weight,
                       max.wgt = Inf)

## ----outputresults, eval=TRUE--------------------------------------------
paste0("EY1: ",round(out$EY1,3),"; 95% CI: (",round(out$ci1[1],3),",",round(out$ci1[2],3) ,")")
paste0("EY0: ",round(out$EY0,3),"; 95% CI: (",round(out$ci0[1],3),",",round(out$ci0[2],3) ,")")

## ----outputresults2, eval=TRUE-------------------------------------------
paste0("Risk Difference: ",round(out$est.RiskDiff,3),"; 95% CI: (",round(out$ci.RiskDiff[1],3),",",round(out$ci.RiskDiff[2],3) ,")")
paste0("RR: ",round(out$est.RR,3),"; 95% CI: (",round(out$ci.rr[1],3),",",round(out$ci.rr[2],3) ,")")
paste0("OR: ",round(out$est.OR,3),"; 95% CI: (",round(out$ci.OR[1],3),",",round(out$ci.OR[2],3) ,")")


## ----data_summary6, eval=TRUE,comment=NA, message=FALSE, warning=FALSE----
out.fitted <- weighted.tmle.predicted.values(w=dat[,SL.W], 
                       a=dat$A, 
                       y=dat$VCD, 
                       delta=dat$delta, 
                       Q.SL.library=c("SL.mean","SL.glm"), 
                       g.SL.library=c("SL.mean"), 
                       wgts=dat$weight,
                       max.wgt = Inf)

## ----data_summary7, eval=TRUE,comment=NA, message=FALSE, warning=FALSE----
head(out.fitted$fitted.values)

## ----data_summary8, eval=TRUE,comment=NA, message=FALSE, warning=FALSE----
# Make sure osDesign is loaded
  if (!require(osDesign)) install.packages('osDesign')
  library(osDesign)
# Phase I counts for use in the tps function
nn0<-c(sum(dat$weight[dat$VCD==0 &dat$VACC==0]),sum(dat$weight[dat$VCD==0 &dat$VACC==1]))
nn1<-c(sum(dat$weight[dat$VCD==1 &dat$VACC==0]),sum(dat$weight[dat$VCD==1 &dat$VACC==1]))

# osDesign requires the stratification variable to be coded 1,2,..K for K strata
dat$S <- dat$VACC + 1
    
# using osDesign package, calculate the BH estimator ("method='ML'")
# using only the phase 2 sample of data, but includeing the phase 1 counts as nn0 and nn1
dat.2 <- dat[dat$delta==1,]
tps.fit <- tps(fmla, data=dat.2, nn0=nn0, nn1=nn1, group=dat.2$VACC, method="ML")
  
## Run the package function to output the corresponding estimates for the BH approach
out.BH <- BH.Estimates(fmla=fmla,dat=dat,tps.fit=tps.fit)
## Print the corresponding estimates for the BH approach
out.BH

