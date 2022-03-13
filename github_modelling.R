library(hBayesDM)
library(rstan)
library(doParallel)

#setwd()
data <- read.csv("alldata.csv",header = T, sep = ";")
head(data)
names(data) <- c('subject', 'subject_ID', 'session', 'outcome1', 'outcome2', 'outcome3', 'outcome4', 'prob1','prob2', 'choice', 'rt', 'Onsettime', 'gender')
attach(data) 
options(digits=5)
options(scipen=999)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


#modelPath <- "/home/kpanidi/HBM_g_all.stan"
modelPath <- "stancode.stan"
m_gamma = rstan::stan_model(modelPath)


##############  preparing dataset to be submitted to rstan  ############## 

raw_data <- c()
raw_data <- data [ which (data$outcome1>0) ,] 

#fit_prw <- function(raw_data){
raw_data1 <- c()
raw_data1 <- raw_data [ which (raw_data$outcome1>0 ),] 
raw_data1 <-cbind(raw_data1, session1=0, session2=0)
raw_data1[which(raw_data1$session==0),]$session1=0
raw_data1[which(raw_data1$session==0),]$session2=0
raw_data1[which(raw_data1$session==1),]$session1=1
raw_data1[which(raw_data1$session==1),]$session2=0
raw_data1[which(raw_data1$session==2),]$session1=0
raw_data1[which(raw_data1$session==2),]$session2=1


subjList <- unique(raw_data1[, "subject_ID"]) 
subjList <- subjList[c(-14,-16)] #excluding two participants with interhemispheric asymmetry in motor threshold
numSubjs=length(subjList)
maxTrials=48+48+48
Tsubj=as.vector(rep(maxTrials, numSubjs))

outcome1 <- array(0, c(numSubjs, maxTrials))
outcome2 <- array(0, c(numSubjs, maxTrials))
outcome3 <- array(0, c(numSubjs, maxTrials))
outcome4 <- array(0, c(numSubjs, maxTrials))
choice0 <- array(0, c(numSubjs, maxTrials))
probs <- array(0, c(numSubjs, maxTrials))
session <- array(0, c(numSubjs, maxTrials))

for (i in 1:numSubjs) {
  
  curSubj <- subjList[i]
  useTrials <- Tsubj[i]
  tmp <- subset(raw_data1, raw_data1$subject_ID == curSubj)
  outcome1[i, 1:useTrials] <- tmp[1:useTrials, "outcome1"]
  outcome2[i, 1:useTrials] <- tmp[1:useTrials, "outcome2"]
  outcome3[i, 1:useTrials] <- tmp[1:useTrials, "outcome3"]
  outcome4[i, 1:useTrials] <- tmp[1:useTrials, "outcome4"]
  choice0[i, 1:useTrials] <- tmp[1:useTrials, "choice"]
  choice<- (choice0-2)*(-1)
  probs[i, 1:useTrials] <- tmp[1:useTrials, "prob1"]
  session[i, 1:useTrials] <- tmp[1:useTrials, "session"]
}


dataList <- list(N = numSubjs, T = maxTrials, Tsubj = Tsubj, 
                 outcome1 = outcome1, outcome2 = outcome2, outcome3 = outcome3, outcome4 = outcome4, choice = choice, probs=probs, session=session)


#this function can be then submitted for parallel computations for a desired number of chains
fun_sample <- function(i){
  s <- rstan::sampling(m_gamma, seed=12345, verbose=TRUE, show_messages=TRUE, data = dataList, pars = NA, warmup = 1000, init = 'random', iter = 5000, chains = 1, chain_id=i, thin = 1, control = list(adapt_delta = 0.999, max_treedepth = 11, stepsize = 1), refresh=10)
  return(s)
}



