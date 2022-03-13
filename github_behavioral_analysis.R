library(dplyr)
library(jtools)

#setwd()

#loading data
data <- read.csv("alldata.csv",header = T, sep = ";")
head(data)
names(data) <- c('subject', 'subject_ID', 'session', 'outcome1', 'outcome2', 'outcome3', 'outcome4', 'prob1','prob2', 'choice', 'rt', 'Onsettime', 'gender', "order")
attach(data) 

subjList <- unique(data[, "subject_ID"]) 
subjList <- subjList[c(-14,-16)]  #excluding two participants with interhemispheric asymmetry in motor threshold
data1 <- data[(subject_ID %in% subjList),]

#loading personality scales and stimulation discomfort data
scales <- read.csv('scales.csv')
ss_key <- c(1,1,2,2,2,2,2,2,1,2,1,2,2,2,1,2)

myf <- function(subj) { 
  s = length(which(scales[which(scales$Subject_ID==subj),33:48]==ss_key))
  return(s)}

scales$SS = mapply(myf, scales$Subject_ID)
data1$SS = mapply(myf, data1$subject_ID)

scales$Drive = scales$B3+scales$B9+scales$B12+scales$B21
scales$Fun = scales$B5+scales$B10+scales$B15+scales$B20
scales$Reward = scales$B4+scales$B7+scales$B14+scales$B18+scales$B23
scales$BIS = scales$B2+scales$B8+scales$B13+scales$B16+scales$B19+scales$B22+scales$B24

myf1 <- function(subj) {
  s1 = scales[which(scales$Subject_ID==subj),'Drive']
  s2 = scales[which(scales$Subject_ID==subj),'Fun']
  s3 = scales[which(scales$Subject_ID==subj),'Reward']
  s4 = scales[which(scales$Subject_ID==subj),'BIS']
  s5 = scales[which(scales$Subject_ID==subj),'Discomf1']
  s6 = scales[which(scales$Subject_ID==subj),'Discomf2']
  s7 = scales[which(scales$Subject_ID==subj),'Discomf3']
  return(cbind(s1,s2,s3,s4,s5,s6,s7))
}

suppl= mapply(myf1, data1$subject_ID)
rownames(suppl)=c('Drive', 'Fun', 'Reward','BIS', 'Discomf1', 'Discomf2', 'Discomf3')
data1 <- cbind(data1, t(suppl))
data1 <- mutate(data1, trial = rep(1:96, 3*length(subjList)))
data1 <- mutate(data1, discomf = ifelse( order==1, Discomf1, ifelse(order==2, Discomf2, Discomf3)))


####### Regression analysis #########
library(lmerTest)
library(lme4)

data1 <- mutate(data1, EV_A = outcome1*prob1 + outcome2*prob2)
data1 <- mutate(data1, EV_B = outcome3*prob1 + outcome4*prob2)
data1 <- mutate(data1, sd_A = sqrt(((outcome1-EV_A)^2)*prob1 + ((outcome2-EV_A)^2)*prob2))
data1 <- mutate(data1, sd_B = sqrt(((outcome3-EV_B)^2)*prob1 + ((outcome4-EV_B)^2)*prob2))

data1 <- mutate(data1, diff_m = ifelse( abs(outcome1)>abs(outcome3), EV_A-EV_B, EV_B-EV_A))
data1 <- mutate(data1, diff_sd = ifelse(abs(outcome1)>abs(outcome3), sd_A-sd_B, sd_B-sd_A))

data1 <- mutate(data1, risky_sd = ifelse( (sd_A>sd_B & choice==1) | ( sd_A<sd_B & choice==2), 1, 0))
data1 <- mutate(data1, choice_rat = ifelse(EV_A>EV_B, 1, 2))
data1 <- mutate(data1, choice_is_rat = ifelse(choice==choice_rat, 1, 0))

mo1<- lme4::glmer(risky_sd ~ as.factor(session)*discomf + trial + diff_sd*diff_m + (1|subject_ID), family = 'binomial', data=data1[which(outcome1>0),]) 
mo2<- lme4::glmer(risky_sd ~ as.factor(session)*discomf + trial + diff_sd*diff_m + (1|subject_ID), family = 'binomial', data=data1[which(outcome1<0),]) 
mo3<- lme4::glmer(choice_is_rat ~ as.factor(session)*discomf + trial + diff_sd*diff_m + (1|subject_ID), family = 'binomial', data=data1[which(outcome1>0),]) 
mo4<- lme4::glmer(choice_is_rat ~ as.factor(session)*discomf + trial + diff_sd*diff_m + (1|subject_ID), family = 'binomial', data=data1[which(outcome1<0),]) 


### recalculation of gradient to check convergence issues (abs(grad)<0.002)
library(numDeriv)
model = mo1 #substitute with model name from above
dd <- update(model,devFunOnly=TRUE)
pars <- unlist(getME(model,c("theta","fixef")))
grad2 <- grad(dd,pars)
hess2 <- hessian(dd,pars)
sc_grad2 <- solve(hess2,grad2)
max(pmin(abs(sc_grad2),abs(grad2)))
#####


library(stargazer)
stargazer(mo1, mo2, mo3, mo4, type="html",
          out="model_output.doc", intercept.bottom = F, star.cutoffs = c(0.05, 0.01, 0.001), intercept.top = T, report = "vc*s"  , digits=3)
