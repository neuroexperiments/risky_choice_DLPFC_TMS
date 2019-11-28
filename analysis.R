library(plyr)
library(dplyr)
library(maxLik)
library(numDeriv)
library(dfoptim)
library(MASS)
library(ggplot2)
library(ggpubr)
library("BSDA") 
library(gridExtra)
library(grid)
library(lattice)
library(lemon)
library(ggplotify)
library(gtable)
library(data.table)
library(dplyr)
library(formattable)
library(tidyr)


#for csv created in Windows use this:
data <- read.csv(file.choose(), header = T, sep = ";")

#for csv created on Mac (from "Numbers" application) use this:
data <- read.csv(file.choose(), header = T, sep = ",")

names(data) <- c('subject', 'subject_ID', 'session', 'outcome1', 'outcome2', 'outcome3', 'outcome4', 'prob1','prob2', 'choice', 'rt', 'Onsettime', 'gender', 'session_order')

attach(data) 
options(digits=5)
options(scipen=999)


########################################################
### Estimation procedures for Kahneman-Tversky specification
########################################################

N=(1:28)

coefs=c()

# j is for subjects
# k is for sessions

# session = 0 sham right DLPFC
# session = 1 right DLPFC
# session = 2 left DLPFC

gains = 1  #set gains = -1 for loss domain

for (j in N) {
  for (k in c(0,1,2)) {
    
    data1 <- data [ which (subject == j & gains*outcome1>0 & session == k),] 
   
    log_lik <- function(pars) {
      r <- pars[1]
      gamma <- pars[2]
      mu <- pars[3]
    
      p <- data1$prob1
      outcome1 <- data1$outcome1
      outcome2 <- data1$outcome2
      outcome3 <- data1$outcome3
      outcome4 <- data1$outcome4
      
      EU <- c()
      
      for(i in 1:nrow(data1)) {
        
        if (gains==1){
        EU_a_i <- ((p[i]^gamma)/(((p[i])^gamma+(1-p[i])^gamma)^(1/gamma))) * ((( outcome1[i])^(1-r))/(1-r))+ 
          (1 - ((p[i]^gamma)/(((p[i])^gamma+(1-p[i])^gamma)^(1/gamma)))) *(((  outcome2[i])^(1-r))/(1-r))
        
        EU_b_i <- ((p[i]^gamma)/(((p[i])^gamma+(1-p[i])^gamma)^(1/gamma))) *(((  outcome3[i])^(1-r))/(1-r))+ 
          (1 - ((p[i]^gamma)/(((p[i])^gamma+(1-p[i])^gamma)^(1/gamma)))) * (((  outcome4[i])^(1-r))/(1-r))
        
        EU <- c(EU, 1/(1 + (EU_b_i/EU_a_i)^(1/mu))) }
      
       else {
        
        
        EU_a_i <- (1-((1-p[i])^gamma)/((((1-p[i]))^gamma+(1-(1-p[i]))^gamma)^(1/gamma))) * (-(((- outcome1[i])^(1-r))/(1-r))) + 
          ((((1-p[i])^gamma)/((((1-p[i]))^gamma+(1-(1-p[i]))^gamma)^(1/gamma)))) *(-(((- outcome2[i])^(1-r))/(1-r))) 
        
        EU_b_i <- (1-((1-p[i])^gamma)/((((1-p[i]))^gamma+(1-(1-p[i]))^gamma)^(1/gamma))) * (-(((- outcome3[i])^(1-r))/(1-r))) + 
          ((((1-p[i])^gamma)/((((1-p[i]))^gamma+(1-(1-p[i]))^gamma)^(1/gamma)))) * (-(((- outcome4[i])^(1-r))/(1-r)))
        
        
        EU <- c(EU, 1/(1 + (EU_a_i/EU_b_i)^(1/mu))) }
    }
      
      
      
      df <- data.frame( EU = EU, y = data1$choice)
      
      df$optim <- 0
      
      for(i in 1:nrow(data1)) {
        if(df$y[i] == 1) {
          df$optim[i] <- log(EU[i])
        } else {
          df$optim[i] <- log(1 - EU[i])
        }
      }
      
      return(-sum(df$optim))
     
    }
    
    theta = c(0.0001, 0.0001, 0.0001)
    
    par <- nmk(theta, f = log_lik)$par
    
    par <- cbind(j, k, t(par))   
    coefs=rbind(coefs, par)
  }
}

colnames(coefs) <- c("subject","session", "RA", "gamma", "mu")



#removing outliers 

s0 = IQR(coefs[which(coefs[,2]==0),3])  #computes interquantile range
s1 = IQR(coefs[which(coefs[,2]==1),3])
s2 = IQR(coefs[which(coefs[,2]==2),3])

s0_25=quantile(coefs[which(coefs[,2]==0),3],1/4)  #computes 25st quantile
s0_75=quantile(coefs[which(coefs[,2]==0),3],3/4)

s1_25=quantile(coefs[which(coefs[,2]==1),3],1/4)  #computes 25st quantile
s1_75=quantile(coefs[which(coefs[,2]==1),3],3/4)

s2_25=quantile(coefs[which(coefs[,2]==2),3],1/4)  #computes 25st quantile
s2_75=quantile(coefs[which(coefs[,2]==2),3],3/4)

RA_0_lower_fence = s0_25-3*s0
RA_0_upper_fence = s0_75+3*s0

RA_1_lower_fence = s1_25-3*s1
RA_1_upper_fence = s1_75+3*s1

RA_2_lower_fence = s2_25-3*s2
RA_2_upper_fence = s2_75+3*s2

s0 = IQR(coefs[which(coefs[,2]==0),4])  #computes interquantile range
s1 = IQR(coefs[which(coefs[,2]==1),4])
s2 = IQR(coefs[which(coefs[,2]==2),4])

s0_25=quantile(coefs[which(coefs[,2]==0),4],1/4)  #computes 25st quantile
s0_75=quantile(coefs[which(coefs[,2]==0),4],3/4)

s1_25=quantile(coefs[which(coefs[,2]==1),4],1/4)  #computes 25st quantile
s1_75=quantile(coefs[which(coefs[,2]==1),4],3/4)

s2_25=quantile(coefs[which(coefs[,2]==2),4],1/4)  #computes 25st quantile
s2_75=quantile(coefs[which(coefs[,2]==2),4],3/4)

gamma_0_lower_fence = s0_25-3*s0
gamma_0_upper_fence = s0_75+3*s0

gamma_1_lower_fence = s1_25-3*s1
gamma_1_upper_fence = s1_75+3*s1

gamma_2_lower_fence = s2_25-3*s2
gamma_2_upper_fence = s2_75+3*s2

s0_25=quantile(coefs[which(coefs[,2]==0),5],1/4)  #computes 25st quantile
s0_75=quantile(coefs[which(coefs[,2]==0),5],3/4)

s1_25=quantile(coefs[which(coefs[,2]==1),5],1/4)  #computes 25st quantile
s1_75=quantile(coefs[which(coefs[,2]==1),5],3/4)

s2_25=quantile(coefs[which(coefs[,2]==2),5],1/4)  #computes 25st quantile
s2_75=quantile(coefs[which(coefs[,2]==2),5],3/4)

mu_0_lower_fence = s0_25-3*s0
mu_0_upper_fence = s0_75+3*s0

mu_1_lower_fence = s1_25-3*s1
mu_1_upper_fence = s1_75+3*s1

mu_2_lower_fence = s2_25-3*s2
mu_2_upper_fence = s2_75+3*s2


Fences_matrix = matrix( c(RA_0_lower_fence, RA_1_lower_fence, RA_2_lower_fence, RA_0_upper_fence, RA_1_upper_fence, RA_2_upper_fence, gamma_0_lower_fence, gamma_1_lower_fence, gamma_2_lower_fence, gamma_0_upper_fence, gamma_1_upper_fence, gamma_2_upper_fence, mu_0_lower_fence, mu_1_lower_fence, mu_2_lower_fence, mu_0_upper_fence, mu_1_upper_fence, mu_2_upper_fence), ncol=6, nrow=3)

K = c()

for (i in N){
  coefs1 = coefs[which(coefs[,1]==i),]
  if (  ( all(coefs1[,3] > Fences_matrix[,1]) ==TRUE) & ( all(coefs1[,3] < Fences_matrix[,2]) ==TRUE) & ( all(coefs1[,4] > Fences_matrix[,3]) ==TRUE) & ( all(coefs1[,4] < Fences_matrix[,4]) ==TRUE)  & ( all(coefs1[,5] > Fences_matrix[,5]) ==TRUE) & ( all(coefs1[,5] < Fences_matrix[,6]) ==TRUE) ) {K = c(K, i)} 
}


#computing difference in coefficients
coefs_diff_10=c()
coefs_diff_20=c()
coefs_diff_21=c()

for (j in K){
  
  coefs_diff_10 = rbind(coefs_diff_10, coefs[which(coefs[,1]==j & coefs[,2]==1),] - coefs[which(coefs[,1]==j & coefs[,2]==0),] )
  coefs_diff_20 = rbind(coefs_diff_20, coefs[which(coefs[,1]==j & coefs[,2]==2),] - coefs[which(coefs[,1]==j & coefs[,2]==0),] )
  coefs_diff_21 = rbind(coefs_diff_21, coefs[which(coefs[,1]==j & coefs[,2]==2),] - coefs[which(coefs[,1]==j & coefs[,2]==1),] )
  
}

coefs <- coefs[which(is.element(coefs[,1], K) ), ]  #this takes into account removing outliers



########################################################
#creating plots for Kahneman-Tversky specification
########################################################
#risk aversion

df <- data.frame()
a = coefs[which(is.element(coefs[,1], K) & (coefs[,2]==0) ), 3] 
b = coefs[which(is.element(coefs[,1], K) & (coefs[,2]==1) ), 3]
c = coefs[which(is.element(coefs[,1], K) & (coefs[,2]==2) ), 3]

df1 <- data.frame()
df2 <- data.frame()
df3 <- data.frame()

df1 <- data.frame(data=a, session="sham")
df2 <- data.frame(data=b, session="right DLPFC")
df3 <- data.frame(data=c, session="left DLPFC")
df <- rbind(df1, df2, df3)  

stat.test <- compare_means(data~session, method="wilcox.test", data= df, paired=TRUE, p.adjust.method="bonferroni")

myplot1 <- ggboxplot(df,  x="session", y="data", fill = "session", palette = c("#00AFBB", "#E7B800", "#FC4E07"))
myplot1 <- myplot1+theme_light()
my_comparisons <- list( c("sham", "right DLPFC"), c("sham", "left DLPFC"), c("right DLPFC", "left DLPFC") )
myplot1 <- myplot1+ stat_compare_means(method="wilcox.test", paired=TRUE, comparison = my_comparisons, label = "p.format") 
myplot1 <- myplot1+xlab("Stimulation condition") + ylab("Estimated risk aversion (r)") +  labs(fill = "Stimulation")
myplot1 <- myplot1 + geom_point() 
#myplot1

#remove legend from figure
myplot1 <- myplot1 + theme(legend.position='none')

####### here we add adjusted p-values instead of unadjusted

myplot1_2 <- ggplot_build(myplot1) 
#myplot1_2$data

  
fun0 <- function(x){
  
  if (x=="sham-right DLPFC-1"){return (as.character(stat.test[1,5])) }
  if (x=="sham-left DLPFC-2"){return (as.character(stat.test[2,5])) }
  if (x=="right DLPFC-left DLPFC-3"){return (as.character(stat.test[3,5])) }
}  

myplot1_2$data[[2]]$annotation <- apply( as.matrix(paste(myplot1_2$data[[2]]$group)), MARGIN=1, FUN=fun0)

fun1 <- function(x) {
  return(symnum(x, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
  
}
myplot1_2$data[[2]]$annotation <- fun1(as.numeric(as.matrix(paste(myplot1_2$data[[2]]$annotation))))
#plot(ggplot_gtable(myplot1_2))
figure1 <- as.ggplot(ggplot_gtable(myplot1_2))


########################################################
#probability weighting
df <- data.frame()
a = coefs[which(is.element(coefs[,1], K) & (coefs[,2]==0) ), 4] 
b = coefs[which(is.element(coefs[,1], K) & (coefs[,2]==1) ), 4]
c = coefs[which(is.element(coefs[,1], K) & (coefs[,2]==2) ), 4]

df1 <- data.frame()
df2 <- data.frame()
df3 <- data.frame()

df1 <- data.frame(data=a, session="sham")
df2 <- data.frame(data=b, session="right DLPFC")
df3 <- data.frame(data=c, session="left DLPFC")
df <- rbind(df1, df2, df3)  

stat.test <- compare_means(data~session, method="wilcox.test", data= df, paired=TRUE, p.adjust.method="bonferroni")


myplot1 <- ggboxplot(df,  x="session", y="data", fill = "session", palette = c("#00AFBB", "#E7B800", "#FC4E07"))
myplot1 <- myplot1+theme_light()
my_comparisons <- list( c("sham", "right DLPFC"), c("sham", "left DLPFC"), c("right DLPFC", "left DLPFC") )
myplot1 <- myplot1+ stat_compare_means(method="wilcox.test", paired=TRUE, comparison = my_comparisons, label = "p.format") 
myplot1 <- myplot1+xlab("Stimulation condition") + ylab("Estimated probability weighting parameter") +  labs(fill = "Stimulation")
myplot1 <- myplot1 + geom_point()+theme(legend.position='none')
#myplot1

#remove legend from figure
myplot1 <- myplot1 + theme(legend.position='none')

####### here we add adjusted p-values instead of unadjusted

myplot1_3 <- ggplot_build(myplot1) 
#myplot1_3$data


fun0 <- function(x){
  
  if (x=="sham-right DLPFC-1"){return (as.character(stat.test[1,5])) }
  if (x=="sham-left DLPFC-2"){return (as.character(stat.test[2,5])) }
  if (x=="right DLPFC-left DLPFC-3"){return (as.character(stat.test[3,5])) }
}  

myplot1_3$data[[2]]$annotation <- apply( as.matrix(paste(myplot1_2$data[[2]]$group)), MARGIN=1, FUN=fun0)

fun1 <- function(x) {
  return(symnum(x, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
  
}
myplot1_3$data[[2]]$annotation <- fun1(as.numeric(as.matrix(paste(myplot1_3$data[[2]]$annotation))))
figure2 <- as.ggplot(ggplot_gtable(myplot1_3))
#plot(ggplot_gtable(myplot1_3))



#joint plots together

ggarrange(figure1, figure2, ncol=2, nrow=1, labels=c('A', 'B'))



########################################################
#creating tables for Kahneman-Tversky specification
########################################################

df1 <- data.frame(data=coefs[which(coefs[,2]==0 ), 3], session="sham", model="Model 1")
df2 <- data.frame(data=coefs[which(coefs[,2]==1 ), 3], session="right DLPFC", model="Model 1")
df3 <- data.frame(data=coefs[which(coefs[,2]==2 ), 3], session="left DLPFC", model="Model 1")
df <- rbind(df1, df2, df3)
df_table <- compare_means(data~session, method="wilcox.test", data= df, paired=TRUE, p.adjust.method="bonferroni")
Treatments <- c("right vs. sham", "left vs. sham", "left vs. right")
median_diff <- c( median (df2$data-df1$data), median (df3$data-df1$data), median (df3$data-df2$data))
numberobs <- paste0("n=", length(K))

df_table1 <- data.frame(Treatments, median_diff=round(median_diff, digits=3), pval=signif(df_table$p, digits=3), pvaladj=signif(df_table$p.adj, digits=3), numberobs )
#colnames(df_table1) <- c("Treatments", "Median difference in RA (r)", "Wilcoxon p-value", "Wilcoxon adjusted p-value", "Number of obs.")
df_table1 

df4 <- data.frame(data=coefs[which(coefs[,2]==0 ), 4], session="sham", model="Model 1")
df5 <- data.frame(data=coefs[which(coefs[,2]==1 ), 4], session="right DLPFC", model="Model 1")
df6 <- data.frame(data=coefs[which(coefs[,2]==2 ), 4], session="left DLPFC", model="Model 1")
df <- rbind(df4, df5, df6)
df_table <- compare_means(data~session, method="wilcox.test", data= df, paired=TRUE,  p.adjust.method="bonferroni")
Treatments <- c("right vs. sham", "left vs. sham", "left vs. right")
median_diff <- c( median (df5$data-df4$data), median (df6$data-df4$data), median (df6$data-df5$data))
numberobs <- paste0("n=", length(K))


df_table2 <- data.frame(Treatments, median_diff=round(median_diff, digits=3), pval=signif(df_table$p, digits=3), pvaladj=signif(df_table$p.adj, digits=3), numberobs )
#colnames(df_table2) <- c("Treatments", "Median difference in gamma", "Wilcoxon p-value", "Wilcoxon adjusted p-value", "Number of obs.")
df_table2 

#df_table <- rbind(df_table1, c("Treatments", "Median difference in gamma", "Wilcoxon p-value", "Wilcoxon adjusted p-value", "Number of obs.") ,df_table2)
df_names1 <-data.frame(Treatments="Parameter: Risk aversion (r)", median_diff="", pval="", pvaladj="", numberobs="", stringsAsFactors = FALSE)
df_names2 <-data.frame(Treatments="Parameter: gamma", median_diff="", pval="", pvaladj="", numberobs="", stringsAsFactors = FALSE)
df_table_ALL <- rbind(df_names1, df_table1, df_names2, df_table2)
colnames(df_table_ALL) <- c("Treatments", "Median difference in parameter", "Wilcoxon p-value", "Wilcoxon adjusted p-value", "Number of obs.")

df_table_ALL

formattable( df_table_ALL, align=c("l", "c", "c", "c"), caption="Kahneman-Tversky specifcation")







########################################################
### Estimation procedures for Prelec-1 specification
########################################################


N=(1:28)

coefs=c()

# j is for subjects
# k is for sessions

gains = 1  #set gains = -1 for loss domain


for (j in N) {
  for (k in c(0,1,2)) {
    data1 <- data [ which (subject == j & gains*outcome1>0 & session == k),] 
    log_lik <- function(pars) {
      r <- pars[1]
      alpha <- pars[2]
      mu <- pars[3]
      
      p <- data1$prob1
      outcome1 <- data1$outcome1
      outcome2 <- data1$outcome2
      outcome3 <- data1$outcome3
      outcome4 <- data1$outcome4
      
      EU <- c()
      
      for(i in 1:nrow(data1)) {
        
        if (gains == 1){
        EU_a_i <- (exp(-(-log(p[i]))^alpha)) * ((( outcome1[i])^(1-r))/(1-r)) + 
          (1 - (exp(-(-log(p[i]))^alpha))) *((( outcome2[i])^(1-r))/(1-r)) 
        
        EU_b_i <- (exp(-(-log(p[i]))^alpha)) * ((( outcome3[i])^(1-r))/(1-r))  + 
          (1 - (exp(-(-log(p[i]))^alpha))) * ((( outcome4[i])^(1-r))/(1-r))
        
        EU <- c(EU, 1/(1 + (EU_b_i/EU_a_i)^(1/mu)))   }
       
        else {
          EU_a_i <- (1-exp(-(-log(1-p[i]))^alpha)) * (-((- outcome1[i])^(1 - r))) / (1 - r) + 
            ((exp(-(-log(1-p[i]))^alpha))) * (-(( - outcome2[i])^(1 - r))) / (1 - r) 
          
          EU_b_i <- (1-exp(-(-log(1-p[i]))^alpha)) * (-((- outcome3[i])^(1 - r))) / (1 - r) + 
            ((exp(-(-log(1-p[i]))^alpha))) * (-(( - outcome4[i])^(1 - r))) / (1 - r)

          EU <- c(EU, 1/(1 + (EU_a_i/EU_b_i)^(1/mu)))  
        }
        
        
         
      }
      
      

      df <- data.frame( EU = EU, y = data1$choice)
      
      df$optim <- 0
      
      for(i in 1:nrow(data1)) {
        if(df$y[i] == 1) {
          df$optim[i] <- log(EU[i])
        } else {
          df$optim[i] <- log(1 - EU[i])
        }
      }
      
      return(-sum(df$optim))

    }
    
    theta = c(0.0001, 0.0001, 0.0001)
    par <- nmk(theta, f = log_lik)$par
    
    par <- cbind(j, k, t(par))   
    coefs=rbind(coefs, par)
  }
}

colnames(coefs) <- c("subject","session", "RA", "alpha", "mu")




#removing outliers 

s0 = IQR(coefs[which(coefs[,2]==0),3])  #computes interquantile range
s1 = IQR(coefs[which(coefs[,2]==1),3])
s2 = IQR(coefs[which(coefs[,2]==2),3])

s0_25=quantile(coefs[which(coefs[,2]==0),3],1/4)  #computes 25st quantile
s0_75=quantile(coefs[which(coefs[,2]==0),3],3/4)

s1_25=quantile(coefs[which(coefs[,2]==1),3],1/4)  #computes 25st quantile
s1_75=quantile(coefs[which(coefs[,2]==1),3],3/4)

s2_25=quantile(coefs[which(coefs[,2]==2),3],1/4)  #computes 25st quantile
s2_75=quantile(coefs[which(coefs[,2]==2),3],3/4)

RA_0_lower_fence = s0_25-3*s0
RA_0_upper_fence = s0_75+3*s0

RA_1_lower_fence = s1_25-3*s1
RA_1_upper_fence = s1_75+3*s1

RA_2_lower_fence = s2_25-3*s2
RA_2_upper_fence = s2_75+3*s2

s0 = IQR(coefs[which(coefs[,2]==0),4])  #computes interquantile range
s1 = IQR(coefs[which(coefs[,2]==1),4])
s2 = IQR(coefs[which(coefs[,2]==2),4])


s0_25=quantile(coefs[which(coefs[,2]==0),4],1/4)  #computes 25st quantile
s0_75=quantile(coefs[which(coefs[,2]==0),4],3/4)

s1_25=quantile(coefs[which(coefs[,2]==1),4],1/4)  #computes 25st quantile
s1_75=quantile(coefs[which(coefs[,2]==1),4],3/4)

s2_25=quantile(coefs[which(coefs[,2]==2),4],1/4)  #computes 25st quantile
s2_75=quantile(coefs[which(coefs[,2]==2),4],3/4)

alpha_0_lower_fence = s0_25-3*s0
alpha_0_upper_fence = s0_75+3*s0

alpha_1_lower_fence = s1_25-3*s1
alpha_1_upper_fence = s1_75+3*s1

alpha_2_lower_fence = s2_25-3*s2
alpha_2_upper_fence = s2_75+3*s2


s0_25=quantile(coefs[which(coefs[,2]==0),5],1/4)  #computes 25st quantile
s0_75=quantile(coefs[which(coefs[,2]==0),5],3/4)

s1_25=quantile(coefs[which(coefs[,2]==1),5],1/4)  #computes 25st quantile
s1_75=quantile(coefs[which(coefs[,2]==1),5],3/4)

s2_25=quantile(coefs[which(coefs[,2]==2),5],1/4)  #computes 25st quantile
s2_75=quantile(coefs[which(coefs[,2]==2),5],3/4)

mu_0_lower_fence = s0_25-3*s0
mu_0_upper_fence = s0_75+3*s0

mu_1_lower_fence = s1_25-3*s1
mu_1_upper_fence = s1_75+3*s1

mu_2_lower_fence = s2_25-3*s2
mu_2_upper_fence = s2_75+3*s2


Fences_matrix = matrix( c(RA_0_lower_fence, RA_1_lower_fence, RA_2_lower_fence, RA_0_upper_fence, RA_1_upper_fence, RA_2_upper_fence, alpha_0_lower_fence, alpha_1_lower_fence, alpha_2_lower_fence, alpha_0_upper_fence, alpha_1_upper_fence, alpha_2_upper_fence, mu_0_lower_fence, mu_1_lower_fence, mu_2_lower_fence, mu_0_upper_fence, mu_1_upper_fence, mu_2_upper_fence), ncol=6, nrow=3)

K = c()

for (i in N){
  coefs1 = coefs[which(coefs[,1]==i),]
  if (  ( all(coefs1[,3] > Fences_matrix[,1]) ==TRUE) & ( all(coefs1[,3] < Fences_matrix[,2]) ==TRUE) & ( all(coefs1[,4] > Fences_matrix[,3]) ==TRUE) & ( all(coefs1[,4] < Fences_matrix[,4]) ==TRUE) & ( all(coefs1[,5] > Fences_matrix[,5]) ==TRUE) & ( all(coefs1[,5] < Fences_matrix[,6]) ==TRUE) ) {K = c(K, i)} 
}


#computing difference in coefficients
coefs_diff_10=c()
coefs_diff_20=c()
coefs_diff_21=c()

for (j in K){
  
  coefs_diff_10 = rbind(coefs_diff_10, coefs[which(coefs[,1]==j & coefs[,2]==1),] - coefs[which(coefs[,1]==j & coefs[,2]==0),] )
  coefs_diff_20 = rbind(coefs_diff_20, coefs[which(coefs[,1]==j & coefs[,2]==2),] - coefs[which(coefs[,1]==j & coefs[,2]==0),] )
  coefs_diff_21 = rbind(coefs_diff_21, coefs[which(coefs[,1]==j & coefs[,2]==2),] - coefs[which(coefs[,1]==j & coefs[,2]==1),] )
  
}

coefs <- coefs[which(is.element(coefs[,1], K) ), ]  #this takes into account removing outliers



########################################################
#creating plots for Prelec-1 specification
########################################################
#risk aversion

df <- data.frame()
a = coefs[which(is.element(coefs[,1], K) & (coefs[,2]==0) ), 3] 
b = coefs[which(is.element(coefs[,1], K) & (coefs[,2]==1) ), 3]
c = coefs[which(is.element(coefs[,1], K) & (coefs[,2]==2) ), 3]

df1 <- data.frame()
df2 <- data.frame()
df3 <- data.frame()

df1 <- data.frame(data=a, session="sham")
df2 <- data.frame(data=b, session="right DLPFC")
df3 <- data.frame(data=c, session="left DLPFC")
df <- rbind(df1, df2, df3)  

stat.test <- compare_means(data~session, method="wilcox.test", data= df, paired=TRUE, p.adjust.method="bonferroni")

myplot1 <- ggboxplot(df,  x="session", y="data", fill = "session", palette = c("#00AFBB", "#E7B800", "#FC4E07"))
myplot1 <- myplot1+theme_light()
my_comparisons <- list( c("sham", "right DLPFC"), c("sham", "left DLPFC"), c("right DLPFC", "left DLPFC") )
myplot1 <- myplot1+ stat_compare_means(method="wilcox.test", paired=TRUE, comparison = my_comparisons, label = "p.format") 
myplot1 <- myplot1+xlab("Stimulation condition") + ylab("Estimated risk aversion (r)") +  labs(fill = "Stimulation")
myplot1 <- myplot1 + geom_point() 
#myplot1

#remove legend from figure
myplot1 <- myplot1 + theme(legend.position='none')

####### here we add adjusted p-values instead of unadjusted

myplot1_2 <- ggplot_build(myplot1) 
#myplot1_2$data


fun0 <- function(x){
  
  if (x=="sham-right DLPFC-1"){return (as.character(stat.test[1,5])) }
  if (x=="sham-left DLPFC-2"){return (as.character(stat.test[2,5])) }
  if (x=="right DLPFC-left DLPFC-3"){return (as.character(stat.test[3,5])) }
}  

myplot1_2$data[[2]]$annotation <- apply( as.matrix(paste(myplot1_2$data[[2]]$group)), MARGIN=1, FUN=fun0)

fun1 <- function(x) {
  return(symnum(x, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
  
}
myplot1_2$data[[2]]$annotation <- fun1(as.numeric(as.matrix(paste(myplot1_2$data[[2]]$annotation))))
figure1 <- as.ggplot(ggplot_gtable(myplot1_2))

#plot(ggplot_gtable(myplot1_2))


########################################################
#probability weighting
df <- data.frame()
a = coefs[which(is.element(coefs[,1], K) & (coefs[,2]==0) ), 4] 
b = coefs[which(is.element(coefs[,1], K) & (coefs[,2]==1) ), 4]
c = coefs[which(is.element(coefs[,1], K) & (coefs[,2]==2) ), 4]

df1 <- data.frame()
df2 <- data.frame()
df3 <- data.frame()

df1 <- data.frame(data=a, session="sham")
df2 <- data.frame(data=b, session="right DLPFC")
df3 <- data.frame(data=c, session="left DLPFC")
df <- rbind(df1, df2, df3)  

stat.test <- compare_means(data~session, method="wilcox.test", data= df, paired=TRUE, p.adjust.method="bonferroni")


myplot1 <- ggboxplot(df,  x="session", y="data", fill = "session", palette = c("#00AFBB", "#E7B800", "#FC4E07"))
myplot1 <- myplot1+theme_light()
my_comparisons <- list( c("sham", "right DLPFC"), c("sham", "left DLPFC"), c("right DLPFC", "left DLPFC") )
myplot1 <- myplot1+ stat_compare_means(method="wilcox.test", paired=TRUE, comparison = my_comparisons, label = "p.format") 
myplot1 <- myplot1+xlab("Stimulation condition") + ylab("Estimated probability weighting parameter") +  labs(fill = "Stimulation")
myplot1 <- myplot1 + geom_point()+theme(legend.position='none')
#myplot1

#remove legend from figure
myplot1 <- myplot1 + theme(legend.position='none')

####### here we add adjusted p-values instead of unadjusted

myplot1_3 <- ggplot_build(myplot1) 
#myplot1_3$data


fun0 <- function(x){
  
  if (x=="sham-right DLPFC-1"){return (as.character(stat.test[1,5])) }
  if (x=="sham-left DLPFC-2"){return (as.character(stat.test[2,5])) }
  if (x=="right DLPFC-left DLPFC-3"){return (as.character(stat.test[3,5])) }
}  

myplot1_3$data[[2]]$annotation <- apply( as.matrix(paste(myplot1_2$data[[2]]$group)), MARGIN=1, FUN=fun0)

fun1 <- function(x) {
  return(symnum(x, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
  
}
myplot1_3$data[[2]]$annotation <- fun1(as.numeric(as.matrix(paste(myplot1_3$data[[2]]$annotation))))
figure2 <- as.ggplot(ggplot_gtable(myplot1_3))

#plot(ggplot_gtable(myplot1_3))



#joint plots together
ggarrange(figure1, figure2, ncol=2, nrow=1, labels=c('A', 'B'))


########################################################
#creating tables for Prelec-1 specification
########################################################

df1 <- data.frame(data=coefs[which(coefs[,2]==0 ), 3], session="sham", model="Model 2")
df2 <- data.frame(data=coefs[which(coefs[,2]==1 ), 3], session="right DLPFC", model="Model 2")
df3 <- data.frame(data=coefs[which(coefs[,2]==2 ), 3], session="left DLPFC", model="Model 2")
df <- rbind(df1, df2, df3)
df_table <- compare_means(data~session, method="wilcox.test", data= df, paired=TRUE, p.adjust.method="bonferroni")
Treatments <- c("right vs. sham", "left vs. sham", "left vs. right")
median_diff <- c( median (df2$data-df1$data), median (df3$data-df1$data), median (df3$data-df2$data))
numberobs <- paste0("n=", length(K))

df_table1 <- data.frame(Treatments, median_diff=round(median_diff, digits=3), pval=signif(df_table$p, digits=3), pvaladj=signif(df_table$p.adj, digits=3), numberobs )
#colnames(df_table1) <- c("Treatments", "Median difference in RA (r)", "Wilcoxon p-value", "Wilcoxon adjusted p-value", "Number of obs.")
df_table1 

df4 <- data.frame(data=coefs[which(coefs[,2]==0 ), 4], session="sham", model="Model 2")
df5 <- data.frame(data=coefs[which(coefs[,2]==1 ), 4], session="right DLPFC", model="Model 2")
df6 <- data.frame(data=coefs[which(coefs[,2]==2 ), 4], session="left DLPFC", model="Model 2")
df <- rbind(df4, df5, df6)
df_table <- compare_means(data~session, method="wilcox.test", data= df, paired=TRUE,  p.adjust.method="bonferroni")
Treatments <- c("right vs. sham", "left vs. sham", "left vs. right")
median_diff <- c( median (df5$data-df4$data), median (df6$data-df4$data), median (df6$data-df5$data))
numberobs <- paste0("n=", length(K))


df_table2 <- data.frame(Treatments, median_diff=round(median_diff, digits=3), pval=signif(df_table$p, digits=3), pvaladj=signif(df_table$p.adj, digits=3), numberobs )
#colnames(df_table2) <- c("Treatments", "Median difference in gamma", "Wilcoxon p-value", "Wilcoxon adjusted p-value", "Number of obs.")
df_table2 

#df_table <- rbind(df_table1, c("Treatments", "Median difference in gamma", "Wilcoxon p-value", "Wilcoxon adjusted p-value", "Number of obs.") ,df_table2)
df_names1 <-data.frame(Treatments="Parameter: Risk aversion (r)", median_diff="", pval="", pvaladj="", numberobs="", stringsAsFactors = FALSE)
df_names2 <-data.frame(Treatments="Parameter: gamma", median_diff="", pval="", pvaladj="", numberobs="", stringsAsFactors = FALSE)
df_table_ALL <- rbind(df_names1, df_table1, df_names2, df_table2)
colnames(df_table_ALL) <- c("Treatments", "Median difference in parameter", "Wilcoxon p-value", "Wilcoxon adjusted p-value", "Number of obs.")

df_table_ALL

formattable( df_table_ALL, align=c("l", "c", "c", "c"), caption="Prelec-1 specifcation")







