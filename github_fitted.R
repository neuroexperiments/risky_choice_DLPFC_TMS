library(rstan)
library(hBayesDM)
library(bayestestR)
library(bayesplot)

#setwd()


#fit <- readRDS()  #a variable containing fitted model result
fitname='fit'
sumfit <- rstan::summary(fit)
r <- rownames(sumfit$summary)
r1 <- which(r=='mu_rho')
r2 <- which(r=='mu_delta2_prw')
sumfit$summary[r1:r2,"Rhat"]
sumfit$summary[r1:r2,"n_eff"]

sims <- as.matrix(fit) #gives a matrix of posterior draws
#computing R-hats and effective sample size
r <- colnames(sims)
r1 <- which(r=='mu_rho')
r2 <- which(r=='mu_delta2_prw')
cols = c("mu_rho", "mu_tau", "mu_prw", "mu_delta1_rho", "mu_delta1_tau", "mu_delta1_prw", "mu_delta2_rho", "mu_delta2_tau", "mu_delta2_prw")
Rh <- c()
Bulk <- c()
Tail <- c()
for (i in 1:length(cols)) {
  Rh <- cbind(Rh, Rhat(sims[,cols[i]]))
  Bulk <- cbind(Bulk, ess_bulk(sims[,cols[i]]))
  Tail <- cbind(Tail, ess_tail(sims[,cols[i]]))
}
Rh
Bulk
Tail


### plotting traceplots for visual analysis
pars <- data.frame(rho0= parvals[['mu_rho']], prw0=parvals[['mu_prw']], tau0 = parvals[['mu_tau']], rho1=parvals[['mu_delta1_rho']], prw1=parvals[['mu_delta1_prw']], tau1=parvals[['mu_delta1_tau']], rho2=parvals[['mu_delta2_rho']], prw2=parvals[['mu_delta2_prw']], tau2=parvals[['mu_delta2_tau']])
View(sumfit$summary[r1:r2,])
mean(sumfit$summary[r1:r2,"Rhat"])
check_hmc_diagnostics(fit)
traceplot(fit, pars=r[r1:r2])

### calculating WAIC 
colVars <- function (a){
  diff <- a - matrix (colMeans(a), nrow(a), ncol(a), byrow=TRUE)
  vars <- colMeans (diff^2)*nrow(a)/(nrow(a)-1)
  return (vars)
}

waic <- function (stanfit){
  log_lik <- rstan::extract (stanfit, "log_lik")$log_lik
  lppd <- sum (log (colMeans(exp(log_lik))))
  p_waic_1 <- 2*sum (log(colMeans(exp(log_lik))) - colMeans(log_lik))
  p_waic_2 <- sum (colVars(log_lik))
  waic_2 <- -2*lppd + 2*p_waic_2
  return (list (waic=waic_2, p_waic=p_waic_2, lppd=lppd, p_waic_1=p_waic_1))
}


### calculating posterior predictive checks

parvals <- rstan::extract(fit, permuted=TRUE)
pred0=parvals[['y_pred']]
raw_data <- data [ which (data$outcome1>0 ),] 

LB<- c()
HB <- c()

v1 <- sample.int(32000, 8000, replace=FALSE)

mean_ppp <- c()
for (v in v1){
  # curSubj <- subjList[i]
  #tmp <- raw_data [ which (raw_data$subject!=14 & raw_data$subject!=16 ), ]
  tmp <- raw_data [ which (raw_data$subject!=14 & raw_data$subject!=16), ] #excluding two participants with high interhemispheric difference in motor thresholds
  choice0 <- tmp[1:nrow(tmp), "choice"]
  choice<- (choice0-2)*(-1)
  ppp_ij <- c()
  pred=pred0[v,,]
  pred <- as.vector(t(pred))
  for (j in 1:nrow(tmp)){
    if (choice[j]==pred[j]){ppp_ij <- cbind(ppp_ij, 1)  
    } else{ppp_ij <- cbind(ppp_ij, 0)}
  }
  
  mean_ppp <- cbind(mean_ppp, mean(ppp_ij))
}

LB = rbind(LB, HDIofMCMC(mean_ppp)[1])
HB = rbind(HB, HDIofMCMC(mean_ppp)[2])

cbind(LB, HB)
plotHDI(t(mean_ppp))
posterior=data.frame(ppp=t(mean_ppp))
names(posterior)[1]="Mean PPP"


### plotting posterior distributions

parvals <- parvalsgain
posteriors = data.frame(mu_rho=parvals[['mu_rho']], 
                        mu_tau=parvals[['mu_tau']], 
                        mu_prw=parvals[['mu_prw']], 
                        mu_delta1_rho=parvals[['mu_delta1_rho']],
                        mu_delta1_tau=parvals[['mu_delta1_tau']],
                        mu_delta1_prw=parvals[['mu_delta1_prw']],
                        mu_delta2_rho=parvals[['mu_delta2_rho']],
                        mu_delta2_tau=parvals[['mu_delta2_tau']],
                        mu_delta2_prw=parvals[['mu_delta2_prw']])
post = as.list(posteriors)
library(stringr) #to wrap long text in labels
labels= c("Risk aversion","Prob. weighting", "Consistency", "\u0394 risk aversion", "\u0394 prob. weighting", "\u0394 consistency","\u0394 risk aversion", "\u0394 prob. weighting", "\u0394 consistency" )
plot_labels = data.frame( vars = names(post)[1:9], labels=str_wrap(labels, width=20))

library(bayesplot)
draw_CI <- function(x){
  par_plot <- mcmc_dens(as.data.frame(post[[x]]))
  ci_int_89 <- bayestestR::ci(post[[x]], ci=0.89, method='ETI')
  ci_int_95 <- bayestestR::ci(post[[x]], ci=0.95, method='ETI')
  xlab0 <- plot_labels[x,2]
  par_plot <-par_plot+geom_segment(aes(x=ci_int_89[,2],xend=ci_int_89[,3],y=0,yend=0), colour="black", size=4)+geom_segment(aes(x=ci_int_95[,2],xend=ci_int_95[,3],y=0,yend=0), colour="black", size=1)+yaxis_text(on=TRUE)+ylab("density") + xlab(xlab0)
  s = par_plot
  return(s)
}

my_color_scheme <- c("skyblue", "skyblue1",
                     "skyblue2", "skyblue3",
                     "skyblue4", "lightcyan4")
color_scheme_set(my_color_scheme)

pars = as.list(c(1:9))
posterior_plots <- lapply(pars, FUN=draw_CI)


