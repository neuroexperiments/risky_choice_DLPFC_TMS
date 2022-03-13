data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsubj[N];
  real<lower=0> outcome1[N, T];
  real<lower=0> outcome2[N, T];
  real<lower=0> outcome3[N, T];
  real<lower=0> outcome4[N, T];
  real<lower=0, upper=1> probs[N, T]; //probability of best outcome
  int<lower=0, upper=1> choice[N, T];
  int session [N,T];
  
}

parameters {
  // Hyper(group)-parameters
  real mu_pr[9];
  // subject-level raw parameters
  real rho_pr[N];
  real tau_pr[N];
  real prw_pr[N];
  real delta1_rho_pr[N];
  real delta1_tau_pr[N];
  real delta1_prw_pr[N];
  real delta2_rho_pr[N];
  real delta2_tau_pr[N];
  real delta2_prw_pr[N];
  real sigma_normal[9]; 
}

transformed parameters {
   real rho0[N];
   real tau0[N];
   real prw0[N];
   real delta1_rho[N];
   real delta1_tau[N];
   real delta1_prw[N];
   real delta2_rho[N];
   real delta2_tau[N];
   real delta2_prw[N];
   real sigma[9];
   
    for (i in 1:9){
      sigma[i] = exp(0+3*sigma_normal[i]);
    }
    
    for (i in 1:N){
      rho0[i] = 0+5*Phi_approx(mu_pr[1]+sigma[1]*rho_pr[i]);
      tau0[i]  = exp(mu_pr[2]+sigma[2]*tau_pr[i]); 
      prw0[i] = 6*Phi_approx(mu_pr[3]+sigma[3]*prw_pr[i]);
      
      delta1_rho[i] = -2+4*Phi_approx(mu_pr[4]+sigma[4]*delta1_rho_pr[i]);
      delta1_tau[i] = -5+10*Phi_approx(mu_pr[5]+sigma[5]*delta1_tau_pr[i]);
      delta1_prw[i] = -1.5+3*Phi_approx(mu_pr[6]+sigma[6]*delta1_prw_pr[i]);
      
      delta2_rho[i] = -2+4*Phi_approx(mu_pr[7]+sigma[7]*delta2_rho_pr[i]);
      delta2_tau[i] = -5+10*Phi_approx(mu_pr[8]+sigma[8]*delta2_tau_pr[i]);
      delta2_prw[i] = -1.5+3*Phi_approx(mu_pr[9]+sigma[9]*delta2_prw_pr[i]);
      
    }
    

model {
  // Hyper(group)-parameters
  mu_pr[1] ~normal(0, 5.0);
  mu_pr[2] ~normal(0, 10.0);
  mu_pr[3] ~normal(1, 3.0);
  mu_pr[4:9]  ~ normal(0, 1.0);

  sigma_normal ~ normal(0, 1);
  // Subject-level raw parameters (for Matt trick)
   rho_pr    ~ normal(0,1);
   tau_pr    ~ normal(0,1);
   prw_pr   ~ normal(0,1);
   delta1_rho_pr ~ normal(0,1);
   delta1_tau_pr ~ normal(0,1);
   delta1_prw_pr ~ normal(0,1);
   delta2_rho_pr ~ normal(0,1);
   delta2_tau_pr ~ normal(0,1);
   delta2_prw_pr ~ normal(0,1);

  for (i in 1:N) {
    for (t in 1:Tsubj[i]) {
        real EU_a_i;
        real EU_b_i;
        real pChooseA; //probability to choose lottery A over lottery B
        real prw;
        real rho;
        real tau;

        if (session[i,t]==0) {
        rho = rho0[i];
        tau = tau0[i];
        prw = prw0[i];
        
      } else if (session[i,t]==1){
        rho = rho0[i]+delta1_rho[i];
        tau = tau0[i]+delta1_tau[i];
        prw = prw0[i]+delta1_prw[i];
        
      } else {
        rho = rho0[i]+delta2_rho[i];
        tau = tau0[i]+delta2_tau[i];
        prw = prw0[i]+delta2_prw[i];
      }

      EU_a_i = ((probs[i, t]^prw)/(((probs[i, t])^prw+(1-probs[i, t])^prw)^(1/prw))) * ((outcome1[i, t])^rho)+ 
          (1 - ((probs[i, t]^prw)/(((probs[i, t])^prw+(1-probs[i, t])^prw)^(1/prw)))) *((outcome2[i, t])^rho);
        
      EU_b_i = ((probs[i, t]^prw)/(((probs[i, t])^prw+(1-probs[i, t])^prw)^(1/prw))) *((outcome3[i, t])^rho)+ 
          (1 - ((probs[i, t]^prw)/(((probs[i, t])^prw+(1-probs[i, t])^prw)^(1/prw)))) * ((outcome4[i, t])^rho);
      
         pChooseA  = tau * (log( ((EU_a_i)^(1/rho))/((EU_b_i)^(1/rho))));
         choice[i, t] ~ bernoulli_logit(pChooseA);
      
      
    }
  }
}
generated quantities {
  real<lower=0, upper=5> mu_rho;
  real<lower=0> mu_tau;
  real<lower=0, upper=6> mu_prw;
  real<lower=-2, upper=2> mu_delta1_rho;
  real<lower=-5, upper=5> mu_delta1_tau;
  real<lower=-1.5, upper=1.5> mu_delta1_prw;
  real<lower=-2, upper=2> mu_delta2_rho;
  real<lower=-5, upper=5> mu_delta2_tau;
  real<lower=-1.5, upper=1.5> mu_delta2_prw;
  int j_t;

  real log_lik[N*T];

  // For posterior predictive check
  real y_pred[N, T];

  for (i in 1:N) {
    for (t in 1:T) {
      y_pred[i, t] = -1;
    }
  }

  mu_rho    = Phi_approx(mu_pr[1]) * 5;
  mu_tau    = exp(mu_pr[2]);
  mu_prw    = Phi_approx(mu_pr[3]) *6;
  
  mu_delta1_rho = Phi_approx(mu_pr[4]) * (4) -2 ;
  mu_delta1_tau = Phi_approx(mu_pr[5]) * (10) -5 ;
  mu_delta1_prw = Phi_approx(mu_pr[6]) * (3) -1.5 ;
  mu_delta2_rho = Phi_approx(mu_pr[7]) * (4) -2 ;
  mu_delta2_tau = Phi_approx(mu_pr[8]) * (10) -5 ;
  mu_delta2_prw = Phi_approx(mu_pr[9]) * (3) -1.5 ;
  
  j_t=1;
  { // local section, this saves time and space
    for (i in 1:N) {
      //log_lik[i] = 0;
      for (t in 1:Tsubj[i]) {
        real EU_a_i;
        real EU_b_i;
        real pChooseA; //probability to choose lottery A over lottery B
        real prw;
        real rho;
        real tau;
 
        if (session[i,t]==0) {
        rho = rho0[i];
        tau = tau0[i];
        prw = prw0[i];
        
      } else if (session[i,t]==1){
        rho = rho0[i]+delta1_rho[i];
        tau = tau0[i]+delta1_tau[i];
        prw = prw0[i]+delta1_prw[i];
        
      } else {
        rho = rho0[i]+delta2_rho[i];
        tau = tau0[i]+delta2_tau[i];
        prw = prw0[i]+delta2_prw[i];
      }
      
       EU_a_i = ((probs[i, t]^prw)/(((probs[i, t])^prw+(1-probs[i, t])^prw)^(1/prw))) * ((outcome1[i, t])^rho)+ 
          (1 - ((probs[i, t]^prw)/(((probs[i, t])^prw+(1-probs[i, t])^prw)^(1/prw)))) *((outcome2[i, t])^rho);
        
      EU_b_i = ((probs[i, t]^prw)/(((probs[i, t])^prw+(1-probs[i, t])^prw)^(1/prw))) *((outcome3[i, t])^rho)+ 
          (1 - ((probs[i, t]^prw)/(((probs[i, t])^prw+(1-probs[i, t])^prw)^(1/prw)))) * ((outcome4[i, t])^rho);
 
          pChooseA  = tau * (log( ((EU_a_i)^(1/rho))/((EU_b_i)^(1/rho))));
        
        log_lik[j_t] = bernoulli_logit_lpmf(choice[i, t] | pChooseA);
        j_t=j_t+1;
        // generate posterior prediction for current trial
        y_pred[i, t] = bernoulli_logit_rng(pChooseA);
      }
    }
  }
}
