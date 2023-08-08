//--- Simple SEIR model ---//

data {
  
  int Nt; // N time steps
  int pop; // population size
  int cases[Nt]; // reported cases
  real alpha; // incubation period
  real sigma; // infectious period
  real seed; // N initial infections
  vector[Nt] imports; // N imports throught time series

  
}


parameters {
  
  vector<lower=0>[Nt] Rt;
  real <lower=0> phi; // over-dispersion
  real <lower=0,upper=1> rho; // reporting rate
  
}


transformed parameters {
  
  vector<lower=0,upper=1>[Nt] S; // susceptible
  vector<lower=0,upper=1>[Nt] E; // exposed
  vector<lower=0,upper=1>[Nt] I; // infectious
  vector<lower=0,upper=1>[Nt] R; // recovered & seroneg
  vector[Nt] Bt = Rt*sigma;
  vector[Nt] inc; // incidence
  vector[Nt] ecases; // predicted cases
  
  // initial states
      S[1] = 1 - seed;
      E[1] = 0;
      I[1] = seed;
      R[1] = 0;
  
  // update compartments
  for(t in 2:Nt){
    S[t] = S[t-1] - S[t-1]*Bt[t-1]*I[t-1];
    E[t] = E[t-1] + S[t-1]*Bt[t-1]*I[t-1] - E[t-1]*alpha;
    I[t] = I[t-1] + E[t-1]*alpha - I[t-1]*sigma + imports[t];
    R[t] = R[t-1] + I[t-1]*sigma;
  }
  
  // expected cases
  for(t in 1:Nt) inc[t] = S[t]*Bt[t]*I[t];
  for(t in 1:Nt) ecases[t] = rho*pop*inc[t];

    
}


model {
  
  // Priors
  phi ~ normal(2,0.5);
  rho ~ beta(1,5);
  Rt[1] ~ normal(1.5,0.5);
  for(t in 2:Nt) Rt[t] ~ normal(Rt[t-1], 0.1);
 
  // model 
  cases ~ neg_binomial_2(ecases, phi);
  
}

generated quantities {
  
  vector[Nt] predcases;
  for(t in 1:Nt) predcases[t] = neg_binomial_2_rng(ecases[t], phi);
  
}




