// CJS Model with time varying survival (s) and capture probability (p)
    
data{
  int<lower = 1> NsumCH; // number of capture histories
  int<lower = 1> n_occasions; // number of capture occasions
  int<lower = 1, upper = 2> sumCH[NsumCH, n_occasions]; // capture histories (1,2)
  int<lower = 1> sumf[NsumCH]; // first capture occasion
  int<lower = 1> sumFR[NsumCH]; // frequency of each capture history
  int<lower = 1> n_missing; // number of unsampled years
  int<lower = 1> missing[n_missing]; // unsampled years
  int<lower = 1> n_observed; // number of sampled years
  int<lower = 1> observed[n_observed]; // sampled years
  vector[n_occasions-1] temp; // temperature each year
  vector[n_occasions-1] prec; // precipitation each year
  int n_lik; // number of data to calculate a likelihood for
  int<lower = 0, upper = 1> CH[NsumCH, n_occasions]; // capture histories (1,0)
  int S; // number of species
  int species[NsumCH]; // species associated with each capture history
}

parameters{
  real<lower = 0, upper = 1> free_p[S, n_occasions - 1]; // detection probability for all years
  vector[S] phi_0_raw; // intercept for survival model
  vector[S] phi_prec_raw; // effect of prec on survival
  real mu_0; // mean for phi_0 distribution
  real<lower = 0> sig_0; // sd for phi_0 distribution
  real mu_prec; // mean for phi_prec distribution
  real<lower = 0> sig_prec; // sd for phi_prec distribution
}

transformed parameters{
  simplex[2] tr[S, 2,n_occasions - 1]; // survival likelihoods for marginalization
  simplex[2] pmat[S, 2,n_occasions - 1]; // detection likelihoods for marginalization
  // real mu[n_occasions-1];
  vector[S] phi_0; // intercept for survival model
  vector[S] phi_prec; // effect of prec on survival
  real<lower = 0, upper = 1> p[S, n_occasions - 1]; // detection probability with unsampled years = 0
  real<lower=0,upper=1> phi[S, n_occasions-1]; // survival probability

  phi_0 = mu_0 + sig_0 * phi_0_raw;
  phi_prec = mu_prec + sig_prec * phi_prec_raw;
  
for(s in 1:S){
    for(i in missing){
    p[s, i] = 0; // set unsampled p to 0
  }
  
  for(i in observed){
    p[s, i] = free_p[s, i]; // fill in sampled p with a real estimate
  }
  
  for(i in 1:n_occasions - 1){
   phi[s, i] = inv_logit(phi_0[s]+prec[i]*phi_prec[s]); // survival is a function of predictor vars
  }
  
  for(k in 1:n_occasions - 1){
    tr[s,1,k,1] = phi[s,k];
    tr[s,1,k,2] = (1 - phi[s,k]);
    tr[s,2,k,1] = 0;
    tr[s,2,k,2] = 1;
    
    pmat[s,1,k,1] = p[s,k];
    pmat[s,1,k,2] = (1 - p[s,k]);
    pmat[s,2,k,1] = 0;
    pmat[s,2,k,2] = 1;
  }
}
}
  
model{
    vector[2] pz[S,n_occasions]; // marginalized likelihoods
    
    // priors
    phi_0_raw ~ std_normal();
    phi_prec_raw ~ std_normal();
    
    for(i in 1:NsumCH){  
      pz[species[i],sumf[i],1] = 1;
      pz[species[i],sumf[i],2] = 0;
      
      for(k in (sumf[i] + 1):n_occasions){ 
        pz[species[i],k,1] = pz[species[i],k-1,1] * tr[species[i],1,k-1,1] * pmat[species[i],1,k-1,sumCH[i,(k)]];
        pz[species[i],k,2] = (pz[species[i],k-1, 1] * tr[species[i],1,k-1,2] + pz[species[i],k-1, 2]) * pmat[species[i],2,k-1,sumCH[i,(k)]];
      }  
      
      target += sumFR[i] * log(sum(pz[species[i],n_occasions])); 
      
    }
    
}

generated quantities {
  vector[n_lik] log_lik;
  vector[2] pz[S,n_occasions]; // marginalized likelihoods
  vector[n_lik] y; // CH values
  vector[n_lik] y_hat; // predicted observation probabilities from pz

  int counter;
  
  counter = 1;
  for(i in 1:NsumCH){ // doesn't include critters observed first on last occasion
      pz[species[i],sumf[i],1] = 1;
      pz[species[i],sumf[i],2] = 0;
      
      for(k in (sumf[i] + 1):n_occasions){ 
        pz[species[i],k,1] = pz[species[i],k-1,1] * tr[species[i],1,k-1,1] * pmat[species[i],1,k-1,sumCH[i,(k)]];
        pz[species[i],k,2] = (pz[species[i],k-1, 1] * tr[species[i],1,k-1,2] + pz[species[i],k-1, 2]) * pmat[species[i],2,k-1,sumCH[i,(k)]];
      }
      
    if (sumf[i] >= n_occasions){
    }
    else{
          for(j in sumf[i]+1:n_occasions){
      for(k in 1:sumFR[i]){
        y[counter] = CH[i,j];
        y_hat[counter] = pz[species[i], j, sumCH[i,j]];
        log_lik[counter] = bernoulli_lpmf(CH[i,j]|phi[species[i],j-1]*p[species[i],j-1]);
        counter += 1;
      }
    }
    }
    }
    
}