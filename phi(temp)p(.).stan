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
}

parameters{
  real<lower = 0, upper = 1> real_p; // detection probability for all years
  real phi_0; // intercept for survival model
  real phi_temp; // effect of temp on survival
}

transformed parameters{
  simplex[2] tr[2,n_occasions - 1]; // survival likelihoods for marginalization
  simplex[2] pmat[2,n_occasions - 1]; // detection likelihoods for marginalization
  // real mu[n_occasions-1];
  real<lower = 0, upper = 1> p[n_occasions - 1]; // detection probability with unsampled years = 0
  real<lower=0,upper=1> phi[n_occasions-1]; // survival probability

  for(i in missing){
    p[i] = 0; // set unsampled p to 0
  }
  
  for(i in observed){
    p[i] = real_p; // fill in sampled p with a real estimate
  }
  
  for(i in 1:n_occasions - 1){
   phi[i] = inv_logit(phi_0+temp[i]*phi_temp); // survival is a function of predictor vars
  }
  
  for(k in 1:n_occasions - 1){
    tr[1,k,1] = phi[k];
    tr[1,k,2] = (1 - phi[k]);
    tr[2,k,1] = 0;
    tr[2,k,2] = 1;
    
    pmat[1,k,1] = p[k];
    pmat[1,k,2] = (1 - p[k]);
    pmat[2,k,1] = 0;
    pmat[2,k,2] = 1;
  }
  }
  
model{
    vector[2] pz[n_occasions]; // marginalized likelihoods
    
    phi_0 ~ normal(0, 25);
    phi_temp ~ normal(0, 25);
    
    for(i in 1:NsumCH){  
      pz[sumf[i],1] = 1;
      pz[sumf[i],2] = 0;
      
      for(k in (sumf[i] + 1):n_occasions){ 
        pz[k,1] = pz[k-1,1] * tr[1,k-1,1] * pmat[1,k-1,sumCH[i,(k)]];
        pz[k,2] = (pz[k-1, 1] * tr[1,k-1,2] + pz[k-1, 2]) * pmat[2,k-1,sumCH[i,(k)]];
      }  
      
      target += sumFR[i] * log(sum(pz[n_occasions])); 
      
    }
    
}


generated quantities {
  vector[n_lik] log_lik;
  int counter;
  counter = 1;
  for(i in 1:NsumCH){ // doesn't include critters observed first on last occasion
    if (sumf[i] >= n_occasions){
    }
    else{
          for(j in sumf[i]+1:n_occasions){
      for(k in 1:sumFR[i]){
        log_lik[counter] = bernoulli_lpmf(CH[i,j]|phi[j-1]*p[j-1]);
        counter += 1;
      }
    }
    }
    }
}