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
  real<lower=0,upper=1> phi; // survival probability
  real<lower = 0, upper = 1> free_p[n_occasions - 1]; // detection probability for all years
}

transformed parameters{
  simplex[2] tr[2,n_occasions - 1]; // survival likelihoods for marginalization
  simplex[2] pmat[2,n_occasions - 1]; // detection likelihoods for marginalization
  // real mu[n_occasions-1];
  real<lower = 0, upper = 1> p[n_occasions - 1]; // detection probability with unsampled years = 0

  for(i in missing){
    p[i] = 0; // set unsampled p to 0
  }
  
  for(i in observed){
    p[i] = free_p[i]; // fill in sampled p with a real estimate
  }
  
  for(k in 1:n_occasions - 1){
    tr[1,k,1] = phi;
    tr[1,k,2] = (1 - phi);
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
  vector[2] pz[n_occasions]; // marginalized likelihoods
  vector[n_lik] y; // CH values
  vector[n_lik] y_hat; // predicted observation probabilities from pz

  int counter;
  
  counter = 1;
  for(i in 1:NsumCH){ // doesn't include critters observed first on last occasion
      pz[sumf[i],1] = 1;
      pz[sumf[i],2] = 0;
      
      for(k in (sumf[i] + 1):n_occasions){ 
        pz[k,1] = pz[k-1,1] * tr[1,k-1,1] * pmat[1,k-1,sumCH[i,(k)]];
        pz[k,2] = (pz[k-1, 1] * tr[1,k-1,2] + pz[k-1, 2]) * pmat[2,k-1,sumCH[i,(k)]];
      }
      
    if (sumf[i] >= n_occasions){
    }
    else{
          for(j in sumf[i]+1:n_occasions){
      for(k in 1:sumFR[i]){
        y[counter] = CH[i,j];
        y_hat[counter] = pz[j, sumCH[i,j]];
        log_lik[counter] = bernoulli_lpmf(CH[i,j]|phi*p[j-1]);
        counter += 1;
      }
    }
    }
    }
    
}