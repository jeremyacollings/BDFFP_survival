// CJS Model with time varying survival (s) and capture probability (p)
    
data{
  int<lower = 1> NsumCH;
  int<lower = 1> n_occasions;
  int<lower = 1, upper = 2> sumCH[NsumCH, n_occasions];
  int<lower = 1> sumf[NsumCH];
  int<lower = 1> sumFR[NsumCH];
  int<lower = 1> n_missing;
  int<lower = 1> missing[n_missing];
  int<lower = 1> n_observed;
  int<lower = 1> observed[n_observed];
  vector[n_occasions-1] temp;
  vector[n_occasions-1] prec;
}

parameters{
  real s_int;   // survival intercept
  real p_int;   // detection intercept
  real temp_s; // effect of temperature on survival
  real temp_p; // effect of temperature on detection
  real prec_s; // effect of precipitation on survival
  real prec_p; // effect of precipitation on detection
  real<lower = 0, upper = 1> p;
  real<lower = 0, upper = 1> s[n_occasions - 1];
}

transformed parameters{
  simplex[2] tr[2,n_occasions - 1];
  simplex[2] pmat[2,n_occasions - 1];
  vector[2] pz[n_occasions];
  vector[n_occasions] log_lik;
  
  // Add back in this whole missing data thing later
  // for(i in missing){
  //   p[i] = 1*(10^-10);
  // }
  
  for(k in 1:n_occasions - 1){
    tr[1,k,1] = s[k];
    tr[1,k,2] = (1 - s[k]);
    tr[2,k,1] = 0;
    tr[2,k,2] = 1;
    
    pmat[1,k,1] = p;
    pmat[1,k,2] = (1 - p);
    pmat[2,k,1] = 0;
    pmat[2,k,2] = 1;
  }
  
  for(i in 1:NsumCH){  
    pz[sumf[i],1] = 1;
    pz[sumf[i],2] = 0;
      
    for(k in (sumf[i] + 1):n_occasions){ 
      pz[k,1] = pz[k-1,1] * tr[1,k-1,1] * pmat[1,k-1,sumCH[i,(k)]];
      pz[k,2] = (pz[k-1, 1] * tr[1,k-1,2] + pz[k-1, 2]) * pmat[2,k-1,sumCH[i,(k)]];
    }  
    
    for(n in 1:n_occasions){
      log_lik[n] = sumFR[i] * log(sum(pz[n]));
    }
  }
}
  
model{
  for(n in 1:n_occasions){
    target += log_lik[n];
  }
}

    
