// Building up model from STAN Manual

functions {
  int first_capture(int[] y_i) {
    for (k in 1:size(y_i))
      if (y_i[k])
        return k;
    return 0;
  }
  int last_capture(int[] y_i) {
    for (k_rev in 0:(size(y_i) - 1)) {
      int k;
      k = size(y_i) - k_rev;
      if (y_i[k])
        return k;
    }
    return 0;
  }
  vector prob_uncaptured(int T, vector p, real phi) {
    vector[T] chi;
    chi[T] = 1.0;
    for (t in 1:(T - 1)) {
      int t_curr;
      int t_next;
      t_curr = T - t;
      t_next = t_curr + 1;
      chi[t_curr] = (1 - phi)
                     + phi
                       * (1 - p[t_next])
                       * chi[t_next];
    }
    return chi;
  }
}

data {
  int T; // # of time points
  int<lower=0> I; // # of individuals 
  int<lower=0,upper=1> y[I, T]; // # capture histories
  vector[T] temp; // temperature
  vector[T] prec; // precipitation
  int<lower = 1> n_missing;
  int<lower = 1> missing[n_missing];
  int<lower = 1> n_observed;
  int<lower = 1> observed[n_observed];
}

transformed data {
  int<lower=0,upper=T> first[I];
  int<lower=0,upper=T> last[I];
  vector<lower=0,upper=I>[T] n_captured;
  for (i in 1:I)
    first[i] = first_capture(y[i]);
  for (i in 1:I)
    last[i] = last_capture(y[i]);
  n_captured = rep_vector(0, T);
  for (t in 1:T)
    for (i in 1:I)
      if (y[i, t])
        n_captured[t] = n_captured[t] + 1;
}

parameters {
  real<lower=0,upper=1> phi;
  vector<lower = 0,upper = 1>[T] free_p;
}

transformed parameters {
  vector<lower=0,upper=1>[T] chi;
  vector<lower = 0,upper = 1>[T] p;
  
  for(i in missing){
    p[i] = 0;
  }
  
  for(i in observed){
    p[i] = free_p[i];
  }
  
  chi = prob_uncaptured(T,p,phi);
}

model {
  for (i in 1:I) {
    if (first[i] > 0) {
      for (t in (first[i]+1):last[i]) {
        1 ~ bernoulli(phi);
        y[i, t] ~ bernoulli(p[t]);
      }
      1 ~ bernoulli(chi[last[i]]);
    }
  }
}


generated quantities{  
vector[I] log_lik; // I think I'm missing the climate sensitive phi likelihood 
for (i in 1:I) {
    log_lik[i] = 0;
    if (first[i] > 0) {
      for (t in (first[i] + 1):last[i]) {
        log_lik[i] = log_lik[i] + bernoulli_lpmf(1|phi);
        log_lik[i] = log_lik[i] + bernoulli_lpmf(y[i, t]|p[t - 1]);
      }
      log_lik[i] = log_lik[i] + bernoulli_lpmf(1|chi[last[i]]);
    }
  // log_lik[i] = log_lik[i] + normal_lpdf(phi | mu, .2);
  }

}

