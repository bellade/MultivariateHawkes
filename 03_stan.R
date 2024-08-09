
#### prior on K  ####

stanblock_function_EXP_approxInt <- "
  
  real kernel(real x, real K, real beta){ 
    
    real res;
   
    res = K * beta * exp(-beta * x);
    
    return(res);
    
  }
  
  
  real intensityHawkes(real x, int x_m, real[] ts, int[] ts_m, real[] mu, real[,] K, real[,] beta){
    
    int N_ts;
    
    int m_i;
    
    real temp;
    real temp2[2];
    
    N_ts = size(ts);
    
    temp = mu[x_m];
    
    
    for(j in 1:N_ts){
      
      if(ts[j] < x){
        m_i  = ts_m[j];
        temp = temp + kernel(x - ts[j], K[m_i, x_m], beta[m_i, x_m]);
      }
    }
    
    temp2[1] = 0;
    temp2[2] = temp;
      
  return(max(temp2));
  
  }
  
  real approxInterval(real a, real b, int M, real[] ts, int[] ts_m, real[] mu, real[,] K, real[,] beta){
  
    real res;
    real int_middle;
    // real int_a;
    // real int_b;
    
    res = 0;
    
    for(m in 1:M){
      int_middle = intensityHawkes((a+b)/2, m, ts, ts_m, mu, K, beta);
      // int_a = intensityHawkes(a+0.001, m, ts, ts_m, mu, K, beta);
      // int_b = intensityHawkes(b, m, ts, ts_m, mu, K, beta);
      // res = res + ((int_a + int_b)/2) * (b-a);
      
      res = res + int_middle * (b-a);
    }
    
    return(res);
  }
  
  
  
  real approxIntegral(real[] ts, int[] ts_m, real[] mu, real[,] K, real[,] beta, real maxTime, int M){
  
    int N_ts;
    real res;
    
    N_ts = size(ts);
    
    res = approxInterval(0, ts[1], M, ts, ts_m, mu, K, beta);
    
    for(i in 2:N_ts){
      res = res + approxInterval(ts[i-1], ts[i], M, ts, ts_m, mu, K, beta);
    }
    
    res = res + approxInterval(ts[N_ts], maxTime, M, ts, ts_m, mu, K, beta);
    
    return(res);
    
    
  }
  
 
  real Hawkeslikelihood(real[] ts, int[] ts_m, real[] mu, real[,] K, real[,] beta, real maxTime, int M) {
    
    int N_ts;

    real loglik; 
    real temp; 
    real tt;
    int tt_m;
  
    N_ts = size(ts);

    
    loglik = 0;

    for (i in 1:N_ts) {
    
        tt = ts[i]; 
        tt_m = ts_m[i];
        loglik = loglik + log(intensityHawkes(tt, tt_m, ts, ts_m, mu, K, beta)); 
        
    }
    
    
   loglik = loglik - approxIntegral(ts, ts_m, mu, K, beta, maxTime, M);
   
   
    return(loglik);
  }
  
  "


Stan_fit_multi_Model1_EXP_exin_approxInt  <- function(ts, M, maxTime, ncores = 3, 
                                                      prior_mu = c(0, 3), prior_K = c(0, 0.25),
                                                      iter = 1000, chains = 3){
  # uniform priors for mu and beta
  # fancy priors for K
  ## fits the Hawkes process model (assuming complete data)
  
  # there are only two betas (on and off diagonal, i.e. self and other)
  
  # approximates integral (triangular between each obs)
  
  # data        event times of a Hawkes process
  # maxTime     upper limit for the event times
  # ncores      number of cores used for the estimation
  # prior_...   lower and upper bound for the uniform priors of the respective parameter
  
  
  scode_hawkes <- paste0("functions {", stanblock_function_EXP_approxInt, "} \n", "
  
  data {
    int<lower=0> N;
    int<lower=0> M;
    real<lower=0> ts[N];
    int<lower=0> ts_m[N];
    real<lower=0> maxTime;
    real<lower=0> prior_mu1;
    real<lower=0> prior_mu2;
    real<lower=-1> prior_K1;
    real<lower=-1> prior_K2;
  }
      
  parameters {
   real<lower=0> mu[M];
   real<lower=-1,upper=1> K[M,M];
   real<lower=0,upper=1> beta_vec[2];
  }
  
  transformed parameters{
  
    real<lower=0,upper=1> beta[M,M];
    
    for(m1 in 1:M){
      for(m2 in 1:M){
        
        if(m1 == m2){
          beta[m1, m2] = beta_vec[1];
        } else {
          beta[m1, m2] = beta_vec[2];
        }
        
      }
    }
  }
  
  
  model {
    
    for(m in 1:M){
      mu[m] ~ uniform(prior_mu1,prior_mu2);
    }
    
    
    for(m1 in 1:M){
       for(m2 in 1:M){
         K[m1, m2] ~ normal(prior_K1, prior_K2);
       }
     }
    
       
    for(i in 1:2){
      beta_vec[i] ~ uniform(0, 1);
    }
     
    target += Hawkeslikelihood(ts, ts_m, mu, K, beta, maxTime, M);
  }") 
  
  
  # true data
  hawkes_data <- list(N =length(ts[, 1]),
                      M = M,
                      ts = ts[, 1], 
                      ts_m= ts[, 2], 
                      maxTime = maxTime,
                      prior_mu1 = prior_mu[1], prior_mu2 = prior_mu[2],
                      prior_K1 = prior_K[1], prior_K2 = prior_K[2]) 
  fit_h <- stan(model_code = scode_hawkes, iter=iter, verbose=FALSE,
                data=hawkes_data, chains = chains, cores = ncores)
  fit <- list(mu = extract(fit_h)$mu, K = extract(fit_h)$K, beta = extract(fit_h)$beta_vec)
  return(fit)
  
}







#### prior on K Star ####


##### constant background rate #####


stanblock_function_EXP_approxInt_capat_Kstar <- "

  real[,] Kstar_to_K(real[,] Kstar, int M){
    
    real K[M,M];
    matrix[M,M] K_temp;
    matrix[M,M] Kstar_temp;
    matrix[M,M] diagM;
    
    for(m1 in 1:M){
      for(m2 in 1:M){
        Kstar_temp[m1, m2] = Kstar[m1, m2];
      }
    }
    
    diagM = diag_matrix(rep_vector(1.0, M));
    
    K_temp =  diagM - inverse(Kstar_temp + diagM);
    
    for(m1 in 1:M){
      for(m2 in 1:M){
        K[m1, m2] = K_temp[m1, m2];
      }
    }
    
    return(K);
  }
  
  real[,] K_to_Kstar(real[,] K, int M){
    
    real Kstar[M,M];
    matrix[M,M] K_temp;
    matrix[M,M] Kstar_temp;
    matrix[M,M] diagM;
    
    for(m1 in 1:M){
      for(m2 in 1:M){
        K_temp[m1, m2] = K[m1, m2];
      }
    }
    
    diagM = diag_matrix(rep_vector(1.0, M));
    
    Kstar_temp =  inverse(diagM - K_temp) - diagM;
    
    for(m1 in 1:M){
      for(m2 in 1:M){
        Kstar[m1, m2] = Kstar_temp[m1, m2];
      }
    }
    
    return(Kstar);
  }
  
  real kernel(real x, real K, real beta){ 
    
    real res;
   
    res = K * beta * exp(-beta * x);
    
    return(res);
    
  }
  
  
  real intensityHawkes(real x, int x_m, real[] ts, int[] ts_m, real[] mu, real[,] K, real[,] beta, real capat){
    
    int N_ts;
    
    int m_i;
    
    real temp;
    real temp2[2];
    
    N_ts = size(ts);
    
    temp = mu[x_m];
    
    
    for(j in 1:N_ts){
      
      if(ts[j] < x){
        m_i  = ts_m[j];
        temp = temp + kernel(x - ts[j], K[m_i, x_m], beta[m_i, x_m]);
      }
    }
    
    temp2[1] = capat;
    temp2[2] = temp;
      
  return(max(temp2));
  
  }
  
  real approxInterval(real a, real b, int M, real[] ts, int[] ts_m, real[] mu, real[,] K, real[,] beta, real capat){
  
    real res;
    real int_middle1;
    real int_middle2;
    real int_a;
    real int_b;
    real temp;
    
    res = 0;
    
    for(m in 1:M){
     
      int_a = intensityHawkes(a+0.0001, m, ts, ts_m, mu, K, beta, capat);
      int_b = intensityHawkes(b, m, ts, ts_m, mu, K, beta, capat);
      int_middle1 = intensityHawkes((2*a+b)/3, m, ts, ts_m, mu, K, beta, capat);
      int_middle2 = intensityHawkes((a+2*b)/3, m, ts, ts_m, mu, K, beta, capat);
      
      temp = ((b-a)/8) * (int_a + 3*int_middle1 + 3*int_middle2 + int_b); 
      
      res = res + temp;
      
     
    }
    
    return(res);
  }
  
  
  
  real approxIntegral(real[] ts, int[] ts_m, real[] mu, real[,] K, real[,] beta, real maxTime, int M, real capat){
  
    int N_ts;
    real res;
    
    N_ts = size(ts);
    
    res = approxInterval(0, ts[1], M, ts, ts_m, mu, K, beta, capat);
    
    for(i in 2:N_ts){
      res = res + approxInterval(ts[i-1], ts[i], M, ts, ts_m, mu, K, beta, capat);
    }
    
    res = res + approxInterval(ts[N_ts], maxTime, M, ts, ts_m, mu, K, beta, capat);
    
    return(res);
    
    
  }
  
 
  real Hawkeslikelihood(real[] ts, int[] ts_m, real[] mu, real[,] K, real[,] beta, real maxTime, int M, real capat) {
    
    int N_ts;

    real loglik; 
    real temp; 
    real tt;
    int tt_m;
  
    N_ts = size(ts);

    
    loglik = 0;

    for (i in 1:N_ts) {
    
        tt = ts[i]; 
        tt_m = ts_m[i];
        loglik = loglik + log(intensityHawkes(tt, tt_m, ts, ts_m, mu, K, beta, capat)); 
        
    }
    
    
   loglik = loglik - approxIntegral(ts, ts_m, mu, K, beta, maxTime, M, capat);
   
   
    return(loglik);
  }
  
  "







Stan_fit_multi_Model1_EXP_exin_approxInt_capat_Kstar  <- function(ts, M, maxTime, capat = 0.01,
                                                                  ncores = 3, 
                                                                  prior_mu = c(0, 3), prior_K = c(0, 0.25), prior_beta = c(0, 1),
                                                                  iter = 1000, chains = 3){

  # uniform priors for mu and beta
  # normal prior on K*
  ## fits the Hawkes process model
  
  # approximates integral 
  
  # data        event times of a Hawkes process
  # maxTime     upper limit for the event times
  # ncores      number of cores used for the estimation
  # prior_...   lower and upper bound for the uniform priors of the respective parameter
  
  
  scode_hawkes <- paste0("functions {", stanblock_function_EXP_approxInt_capat_Kstar, "} \n", "
  data {
    int<lower=0> N;
    int<lower=0> M;
    real<lower=0> ts[N];
    int<lower=0> ts_m[N];
    real<lower=0> maxTime;
    real<lower=0> capat;
    real<lower=0> prior_mu1;
    real<lower=0> prior_mu2;
    real<lower=-1> prior_K1;
    real<lower=-1> prior_K2;
    real<lower=-1> prior_beta1;
    real<lower=-1> prior_beta2;
  }
      
  parameters {
   real<lower=0> mu[M];
   real Kstar[M,M];
   real<lower=0,upper=1> beta_vec[2];
  }
  
  transformed parameters{
  
    real<lower=0,upper=1> beta[M,M];
    real K[M,M];
    real<upper=1> K_diag[M];
    
    for(m1 in 1:M){
      for(m2 in 1:M){
        
        if(m1 == m2){
          beta[m1, m2] = beta_vec[1];
        } else {
          beta[m1, m2] = beta_vec[2];
        }
        
      }
    }
    
    K = Kstar_to_K(Kstar, M);
    
    for(m in 1:M){
      K_diag[m] = K[m, m];
    }
  }
  
  
  model {
    
    for(m in 1:M){
      mu[m] ~ uniform(prior_mu1,prior_mu2);
    }
    
    
    for(m1 in 1:M){
       for(m2 in 1:M){
         Kstar[m1, m2] ~ normal(prior_K1, prior_K2);
       }
     }
    
       
    for(i in 1:2){
      beta_vec[i] ~ uniform(prior_beta1, prior_beta2);
    }
     
    target += Hawkeslikelihood(ts, ts_m, mu, K, beta, maxTime, M, capat);
  }") 
  
  
  # true data
  hawkes_data <- list(N =length(ts[, 1]),
                      M = M,
                      ts = ts[, 1], 
                      ts_m= ts[, 2], 
                      maxTime = maxTime,
                      capat = capat,
                      prior_mu1 = prior_mu[1], prior_mu2 = prior_mu[2],
                      prior_K1 = prior_K[1], prior_K2 = prior_K[2],
                      prior_beta1 = prior_beta[1], prior_beta2 = prior_beta[2]) 
  fit_h <- stan(model_code = scode_hawkes, iter=iter, verbose=FALSE,
                data=hawkes_data, chains = chains, cores = ncores)
  extract_fit <- extract(fit_h)
  fit <- list(mu = extract_fit$mu, K = extract_fit$K, beta = extract_fit$beta_vec,  Kstar = extract_fit$Kstar )
  
  return(fit)
  
}



##### background rate has a deterministic part #####



stanblock_function_EXP_approxInt_capat_Kstar_mulookup <- "

  real mu_lookup(real t, real[] mu_lookup_values){
  
    int interval;
    int N;
    real res;
    
    interval = 1;
    N = size(mu_lookup_values);
    
    for(j in 1:(N-1)){
      if(t >= j){
        interval = interval + 1;
        }
    }
    
    res = mu_lookup_values[interval];
    return(res);
    
    
  }

  real[,] Kstar_to_K(real[,] Kstar, int M){
    
    real K[M,M];
    matrix[M,M] K_temp;
    matrix[M,M] Kstar_temp;
    matrix[M,M] diagM;
    
    for(m1 in 1:M){
      for(m2 in 1:M){
        Kstar_temp[m1, m2] = Kstar[m1, m2];
      }
    }
    
    diagM = diag_matrix(rep_vector(1.0, M));
    
    K_temp =  diagM - inverse(Kstar_temp + diagM);
    
    for(m1 in 1:M){
      for(m2 in 1:M){
        K[m1, m2] = K_temp[m1, m2];
      }
    }
    
    return(K);
  }
  
  real[,] K_to_Kstar(real[,] K, int M){
    
    real Kstar[M,M];
    matrix[M,M] K_temp;
    matrix[M,M] Kstar_temp;
    matrix[M,M] diagM;
    
    for(m1 in 1:M){
      for(m2 in 1:M){
        K_temp[m1, m2] = K[m1, m2];
      }
    }
    
    diagM = diag_matrix(rep_vector(1.0, M));
    
    Kstar_temp =  inverse(diagM - K_temp) - diagM;
    
    for(m1 in 1:M){
      for(m2 in 1:M){
        Kstar[m1, m2] = Kstar_temp[m1, m2];
      }
    }
    
    return(Kstar);
  }
  
  real kernel(real x, real K, real beta){ 
    
    real res;
   
    res = K * beta * exp(-beta * x);
    
    return(res);
    
  }
  
  
  real intensityHawkes(real x, int x_m, real[] ts, int[] ts_m, real[] v, real[] mu_lookup_values, real[,] K, real[,] beta, real capat){
    
    int N_ts;
    
    int m_i;
    
    real temp;
    real temp2[2];
    
    N_ts = size(ts);
    
    temp = v[x_m] * mu_lookup(x, mu_lookup_values);
    

    for(j in 1:N_ts){
      
      if(ts[j] < x){
        m_i  = ts_m[j];
        temp = temp + kernel(x - ts[j], K[m_i, x_m], beta[m_i, x_m]);
      }
    }
    
    temp2[1] = capat;
    temp2[2] = temp;
      
  return(max(temp2));
  
  }
  
  real approxInterval(real a, real b, int M, real[] ts, int[] ts_m, real[] v, real[] mu_lookup_values, real[,] K, real[,] beta, real capat){
  
    real res;
    real int_middle1;
    real int_middle2;
    real int_a;
    real int_b;
    real temp;
    
    res = 0;
    
    for(m in 1:M){
     
      int_a = intensityHawkes(a+0.0001, m, ts, ts_m, v, mu_lookup_values, K, beta, capat);
      int_b = intensityHawkes(b, m, ts, ts_m, v, mu_lookup_values, K, beta, capat);
      int_middle1 = intensityHawkes((2*a+b)/3, m, ts, ts_m, v, mu_lookup_values, K, beta, capat);
      int_middle2 = intensityHawkes((a+2*b)/3, m, ts, ts_m, v, mu_lookup_values, K, beta, capat);
      
      temp = ((b-a)/8) * (int_a + 3*int_middle1 + 3*int_middle2 + int_b); 
      
      res = res + temp;
      
     
    }
    
    return(res);
  }
  
  
  
  real approxIntegral(real[] ts, int[] ts_m, real[] v, real[] mu_lookup_values, real[,] K, real[,] beta, real maxTime, int M, real capat){
  
    int N_ts;
    real res;
    
    N_ts = size(ts);
    
    res = approxInterval(0, ts[1], M, ts, ts_m, v, mu_lookup_values, K, beta, capat);
    
    for(i in 2:N_ts){
      res = res + approxInterval(ts[i-1], ts[i], M, ts, ts_m, v, mu_lookup_values, K, beta, capat);
    }
    
    res = res + approxInterval(ts[N_ts], maxTime, M, ts, ts_m, v, mu_lookup_values, K, beta, capat);
    
    return(res);
    
    
  }
  
 
  real Hawkeslikelihood(real[] ts, int[] ts_m, real[] v, real[] mu_lookup_values, real[,] K, real[,] beta, real maxTime, int M, real capat) {
    
    int N_ts;

    real loglik; 
    real temp; 
    real tt;
    int tt_m;
  
    N_ts = size(ts);

    
    loglik = 0;

    for (i in 1:N_ts) {
    
        tt = ts[i]; 
        tt_m = ts_m[i];
        loglik = loglik + log(intensityHawkes(tt, tt_m, ts, ts_m, v, mu_lookup_values, K, beta, capat)); 
        
    }
    
    
   loglik = loglik - approxIntegral(ts, ts_m, v, mu_lookup_values, K, beta, maxTime, M, capat);
   
   
    return(loglik);
  }
  
  "


Stan_fit_multi_Model1_EXP_exin_approxInt_capat_Kstar_mulookup  <- function(ts, M, mu_lookup_values, maxTime, capat = 0.01,
                                                                           ncores = 3, 
                                                                           prior_mu = c(0, 3), prior_K = c(0, 0.25), prior_beta = c(0, 1),
                                                                           iter = 1000, chains = 3){
  # uniform priors for mu and beta
  # fancy priors for K
  ## fits the Hawkes process model (assuming complete data)
  
  # approximates integral (triangular between each obs)
  
  # data        event times of a Hawkes process
  # maxTime     upper limit for the event times
  # ncores      number of cores used for the estimation
  # prior_...   lower and upper bound for the uniform priors of the respective parameter
  
  
  scode_hawkes <- paste0("functions {", stanblock_function_EXP_approxInt_capat_Kstar_mulookup, "} \n", "
  data {
    int<lower=0> N;
    int<lower=0> M;
    int<lower=0> maxTime_int;
    real<lower=0> ts[N];
    int<lower=0> ts_m[N];
    real mu_lookup_values[maxTime_int];
    real<lower=0> maxTime;
    real<lower=0> capat;
    real<lower=0> prior_mu1;
    real<lower=0> prior_mu2;
    real<lower=-1> prior_K1;
    real<lower=-1> prior_K2;
    real<lower=-1> prior_beta1;
    real<lower=-1> prior_beta2;
  }
      
  parameters {
   real<lower=0> v[M];
   real<upper=2> Kstar[M,M];
   real<lower=0> beta_vec[2];
  }
  
  transformed parameters{
  
    real<lower=0> beta[M,M];
    real K[M,M];
    real<upper=1> K_diag[M];
    
    for(m1 in 1:M){
      for(m2 in 1:M){
        
        if(m1 == m2){
          beta[m1, m2] = beta_vec[1];
        } else {
          beta[m1, m2] = beta_vec[2];
        }
        
      }
    }
    

    K = Kstar_to_K(Kstar, M);
    
    for(m in 1:M){
      K_diag[m] = K[m, m];
    }
  }
  
  
  model {
    
    for(m in 1:M){
     v[m] ~ normal(prior_mu1,prior_mu2);
    }
    
    
    for(m1 in 1:M){
       for(m2 in 1:M){
         Kstar[m1, m2] ~ normal(prior_K1, prior_K2);
       }
     }
    
       
    for(i in 1:2){
      beta_vec[i] ~ uniform(prior_beta1, prior_beta2);
    }
     
    target += Hawkeslikelihood(ts, ts_m, v, mu_lookup_values, K, beta, maxTime, M, capat);
  }") 
  
  
  # true data
  hawkes_data <- list(N =length(ts[, 1]),
                      M = M,
                      maxTime_int = maxTime,
                      ts = ts[, 1], 
                      ts_m= ts[, 2], 
                      mu_lookup_values = mu_lookup_values,
                      maxTime = maxTime,
                      capat = capat,
                      prior_mu1 = prior_mu[1], prior_mu2 = prior_mu[2],
                      prior_K1 = prior_K[1], prior_K2 = prior_K[2],
                      prior_beta1 = prior_beta[1], prior_beta2 = prior_beta[2]) 
  fit_h <- stan(model_code = scode_hawkes, iter=iter, verbose=FALSE,
                data=hawkes_data, chains = chains, cores = ncores)
  extract_fit <- extract(fit_h)
  fit <- list(v = extract_fit$v, K = extract_fit$K, beta = extract_fit$beta_vec,  Kstar = extract_fit$Kstar )
  
  return(fit)
  
}







Stan_fit_multi_Model1_EXP_exONLY_approxInt_capat_Kstar_mulookup  <- function(ts, M, mu_lookup_values, maxTime, capat = 0.01,
                                                                             ncores = 3, 
                                                                             prior_mu = c(0, 3), prior_K = c(0, 0.25), prior_beta = c(0, 1),
                                                                             iter = 1000, chains = 3){
  # excitation only (no inhibition)
  
  # uniform priors for mu and beta
  # fancy priors for K
  ## fits the Hawkes process model (assuming complete data)
  
  # approximates integral (triangular between each obs)
  
  # data        event times of a Hawkes process
  # maxTime     upper limit for the event times
  # ncores      number of cores used for the estimation
  # prior_...   lower and upper bound for the uniform priors of the respective parameter
  
  
  scode_hawkes <- paste0("functions {", stanblock_function_EXP_approxInt_capat_Kstar_mulookup, "} \n", "
  data {
    int<lower=0> N;
    int<lower=0> M;
    int<lower=0> maxTime_int;
    real<lower=0> ts[N];
    int<lower=0> ts_m[N];
    real mu_lookup_values[maxTime_int];
    real<lower=0> maxTime;
    real<lower=0> capat;
    real<lower=0> prior_mu1;
    real<lower=0> prior_mu2;
    real<lower=-1> prior_K1;
    real<lower=-1> prior_K2;
    real<lower=-1> prior_beta1;
    real<lower=-1> prior_beta2;
  }
      
  parameters {
   real<lower=0> v[M];
   real<lower=0, upper=0.75> Kstar[M,M];
   real<lower=0> beta_vec[2];
  }
  
  transformed parameters{
  
    real<lower=0> beta[M,M];
    real<lower=-0.025,upper=1.5> K[M,M];
    real<upper=1> K_diag[M];
    
    for(m1 in 1:M){
      for(m2 in 1:M){
        
        if(m1 == m2){
          beta[m1, m2] = beta_vec[1];
        } else {
          beta[m1, m2] = beta_vec[2];
        }
        
      }
    }
    

    K = Kstar_to_K(Kstar, M);
    
    for(m in 1:M){
      K_diag[m] = K[m, m];
    }
  }
  
  
  model {
    
    for(m in 1:M){
     v[m] ~ normal(prior_mu1,prior_mu2);
    }
    
    
    for(m1 in 1:M){
       for(m2 in 1:M){
         Kstar[m1, m2] ~ normal(prior_K1, prior_K2);
       }
     }
    
       
    for(i in 1:2){
      beta_vec[i] ~ uniform(prior_beta1, prior_beta2);
    }
     
    target += Hawkeslikelihood(ts, ts_m, v, mu_lookup_values, K, beta, maxTime, M, capat);
  }") 
  
  
  # true data
  hawkes_data <- list(N =length(ts[, 1]),
                      M = M,
                      maxTime_int = maxTime,
                      ts = ts[, 1], 
                      ts_m= ts[, 2], 
                      mu_lookup_values = mu_lookup_values,
                      maxTime = maxTime,
                      capat = capat,
                      prior_mu1 = prior_mu[1], prior_mu2 = prior_mu[2],
                      prior_K1 = prior_K[1], prior_K2 = prior_K[2],
                      prior_beta1 = prior_beta[1], prior_beta2 = prior_beta[2]) 
  fit_h <- stan(model_code = scode_hawkes, iter=iter, verbose=FALSE,
                data=hawkes_data, chains = chains, cores = ncores)
  extract_fit <- extract(fit_h)
  fit <- list(v = extract_fit$v, K = extract_fit$K, beta = extract_fit$beta_vec,  Kstar = extract_fit$Kstar )
  
  return(fit)
  
}



Stan_fit_multi_Model1_EXP_backgrONLY_approxInt_capat_mulookup  <- function(ts, M, mu_lookup_values, maxTime, capat = 0.01,
                                                                           ncores = 3, 
                                                                           prior_mu = c(0, 3),
                                                                           iter = 1000, chains = 3){
  # scaled background rate only (no excitation nor inhibition)
  
  # uniform priors for mu and beta
  # fancy priors for K
  ## fits the Hawkes process model (assuming complete data)
  
  # approximates integral (triangular between each obs)
  
  # data        event times of a Hawkes process
  # maxTime     upper limit for the event times
  # ncores      number of cores used for the estimation
  # prior_...   lower and upper bound for the uniform priors of the respective parameter
  
  K <- matrix(0, M, M)
  beta <-  matrix(0.5, M, M)
  
  
  scode_hawkes <- paste0("functions {", stanblock_function_EXP_approxInt_capat_Kstar_mulookup, "} \n", "
  data {
    int<lower=0> N;
    int<lower=0> M;
    int<lower=0> maxTime_int;
    real<lower=0> ts[N];
    int<lower=0> ts_m[N];
    real mu_lookup_values[maxTime_int];
    real<lower=0> maxTime;
    real<lower=0> capat;
    real prior_mu1;
    real<lower=0> prior_mu2;
    real K[M,M];
    real<lower=0> beta[M,M];
  }
      
  parameters {
   real<lower=0> v[M];
  }

  model {
    
    for(m in 1:M){
     v[m] ~ normal(prior_mu1,prior_mu2);
    }
    
     
    target += Hawkeslikelihood(ts, ts_m, v, mu_lookup_values, K, beta, maxTime, M, capat);
  }") 
  
  
  # true data
  hawkes_data <- list(N =length(ts[, 1]),
                      M = M,
                      maxTime_int = maxTime,
                      ts = ts[, 1], 
                      ts_m= ts[, 2], 
                      mu_lookup_values = mu_lookup_values,
                      maxTime = maxTime,
                      capat = capat,
                      prior_mu1 = prior_mu[1], prior_mu2 = prior_mu[2],
                      K = K,
                      beta = beta) 
  fit_h <- stan(model_code = scode_hawkes, iter=iter, verbose=FALSE,
                data=hawkes_data, chains = chains, cores = ncores)
  extract_fit <- extract(fit_h)
  fit <- list(v = extract_fit$v )
  
  return(fit)
  
}





