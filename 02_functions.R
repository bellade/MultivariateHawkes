##### multivariate Hawkes Functions ###

excit_EXP <- function(t, K, beta){
  res <- K * beta * exp(- beta * t)
  return(res)
}


#### Pre ####

softplus <- function(x, s){
  
  # link function that maps any input x the positive real lines
  
  s * log(1 + exp(x/s))
}



mu_step_multi <- function(t, params){
  
  # params is a list with two entries
  # first entry is a vector of breakpoints (i.e. when the value of the intensity changes)
  # B.. number of breakpoints
  # M.. number of dimensions in the MVH
  # second entry is a matrix [B +1, M] with values for each interval per dimensions
  
  interval <- sum(t[1] >= params[[1]]) + 1
  return(params[[2]][t[2], interval])
  
}



mu_step_scale_single <- function(t, params){
  
  # params is a list with three entries
  # first entry is a vector of breakpoints (i.e. when the value of the intensity changes)
  # second is a vector of values of those points
  # third is a vector of sclae parameters
  
  interval <- sum(as.numeric(t[1]) >= params[[1]]) + 1
  res <- params[[3]][as.numeric(t[2])] * params[[2]][interval]
  return(res)
  
}



#### Intensity ####



intensity_EXP_multi <- function(t, ts, mu, K, beta, softplus = FALSE, s = NA, mu_function = FALSE, params = NA, capat = NA){
  
  # evaluate the intensity at time t where t is a vector of length 2: t = (timestamp, dimension)
  
  # ts...   all events in a Nx2 matrix where the first column has the time and the second the dimension as an integer
  # mu...   background rate as a M-dimensional vector
  # K,r,p...parameters for influence kernel, each a MxM matrix
  # s...    parameter for the link function
  # mu_function, params ... only needed when mu is a function
  
  N <- length(ts[, 1]) # number of ovservations in ts
  
  res <- ifelse(mu_function, mu(t, params), mu[as.numeric(t[2])]) # background rate
  
  m_t <- as.numeric(t[2]) #dimension that t is from
  
  for(i in 1:N){
    if(ts[i, 1] < t[1]){
      # we check for each event in ts is if happened before t
      # if it did, then it contributes to the intensity at t
      
      m_i <- as.numeric(ts[i, 2]) # dimension that the event in ts cam from
      res <- res + excit_EXP(t = t[1] - ts[i, 1], K = K[m_i, m_t], beta = beta[m_i, m_t]) # contribution to the intensity
    }
  }
  
  
  
  if(softplus){
    return(softplus(res, s = s)) # push through link function to make sure it is positive
  } else if (!is.na(capat)) {
    return(max(res, capat))
  } else {
    return(res)
  }
  
  
  
  
}




#### Integral #### 

approxInterval <- function(a, b, ts, mu = NA, K, beta, mu_function = FALSE, params = NA, capat, one_m = FALSE, which_one_m = NA){
  
  # approximation of the integral between two points a and b
  
  res <- 0
  
  M <- dim(K)[1]
  
  if(one_m){
    m_seq <- which_one_m
  } else {
    m_seq <- 1:M
  }
  
  for(m in m_seq){
    
    
    int_a <- intensity_EXP_multi(t = c(a+0.0001, m), ts = ts, mu = mu, K = K, beta = beta, mu_function = mu_function, params = params, capat = capat)
    int_b <- intensity_EXP_multi(t = c(b, m), ts = ts, mu = mu, K = K, beta = beta, mu_function = mu_function, params = params, capat = capat)
    int_middle1 <- intensity_EXP_multi(t = c((2*a+b)/3, m), ts = ts, mu = mu, K = K, beta = beta, mu_function = mu_function, params = params, capat = capat)
    int_middle2 <- intensity_EXP_multi(t = c((a+2*b)/3, m), ts = ts, mu = mu, K = K, beta = beta, mu_function = mu_function, params = params, capat = capat)
    
    temp <- ((b-a)/8) * (int_a + 3*int_middle1 + 3*int_middle2 + int_b) 
    
    res <- res + temp
    
    
  }
  
  return(res)
}


approxIntegral <- function(ts, maxTime, mu = NA, K, beta, mu_function = FALSE, params = NA, useprevious = FALSE, ts_previous = NA, 
                           capat, starttime_test = NA, one_m = FALSE, which_one_m = NA){
  
  # approximation of the integral across the whole [0, T]
  
  M <- dim(K)[1]
  
  N_ts = length(ts[, 1])
  
  if(useprevious){
    
    L <- length(ts_previous[, 1])
    ts <- rbind(ts_previous, ts)
    
  } else {
    
    L <- 0
    
  }
  
  if(N_ts == 0){
    
    dummy_ts <- cbind(c(maxTime + 1:3), rep(1, 3))
    res = approxInterval(a = ifelse(useprevious, starttime_test, 0), b = maxTime, ts = dummy_ts, mu = mu, K = K, beta = beta, mu_function = mu_function, params = params, 
                   capat = capat, one_m=one_m, which_one_m=which_one_m)
    
  } else {
    res = approxInterval(a = ifelse(useprevious, starttime_test, 0), b = ts[(L+1), 1], ts = ts, mu = mu, K = K, beta = beta, mu_function = mu_function, params = params, 
                         capat = capat, one_m=one_m, which_one_m=which_one_m)
    
    if(N_ts > 1){
      for(i in (L+2):(L+N_ts)){
        res = res + approxInterval(a = ts[i-1, 1], b = ts[i, 1], ts = ts, mu = mu, K = K, beta = beta, mu_function = mu_function, params = params,
                                   capat = capat, one_m=one_m, which_one_m=which_one_m)
      }
    }
    
    
    res = res + approxInterval(a = ts[(L+N_ts), 1], b = maxTime, ts = ts, mu = mu, K = K, beta = beta, mu_function = mu_function, params = params,
                               capat = capat, one_m=one_m, which_one_m=which_one_m)
    
    
  }
  
  
 
  return(res)
}



#### Likelihood ####



Likelihood_EXP_multi<- function(ts, mu, K, beta, maxTime, softplus = FALSE, s = NA, mu_function = FALSE, params = NA, approxInt = FALSE,
                                useprevious = FALSE, ts_previous = NA, capat, starttime_test = NA){
  
  if(!approxInt & !useprevious){ print("not implemented yet, please use approximate")}
  
  
  M <- dim(K)[2]
  
  N <- length(ts[, 1])
  res <- 0
  
  
  if(useprevious){
    
    ts_test <- ts
    L <- length(ts_previous[, 1])
    ts <- rbind(ts_previous, ts)
    
  } else {
    
    L <- 0
    ts_previous <- matrix(NA, 3, 3)
    
  }
  
  
  # evaluate the log intensity at each point of the data set and sum up
  for(i in (L+1):(L+N) ){
    res <- res + log(intensity_EXP_multi(t = ts[i, 1:2], ts = ts[, 1:2], mu = mu, K = K, beta = beta, s = s, softplus = softplus,
                                         mu_function = mu_function, params = params, capat = capat))
    
  }
  
  # integral (EXP only)
  
  if(approxInt == FALSE){ # does not work for previous
    for(m in 1:M){
      temp <- mu[m] * maxTime 
      
      for(i in 1:N){ 
        m_i <- ts[i, 2]
        temp <- temp + K[m_i, m] * (1 - exp(-beta[m_i, m] * (maxTime - ts[i, 1])))
      }
      
      res <- res - temp
    }
  } else {
    

    res <- res - approxIntegral(ts = ts[(L+1):(L+N), 1:2], maxTime = maxTime, mu = mu, K = K, beta = beta, mu_function = mu_function, params = params, 
                                useprevious = useprevious, ts_previous = ts_previous[, 1:2], 
                                capat = capat, starttime_test = starttime_test)
    
  }
  
  
  
  
  return(res)
}


#### Sampler ####



simulate_EXP_multi <- function(M, mu , K, beta, softplus = FALSE, s = NA, maxTime, mu_function = FALSE, params = NA, uppereval = 10, capat = 0) {
  
  # simulate a series of events using thinning
  
  # initialise the ts
  ts <- matrix(c(Inf, 0), ncol = 2) # initialise (will be removed in the end, has no effect)
  
  hawkeseval <- 0
  
  index <- 0 # index of accepted events
  t <- rep(0, times = M) #current time
  
  MM <- uppereval #initiate the intensity for the process with thin from
  
  while(min(t) <= maxTime) {
    #choose homogenous poisson intensity that is greater than max hawkes value (MM)
    
    t_old <- t #save the previous time step in case we need to get back to it
    
    t <- t + rexp(M, MM) # suggest new observation in each dimension
    t_m <- which.min(t) # dimension in which we can find the minimum
    
    hawkeseval <-  max(intensity_EXP_multi(t = c(min(t), t_m), ts = ts, mu = mu, K = K, beta = beta, softplus = softplus, s = s, 
                                           mu_function = mu_function, params = params), capat)  # evaluate the intensity at this minimum
    
    if(hawkeseval > MM){ stop("hawkeseval > MM")}
    
    if (runif(1) < hawkeseval/MM) {
      # usual thinning
      index <- index + 1
      ts <- rbind(ts, c(min(t), t_m))
    }
    
    t <- rep(min(t), M) 
    
    
  }
  
  return(ts[ts[, 1] < maxTime, ])
  # remove the first observation which was only there to initiate the ts
}



#### Maximum Likelihood ####

### General


maxlik_multi_general <- function(ts, M, r, maxTime, s, initval = NA, displayOutput=FALSE){
  
  # get the maxium likelihood estimate
  
  fn <- function(params, ts, M, r, maxTime, s) {
    
    mu <- params[1:M]
    
    K <- matrix(params[(M+1):(M^2 + M)], byrow = TRUE, ncol = M)
    
    p <- matrix(params[(M^2 + M+1):(M^2 + M + M^2)], byrow = TRUE, ncol = M)
    
    
    lik <- -Likelihood_NB_multi(ts = ts, M= M, mu = mu, K = K, p = p, r = r, s = s, maxTime = maxTime, mu_function = FALSE, params = NA)
    
    
    if (displayOutput==TRUE) {
      print(c(params,-lik))
    }
    return(lik)
  }
  
  if (is.na(initval) ) {initval <- c(rep(length(ts[, 1])/(maxTime*M),M), rep(0,M^2),  rep(0.5,M^2))}
  
  temp <- optim(initval, fn, method='L-BFGS-B',
                lower = c(rep(0.01,M), rep(-1,M^2),  rep(0.01,M^2)), 
                upper = c(rep(Inf,M), rep(1,M^2),  rep(1,M^2)),
                ts = ts, M = M, r = r, maxTime = maxTime, s = s) 
  temppar <- temp$par
  return(list(mu =  temppar[1:M], 
              K = matrix(temppar[(M+1):(M^2 + M)], byrow = TRUE, ncol = M), 
              p = matrix(temppar[(M^2 + M+1):(M^2 + M + M^2)], byrow = TRUE, ncol = M), 
              lik= -temp$value))
}






maxlik_multi_EXP_1 <- function(ts, M, maxTime, softplus = FALSE, s = NA, initval = NA, displayOutput=FALSE, maxit = 100){
  
  # get the maxium likelihood estimate when there are only two betas
  
  fn <- function(params, ts, M, mu, beta, maxTime, softplus, s) {
    
    mu <- params[1:M]
    
    K <- matrix(params[(M + 1):(M + M^2)], M, M)
    
    beta_vec  <- params[(M + M^2 + 1):(M + M^2 + 2)]
    
    beta <- matrix(beta_vec[2], M, M)
    diag(beta) <- beta_vec[1]
    
    lik <- -Likelihood_EXP_multi(ts = ts, M = M, mu = mu, K = K, beta = beta, softplus = softplus, s = s, maxTime = maxTime, mu_function = FALSE)
    
    if (displayOutput==TRUE) {
      print(c(params,-lik))
    }
    return(lik)
  }
  
  if (is.na(initval[1]) ) {initval <- c(rep(0.5, times = M), rep(0, times = M^2), rep(0.5, times = 2))}
  
  temp <- optim(par = initval, fn = fn, method='L-BFGS-B',
                lower = c(rep(0.01, M), rep(-0.95, M^2), rep(0,  2)) ,
                upper = c(rep(3, times = M), rep(0.95, times = M^2), rep(0.95, times = 2)), control = list( maxit = maxit), 
                ts = ts, M = M, maxTime = maxTime, softplus = softplus, s = s) 
  par <- temp$par
  return(list(mu = par[1:M],
              K = par[(M + 1):(M + M^2)],
              beta_vec= par[(M + M^2 + 1):(M + M^2 + 2)],
              lik= -temp$value))
}



#### For K Star ####


spectralnorm <- function(mat){
  sqrt(max(eigen(t(mat) %*% mat)$values))
}


K_to_Kstar <- function(K){
  M <- dim(K)[1]
  # t(solve(diag(M) - K) - diag(M))
  solve(diag(M) - K) - diag(M)
}



Kstar_to_K <- function(Kstar){
  M <- dim(Kstar)[1]
  diag(M) - solve(Kstar + diag(M))
}







#### background date #### 

mu_date_intensity <- function(x_ind, theta_mu, background_scale){
  
  # wrapper function to generate background from the scalar and mu values
  
  x_ind <- as.numeric(as.character(x_ind))
  
  # x_ind is a vector of length 3 for weekday, month, christmas
  
  
  if(x_ind[3] == 1){ # if it is christmas
    res <- theta_mu[20]/background_scale
    
  } else { # multiplication between day and month (not christmas)
    
    res <- prod(theta_mu[x_ind[1]], theta_mu[x_ind[2] + 7])/background_scale
  }
  
  return(res)
  
  
}





# make the indicator matrix

make_my_indicator <- function(x_dates){
  
  ## returns a matrix of Nx3 size
  
  
  # construct weekday vector
  
  weekdays <- as.numeric(factor(weekdays(x_dates), levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"), ordered = TRUE))
  
  # construct month vector
  
  months <- month(x_dates)
  
  # construct Christmas vector
  
  xmas <- ifelse(month(x_dates) == 12 & as.numeric(format(x_dates, format = "%d") ) %in% c(24, 25, 26, 27), 1, 0)
  
  # make model matrix
  
  res <- cbind(weekdays, months, xmas)
  return(res)
}




mu_date_intensity_at <- function(x_dates, theta_mu, background_scale){
  
  # give intensity at certain times
  
  N <- length(x_dates)
  x_ind <- make_my_indicator(x_dates)
  
  res <- numeric(N)
  
  for(i in 1:N){
    res[i] <- mu_date_intensity(x_ind = x_ind[i, ], theta_mu = theta_mu, background_scale = background_scale)
  }
  
  return(res)
}



# evaluate the likelihood for given events and parameters

mu_date_likelihood <- function(theta_mu, ts_dates, firstday, lastday, background_scale , displayOutput = FALSE, 
                               use_c_mu = FALSE, c_mu = NA, ts_articles = NA){
  
  
  N <- length(ts_dates)
  ts_ind <- make_my_indicator(ts_dates)
  
  
  ts_dates_all <- seq(as.Date(firstday), as.Date(lastday), by = "days")
  P <- length(ts_dates_all)
  ts_ind_all <- make_my_indicator(ts_dates_all)
  
  
  
  ## evaluate the likelihood
  
  loglik <- 0
  
  # part 1: sum over events
  
  for(i in 1:N){
    loglik <- loglik + ifelse(use_c_mu, c_mu[ts_articles[i]], 1) * log(mu_date_intensity(x_ind = ts_ind[i, ], theta_mu = theta_mu, 
                                                                                         background_scale = background_scale) )
  }
  
  # part 2: integral of intensity
  
  for(j in 1:P){ 
    loglik <- loglik - ifelse(use_c_mu, sum(c_mu), 1) * mu_date_intensity(x_ind = ts_ind_all[j, ], theta_mu = theta_mu, 
                                                                          background_scale = background_scale)
  }
  
  
  # return
  
  if(displayOutput){
    print(c(theta_mu, -loglik))
  }
  
  return(-loglik)
  
}



# maximise likelihood for given events 
mu_date_maxlik <- function(ts_dates, firstday, lastday, initval = NA, displayOutput=FALSE){
  
  
  if (is.na(initval) ) {initval <- rep(0.5, times = 20)}
  
  temp <- optim(initval, mu_date_likelihood, method='L-BFGS-B',
                lower = rep(0.01, 20), 
                upper = rep(Inf, 20),
                ts_dates = ts_dates, firstday = firstday, lastday = lastday, background_scale = 1,
                displayOutput = displayOutput) 
  temppar <- temp$par
  return(list(theta_mu =  temp$par,
              lik= -temp$value))
}





# construct event data set for maximisation & return scaled estimate



mu_date_estimate <- function(ts, firstday = NA, lastday = NA, estimate_at = NA,
                             resample_dim = TRUE, N_resample = NA, rescale = "none"){
  
  #ts has format: date, product
  
  if(is.na(firstday)){ firstday <- min(ts[, 1])}
  
  if(is.na(lastday)){ lastday <- max(ts[, 1])}
  
  
  
  if(is.na(estimate_at)){
    estimate_at <- seq(as.Date(firstday), as.Date(lastday), by = "days")
  }
  
  N <- length(estimate_at)
  ind_estimate_at <- make_my_indicator(estimate_at)
  
  
  
  ts[, 2] <- as.numeric(as.factor(ts[, 2]))
  M <- max(ts[, 2])
  
  
  
  if(resample_dim){ # every dimension contributes the same number of observations
    
    
    ts1 <- as.Date(x = integer(0), origin = "1970-01-01")
    
    for(m in 1:M){
      
      n_resample <- min(sum(ts[, 2] == m), N_resample)
      
      ts_m <- ts[ts[, 2] == m, 1]
      n_m <- length(ts_m)
      
      ind <- sample(x = 1:n_m, size = n_resample, replace = FALSE)
      ts1 <- c(ts1, ts_m[ind])
    }
    
    ts1 <- sort(ts1)
    
  } else {
    
    ts1 <- ts[, 1]
  }
  
  
  # estimate the MLE
  mle_theta_mu <- mu_date_maxlik(ts_dates = ts1, firstday = firstday, lastday = lastday)$theta_mu
  
  
  
  res <- numeric(N)
  
  for(i in 1:N){
    res[i] <- mu_date_intensity(x_ind = ind_estimate_at[i, ], theta_mu = mle_theta_mu, background_scale = 1)
  }
  
  res <- cbind(as.numeric(estimate_at) - min(as.numeric(estimate_at)) , res)
  
  
  if(rescale == "one"){ # estimations to have mean 1
    
    est_mean <- mean(res[, 2])
    res[, 2] <- res[, 2] * 1/est_mean
    
  } else if(rescale == "M"){
    
    res[, 2] <- res[, 2] / M
    
  } else if(rescale == "pseudoM"){
    
    res[, 2] <- res[, 2] / (length(ts1)/N_resample)
    
  } else if(is.numeric(rescale)){
    
    res[, 2] <- res[, 2] / rescale
  }
  
  
  return(list(estimates = res, mle = mle_theta_mu))
  
}



mu_date_estimate_param <- function(x_dates, theta_mu, background_scale = background_scale){
  
  N <- length(x_dates)
  ind_estimate_at <- make_my_indicator(x_dates)
  
  res <- numeric(N)
  
  for(i in 1:N){
    res[i] <- mu_date_intensity(x_ind = ind_estimate_at[i, ], theta_mu = theta_mu, background_scale = background_scale)
  }
  
  return(res)
  
}




#### compare models ####


BIC_DIC_Hawkes <- function(ts, M, theta_mu, fit, bool = NA, maxTime, firstday, lastday, background_scale, S = 100, background_only = FALSE, capat){
  
  # ts has three colums: eventtime, article, date
  
  N_fit <- length(fit$v[,1])
  
  if(is.na(bool[1])){
    bool <- rep(TRUE, N_fit)
  }
  
  c_mu_mean <- apply(fit$v[bool, ], 2, mean)
  
  K_mean <-  as.matrix(apply(fit$K[bool, , ], c(2,3), mean))
  beta_vec_mean <- apply(fit$beta[bool, ], 2, mean)
  
  
  
  # set up beta 
  beta_mean <- matrix(beta_vec_mean[2], nrow = M, ncol = M)
  diag(beta_mean) <- beta_vec_mean[1]
  
  backgroundest <- mu_date_estimate_param(x_dates = seq(from = as.Date(firstday), to = as.Date(lastday), by = "days"), theta_mu, 
                                          background_scale = background_scale)
  
  params <- list(1:(maxTime-1), backgroundest, c_mu_mean)
  
  
  print("calculate loglik")
  
  
  loglik <- Likelihood_EXP_multi(ts = ts, mu = mu_step_scale_single, K = K_mean, beta = beta_mean, maxTime = maxTime, 
                                 mu_function = TRUE, params = params, approxInt = TRUE, capat = capat)
  
  
  
  
  
  # effective number of param
  
  print("start samples")
  
  loglik_samples <- numeric(S)
  
  for(s in 1:S){
    
    if(s %% 10 == 0) print(s)
    
    ind <- sample((1:N_fit)[bool], 1)
    
    c_mu_iter <- fit$v[ind, ]
    K_iter <- as.matrix(fit$K[ind, , ])
    beta_vec_iter <- fit$beta[ind,  ]
    
    beta_iter <- matrix(beta_vec_iter[2], nrow = M, ncol = M)
    diag(beta_iter) <- beta_vec_iter[1]
    
    params_iter <- list(1:(maxTime-1), backgroundest, c_mu_iter)
    
    
    loglik_samples[s] <- Likelihood_EXP_multi(ts = ts, mu = mu_step_scale_single, K = K_iter, beta = beta_iter, maxTime = maxTime, 
                                              mu_function = TRUE, params = params_iter, approxInt = TRUE, capat = capat)
    
    
    
    
    
  }
  
  p_dic <- 2 * (loglik - mean(unlist(loglik_samples)) ) 
  
  if(!background_only){
    k_bic <- M + M*M + 2 
  } else {
    k_bic <- M
  }
  
  
  bic <- k_bic * log(length(ts[, 1])) + 2*loglik
  
  dic <- - 2*loglik + 2*p_dic
  
  return(list(BIC = as.numeric(bic), 
              DIC = as.numeric(dic), 
              loglik = as.numeric(loglik), 
              k_bic = as.numeric(k_bic), 
              p_dic = as.numeric(p_dic)))
  
}



predlik_Hawkes <- function(ts_test, ts_train, M, theta_mu, fit, bool = NA, starttime_test, maxTime_test, firstday_test, lastday_test, background_scale, capat){
  
  # predictive likelihood
  
  if(is.na(bool[1])){
    bool <- rep(TRUE, dim(fit$v)[1])
  }
  
  c_mu_mean <- apply(fit$v[bool,], 2, mean)
  K_mean <-  as.matrix(apply(fit$K[bool, , ], c(2,3), mean))
  beta_vec_mean <- apply(fit$beta[bool, ], 2, mean)
  
  
  # set up beta 
  beta_mean <- matrix(beta_vec_mean[2], nrow = M, ncol = M)
  diag(beta_mean) <- beta_vec_mean[1]
  
  backgroundest <- mu_date_estimate_param(x_dates = seq(from = as.Date(firstday_test), to = as.Date(lastday_test), by = "days"), theta_mu = theta_mu, 
                                          background_scale = background_scale)
  
  params <- list((starttime_test+1):(maxTime_test-1), backgroundest, c_mu_mean)
  
  res <- Likelihood_EXP_multi(ts = ts_test, mu = mu_step_scale_single, K = K_mean, beta = beta_mean, maxTime = maxTime_test, 
                              mu_function = TRUE, params = params, approxInt = TRUE, capat = capat,
                              useprevious = TRUE, ts_previous = ts_train, starttime_test = starttime_test)
  
  return(res)
}





log_number_of_events <- function(ts){
  return(log(length(ts)))
}


interevent_functions <- function(ts){
  
  interevent <- diff(ts)
  res <- numeric(5)
  
  res[1] <- median(interevent)
  res[2] <- var(interevent)
  res[3] <- median(interevent)/mean(interevent)
  res[4] <- mean(interevent[interevent > quantile(interevent, 0.9)])
  res[5] <- mean(interevent[interevent < median(interevent)])
  
  return(res)
  
}



ripleyK <- function(ts, W){
  
  len <- length(ts)
  mat <- matrix(data = ts, ncol = len, nrow = len, byrow = TRUE)
  matsub <- mat - ts
  count <- sum(matsub > 0 & matsub < W) * 2
  
  return(count/(len*(len-1)))
}

summary_statistics <- function(ts, W){
  
  N_sumstats <- length(W) + 6 
  
  res <- numeric(N_sumstats)
  
  res[1] <- log_number_of_events(ts)
  res[2:6] <- interevent_functions(ts)
  
  for(i in 1:length(W)){
    
    res[6 + i] <- ripleyK(ts, W[i])
    
  }
  
  acfs <- summary_acf(ts)
  length_acf <- length(acfs)
  
  if(length_acf < 5){
    acfs <- c(acfs, rep(0, 5-length_acf))
  }
  
  res <- c(res, acfs)
  
  
  names_summstats<- c("log.no.obs", "median.inter","var.inter", "medmean.inter", "meanq90.inter", "meanq50.inter", 
                      paste0("Ripley", W),
                      paste0("ACF", 1:5))
  
  names(res) <- names_summstats
  
  return(res)
}



summary_acf <- function(ts){
  
  interevent <- diff(ts)
  temp <- acf(interevent, 5, plot = FALSE)
  res <- c(temp$acf)[-1]
  return(res)
}


beta_vec_to_mat <- function(beta_vec, M){
  beta_mat <- matrix(beta_vec[2], M, M)
  diag(beta_mat) <- beta_vec[1]
  return(beta_mat)
}



find_xlim <- function(x, truevalue = NA){
  qq <- quantile(x, c(0.025, 0.975))
  
  if(!is.na(truevalue)){
    qq[1] <- min(qq[1], truevalue)
    qq[2] <- max(qq[2], truevalue)
  }
  
  return(qq)
}








