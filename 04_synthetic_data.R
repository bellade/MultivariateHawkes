


#### prior on K star, mu constant ####


M <- 3
capat <- 0
mu <- rep(0.15, times = M)

Kstar <- matrix(0, M, M, byrow = TRUE)

Kstar[1, 2] <- -0.3
Kstar[2, 2] <- 0.3
Kstar[2, 3] <- 0.3
Kstar[3, 2] <- -0.3
Kstar[1, 1] <- 0.3


K <- Kstar_to_K(Kstar)
eigen(apply(K, c(1, 2), max, 0))

beta <- matrix(0.5, M, M)

maxTime <- 500



set.seed(4)
ts <- simulate_EXP_multi(M = M, mu = mu, K = K, beta = beta, maxTime = maxTime, uppereval = 5, capat = capat)
dim(ts)
table(ts[, 2])
par(mfrow=c(1,1))
plot(ts)


params <- list(mu = mu, K = K, Kstar = Kstar, beta = beta, maxTime = maxTime, capat = capat)




set.seed(4)
fit_normal <- Stan_fit_multi_Model1_EXP_exin_approxInt_capat_Kstar(ts = ts, M = M, maxTime = maxTime,
                                                                   prior_K = c(0, 1), prior_beta = c(0, 1), capat = 0.01,
                                                                   iter = 500)


bool_normal <- rep(FALSE, 750)


for(i in 1:750){
  K_temp <- fit_normal$K[i, , ]
  Kpos <- apply(K_temp, c(1, 2), max, 0)
  if(max(abs(eigen(Kpos)$values)) < 1){
    bool_normal[i] <- TRUE
  }
  
 
}


fit <- fit_regHS
bool <- bool_horse

par(mfrow=c(1,M))
for(i in 1:M){
  plot(density(fit_normal$mu[bool_normal, i]), main = paste0("mu ", i), xlim = c(0, 0.6), col = "blue")
  abline(v = mu[i])
  
  
}

par(mfrow=c(M, M))
for(i in 1:M){
  for(j in 1:M){
    plot(density(fit_normal$Kstar[bool_normal, i, j]), main = paste0("Kstar ", i, j), xlim = c(-1, 1), col = "blue")
    abline(v = Kstar[i,j])
    
  }
  
}

par(mfrow=c(1,2))
for(i in 1:2){
  plot(density(fit_normal$beta[bool_normal, i]), main = paste0("beta", i), xlim = c(0, 1), col = "blue")
  abline(v = beta[1, i])
  
}





