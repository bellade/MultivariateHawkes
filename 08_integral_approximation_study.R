

mus <- 0.25
Ks <- c(-0.4, -0.25, 0, 0.25, 0.5, 0.75)
betas <- c(0.2, 0.35, 0.5)
scenarios <- c("beta_low", "beta_medium", "beta_high")

maxTime <- 100

df <- data.frame(mu = mus, K = rep(Ks, each = 3), beta = betas, scenario = scenarios,
                 Int_true = NA, Int_app = NA, Error_abs = NA, Error_percentage = NA)

N_each <- 25 # was 10


df <- df[rep(seq_len(nrow(df)), each = N_each), ]

N <- dim(df)[1]


set.seed(4)
for(i in 1:N){
  
  print(i)
  
  mu_iter <-df$mu[i]
  K_iter <- as.matrix(df$K[i])
  beta_iter <- as.matrix(df$beta[i])
  
  
  some_negative <- TRUE
  counter <- 0
  
  while(some_negative){
    
    ts <- simulate_EXP_multi(M = 1, mu = mu_iter, K = K_iter, beta = beta_iter, maxTime = maxTime, uppereval = ifelse(K_iter>0.3, 10, 5))
    
    N_ts <- dim(ts)[1]
   
    some_negative <- FALSE
    
    for(j in 1:N_ts){
      
      intensity_at_t <- intensity_EXP_multi(t = ts[j, ] + c(0.00001, 0), ts = ts, 
                                            mu = mu_iter, K = K_iter, beta = beta_iter, capat = NA)
      
      some_negative <- some_negative | (intensity_at_t < 0)
      
      
    }
    
    counter <- counter+1
    if(counter > 500) stop("counter exceeded")
    
  }
  
  
  df$Int_app[i] <- approxIntegral(ts = ts, maxTime = maxTime,
                 mu = mu_iter, K = K_iter, beta = beta_iter, capat = 0)
  
  mu_iter <-df$mu[i]
  K_iter <- df$K[i]
  beta_iter <- df$beta[i]
  
  df$Int_true[i] <- mu_iter * maxTime + K_iter * sum(1 - exp(-beta_iter * (maxTime - ts[, 1])))
  
}





df$Error_abs <- df$Int_app - df$Int_true
df$Error_percentage <- df$Error_abs/df$Int_true

df$scenario <- factor(df$scenario, ordered = TRUE, levels = c("beta_low", "beta_medium", "beta_high"), 
                      labels = c("beta == 0.2", "beta == 0.35", "beta == 0.5"))


colorBlind <- c("#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



p1 <- ggplot(data = df, aes(x = K, y = Error_percentage, group = scenario, col = scenario)) +
  geom_point(alpha = 0.5, size = 2) + 
  scale_color_manual(values = colorBlind) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(text = element_text(size=20)) +
  labs(x = "K", y = "(Approximation - True)/True") 


p2 <- p1 +  facet_wrap(~scenario, 
                 labeller = label_parsed) + 
  scale_y_continuous(labels =scales::percent_format(accuracy = 1), limits = c(-0.03, 0.03)) +
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 4,
    shape = 18
  )


pdf("approximation.pdf", height = 4, width = 12)
print(p2)
dev.off() 

