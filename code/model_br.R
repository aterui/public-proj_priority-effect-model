model {
  
  # likelihood --------------------------------------------------------------
  
  ## observation
  for (n in 1:Nsample) {
    loglik[n] <- logdensity.pois(N[n], lambda[Year[n], Species[n]])
    N[n] ~ dpois(lambda[Year[n], Species[n]])
  }
  
  for (t in Nyr1:Nyr) {
    for (i in 1:Nsp) {
      lambda[t, i] <- d_obs[t, i]
      log(d_obs[t, i]) <- log_d_obs[t, i]
      log_d_obs[t, i] ~ dnorm(log_d[t, i], tau_obs[i])
      log(d[t, i]) <- log_d[t, i]
    }
  }  
  
  ## state
  for (t in (Nyr1 + Q):Nyr) {
    for(i in 1:Nsp) {
      log_d[t, i] ~ dnorm(log_mu_d[t, i], tau[i])
      
      log_mu_d[t, i] <- 
        log_d[t - Q, i] + 
        log_r[i] - 
        x[t - Q, i]
      
      x[t - Q, i] <-
        b[1, i] * d[t - Q, i] + 
        b[2, i] * (sum(d[t - Q, ]) - d[t - Q, i])
    }
  }
  
  # prior -------------------------------------------------------------------
  
  tau0 <- 1 / 2.5
  scale0 <- 2.5
  df0 <- 6
  
  ## low-level parameters ####
  mu_log_r ~ dnorm(0, tau0)
  tau_log_r ~ dscaled.gamma(scale0, df0)
  
  for (i in 1:Nsp) {
    log_d[Nyr1, i] ~ dnorm(0, 0.1)
    log_r[i] ~ dnorm(mu_log_r, tau_log_r)
    
    tau[i] ~ dscaled.gamma(scale0, df0)
    sigma[i] <- sqrt(1 / tau[i])
    
    tau_obs[i] ~ dscaled.gamma(scale0, df0)
    sigma_obs[i] <- sqrt(1 / tau_obs[i])
  }  
  
  ## sparse prior for alpha[i, j] for i != j ####
  for (i in 1:Nsp) {
    for (k in 1:2) {
      b[k, i] ~ dnorm(0, 1)T(0, )
    }
    log_br[i] <- log(b[2, i]) - log(b[1, i])
  }
  
  mu_log_br <- mean(log_br[])
}
