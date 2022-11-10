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
        log(1 + inprod(alpha[i, ], d[t - Q, ]))
    }
  }
  
  # prior -------------------------------------------------------------------
  
  tau0 <- 1 / 50
  scale0 <- 2.5
  df0 <- 6
  
  ## low-level parameters ####
  log_r0 ~ dnorm(0, tau0)
  
  tau_r ~ dscaled.gamma(scale0, df0)
  sigma_r <- sqrt(1 / tau_r)
  
  for (i in 1:Nsp) {
    log_d[Nyr1, i] ~ dnorm(0, 0.1)
    log_r[i] ~ dnorm(log_r0, tau_r)
    
    tau[i] ~ dscaled.gamma(scale0, df0)
    sigma[i] <- sqrt(1 / tau[i])
    
    tau_obs[i] ~ dscaled.gamma(scale0, df0)
    sigma_obs[i] <- sqrt(1 / tau_obs[i])
  }  
  
  ## interaction
  shape <- 0.5
  p ~ dbeta(shape, shape)
  
  ### select neutral or niche-structured
  for (i in 1:Nsp) {
    q0[i] ~ dnorm(0, 1)T(0,)
    
    for (j in 1:Nsp) {
      alpha[i, j] <- alpha_prime[i, j] * q0[i]
      alpha_prime[i, j] <- W[i, j] + (1 - W[i, j]) * alpha1[i, j]
      alpha1[i, j] ~ dnorm(mu_alpha[i, j], tau_alpha[i, j])T(0, 2)
      
      mu_alpha[i, j] <- z[i, j] * 1 + (1 - z[i, j]) * 0
      tau_alpha[i, j] <- 100#z[i, j] * 100 + (1 - z[i, j]) * 1
    }
  }  
  
  
  for (i in 1:Nsp) {
    
    for (j in (i + 1):Nsp) {
      z[i, j] ~ dbern(p)
      z[j, i] <- z[i, j]
    }
    
    for (i_prime in i) {
      z[i, i_prime] <- 1
    }
    
  }
  
}
