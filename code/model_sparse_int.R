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
        (1 + inprod(alpha[i, ], d[t - Q, ]))
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
    q[i] ~ dnorm(0, 100)T(0, )
    for (j in 1:Nsp) {
      # alpha[i, j] ~ dexp(theta[i])
      # alpha0[i, j] <- alpha[i, j] / alpha[i, i]

      alpha[i, j] <- q[i] * alpha0[i, j]
      alpha0[i, j] <- W[i, j] + (1 - W[i, j]) * alpha_prime[i, j]
      alpha_prime[i, j] ~ dnorm(1, tau_alpha[i, j])T(0, )
      tau_alpha[i, j] <- z[i, j] * q1 + (1 - z[i, j]) * q0

      z[i, j] ~ dbern(p[i, j])
      p[i, j] <- W[i, j] * p0[1] + (1 - W[i, j]) * p0[2]
    }
  }  
  
  q0 <- pow(10, 2)
  q1 <- pow(1, -2) # sd = 1.5
  p0[1] <- 1
  p0[2] ~ dbeta(1, 1)
  #theta ~ dgamma(0.1, 0.1)
}
