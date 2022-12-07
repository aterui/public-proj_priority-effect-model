model {
  
  # likelihood --------------------------------------------------------------
  
  ## observation
  for (n in 1:Nsample) {
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
        inprod(alpha[i, ], d[t - Q, ])
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
    #q[i] ~ dnorm(0, 1)T(0, )
    q[i] ~ dnorm(-6, 3)
    psi[i] ~ dscaled.gamma(1, 1)
    for (j in 1:Nsp) {
      log(alpha[i, j]) <- log_alpha[i, j]
      log_alpha[i, j] ~ dnorm(q[i], 1 / (psi[i] * phi))
      z[i, j] <- step(alpha[i, j] - alpha[i, i])
      # alpha[i, j] <- q[i] * alpha0[i, j]
      # alpha0[i, j] <- W[i, j] + (1 - W[i, j]) * b[i, j]
      # b[i, j] ~ dnorm(0, 1)T(0, )
      #z[i, j] <- step(b[i, j] - 0.8)
      # b[i, j] <- z[i, j] * b1 + (1 - z[i, j]) * b0
      # z[i, j] ~ dbern(p[i, j])
      # p[i, j] <- W[i, j] * p0[1] + (1 - W[i, j]) * p0[2]
    }
  }  
  
  phi ~ dscaled.gamma(1, 1)
  
  b1 ~ dnorm(1, pow(1, -2))T(0, )
  b0 <- 0
  
  p0[1] <- 1
  p0[2] ~ dbeta(1, 1)
}
