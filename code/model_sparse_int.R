model {
  
  # likelihood --------------------------------------------------------------
  
  ## observation
  for (n in 1:Nsample) {
    loglik[n] <- logdensity.pois(N[n], lambda[Year[n], Species[n]])
    N[n] ~ dpois(lambda[Year[n], Species[n]])
  }
  
  for (t in Nyr1:Nyr) {
    for (i in 1:Nsp) {
      lambda[t, i] <- d[t, i]
      #log(d_obs[t, i]) <- log_d_obs[t, i]
      #log_d_obs[t, i] ~ dnorm(log_d[t, i], tau_obs[i])
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
        inprod(alpha[i, ], d[t - Q, ]) +
        inprod(xi[ , i], eta[t - Q, ])
    }
  }
  
  # prior -------------------------------------------------------------------
  
  tau0 <- 1 / 50
  scale0 <- 2.5
  df0 <- 6
  a1 <- 2
  a2 <- 3
  
  ## low-level parameters ####
  log_r0 ~ dnorm(0, tau0)
  
  for (i in 1:Nsp) {
    log_d[Nyr1, i] ~ dt(-2, 0.5, df0)
    log_r[i] <- Psi * log_r0 + (1 - Psi) * log_r_prime[i]
    log_r_prime[i] ~ dnorm(0, tau0)
    
    tau_obs[i] ~ dscaled.gamma(scale0, df0)
    sigma_obs[i] <- sqrt(1 / tau_obs[i])
  }  
  
  ## sparse prior for alpha[i, j] for i != j ####
  alpha0 ~ dnorm(0, 1)
  
  ### select neutral or niche-structured
  for(i in 1:Nsp) {
    for(j in 1:Nsp) {
      alpha[i, j] <- Psi * alpha0 + (1 - Psi) * alpha_prime[i, j] 
    }
  }
  
  for (i in 1:Nsp) {
    for (j in 1:Nsp) {
      alpha_prime[i, j] ~ dnorm(0, tau_alpha[i, j])T(0,)
      tau_alpha[i, j] <- (1 - z[i, j]) * q0 + z[i, j] * q1
      z[i, j] ~ dbern(p[i, j])
      p[i, j] <- W[i, j] * p0[1] + (1 - W[i, j]) * p0[2]
    }
  }  
  
  p0[1] ~ dbeta(1, 1)
  p0[2] ~ dbeta(1, 1)  
  q0 <- 100
  q1 <- 1
  
  ## factor analysis for epsilon ###
  ### latent variables
  for(t in Nyr1:(Nyr - Q)) {
    for (k in 1:Nf) {
      eta[t, k] ~ dnorm(0, 1)
    }
  }
  
  for (k in 1:Nf) {
    for(i in 1:Nsp) {
      xi[k, i] ~ dnorm(0, phi[k, i] * theta[k])
      phi[k, i] ~ dgamma(1.5, 1.5)
    }
  }
  
  ### multiplicative gamma prior for factor loading
  delta[1] ~ dgamma(a1, 1)
  theta[1] <- delta[1]
  
  for(k in 2:Nf) {
    delta[k] ~ dgamma(a2, 1)
    theta[k] <- theta[k - 1] * delta[k]
  }
  
  ### var-covar matrix
  OMEGA[1:Nsp, 1:Nsp] <- t(xi[ , ]) %*% xi[ , ] + diag_omega[ , ]
  TAU[1:Nsp, 1:Nsp] <- inverse(OMEGA[ , ])
  
  for(i in 1:Nsp) {
    #### unique error for species-level time series
    tau[i] ~ dscaled.gamma(scale0, df0)
    sigma[i] <- sqrt(1 / tau[i])
    sigma_time[i] <- sqrt(OMEGA[i, i])
    
    for(j in 1:Nsp) {
      diag_omega[i, j] <- W[i, j] * pow(sigma[j], 2)
      rho[i, j] <- OMEGA[i, j] / sqrt(OMEGA[i, i] * OMEGA[j, j])
    }
  }

  
  # posterior predictive checking -------------------------------------------
  
  for (n in 1:Nsample) {
    y[n] ~ dpois(lambda[Year[n], Species[n]])
    chi_y[n] <- (y[n] - lambda[Year[n], Species[n]])^2 / lambda[Year[n], Species[n]]
    chi_n[n] <- (N[n] - lambda[Year[n], Species[n]])^2 / lambda[Year[n], Species[n]]
  }
  
  sum_chi_y <- sum(chi_y[])
  sum_chi_n <- sum(chi_n[])
  
  bp <- step(sum_chi_y - sum_chi_n)
}
