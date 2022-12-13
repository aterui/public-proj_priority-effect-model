model {
  
  scale0 <- 2.5
  df0 <- 6
  
  # liklihood ---------------------------------------------------------------
  for (i in 1:Nsample) {
    lh[i] <- log(dnorm(Y[i], mu[Year[i], Species[i]], tau))
    Y[i] ~ dnorm(mu[Year[i], Species[i]], tau)
  }
  
  for (t in 1:Nyear) {
    for (j in 1:Nsp) {
      mu[t, j] <- Z * mu1[t, j] + (1 - Z) * mu0[t, j]
      mu0[t, j] <- b0[j, 1] + b0[j, 2] * scl_x_t[t]
      mu1[t, j] <- b1[j, 1] + b1[j, 2] * scl_x_t[t]
    }
  }
  
  # prior -------------------------------------------------------------------
  s <- rep(2.5, K)
  
  for (j in 1:Nsp) {
    # model H0
    b0[j, 1] ~ dnorm(mu_b0[1], tau_b0[1])
    b0[j, 2] <- mu_b0[2]
    
    # model H1
    b1[j, 1:K] ~ dmnorm(mu_b1[], TAU[,])
  }
  
  TAU[1:K, 1:K] ~ dscaled.wishart(s, 2)
  SIGMA[1:K, 1:K] <- inverse(TAU[,])
  tau ~ dscaled.gamma(scale0, df0)
  sigma <- sqrt(1 / tau)
  
  for (k in 1:K) {
    # model H0
    mu_b0[k] ~ dnorm(0, pow(scale0, -2))
    tau_b0[k] ~ dscaled.gamma(scale0, df0)
    sigma_b0[k] <- sqrt(1 / tau_b0[k])
    
    # model H1
    mu_b1[k] ~ dnorm(0, pow(scale0, -2))
    tau_b1[k] <- TAU[k, k]
    sigma_b1[k] <- sqrt(tau_b1[k])
  }
}

data {
  
  for (t in 1:Nyear) {
    scl_x_t[t] <- (x_t[t] - mean(x_t[])) / sd(x_t[])
    x_t[t] <- sum(x[t, ])
  }
  
  for (i in 1:Nsample) {
    x[Year[i], Species[i]] <- X[i]
  }
}