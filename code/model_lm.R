model {
  
  scale0 <- 2.5
  df0 <- 6
  
  # liklihood ---------------------------------------------------------------
  for (i in 1:Nsample) {
    lh[i] <- log(dnorm(Y[i], mu[Year[i], Species[i]], tau))
    Y[i] ~ dnorm(mu[Year[i], Species[i]], tau)
  }
  
  for (j in 1:Nsp) {
    for (t in 1:Nyear) {
      mu[t, j] <- Z * mu1[t, j] + (1 - Z) * mu0[t, j]
      mu0[t, j] <- b0[j, 1] - b0[j, 2] * x_t[t, j]
      mu1[t, j] <- b1[j, 1] - b1[j, 2] * x_t[t, j]
    }
  }
  
  # prior -------------------------------------------------------------------
  s <- rep(2.5, K)
  
  for (j in 1:Nsp) {
    b0[j, 1] ~ dnorm(mu_b0[1], tau_b0[1])
    b0[j, 2] <- mu_b0[2]
    
    b1[j, 1:K] ~ dmnorm(mu_b1[], TAU[,])
  }
  
  TAU[1:K, 1:K] ~ dscaled.wishart(s, 2)
  SIGMA[1:K, 1:K] <- inverse(TAU[,])
  tau ~ dscaled.gamma(2.5, 3)
  
  for (k in 1:K) {
    mu_b0[k] ~ dnorm(0, pow(2.5, -2))T(0,)
    tau_b0[k] ~ dscaled.gamma(2.5, 3)
    
    mu_b1[k] ~ dnorm(0, pow(2.5, -2))T(0,)
    tau_b1[k] <- TAU[k, k]
    sigma_b1[k] <- sqrt(1 / tau_b1[k])
  }
}

data {
  
  for (t in 1:Nyear) {
    for(j in 1:Nsp) {
      x_t[t, j] <- sum(x[t, ])
      x1[t, j] <- x[t, j]
      x2[t, j] <- sum(x[t, ]) - x[t, j]
    }
  }
  
  for (i in 1:Nsample) {
    x[Year[i], Species[i]] <- X[i]
  }
}