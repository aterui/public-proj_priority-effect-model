model {
  
  scale0 <- 2.5
  df0 <- 6
  
  # liklihood ---------------------------------------------------------------
  for (i in 1:Nsample) {
    Y[i] ~ dnorm(mu[Year[i], Species[i]], tau)
  }
  
  for (t in 1:Nyear) {
    for (j in 1:Nsp) {
      mu[t, j] <- b[j, 1] + b[j, 2] * scl_x_t[t, j]
    }
  }
  
  # prior -------------------------------------------------------------------
  s <- rep(2.5, K)
  
  for (j in 1:Nsp) {
    b[j, 1:K] ~ dmnorm(mu_b[], TAU[,])
  }
  
  TAU[1:K, 1:K] ~ dscaled.wishart(s, 2)
  SIGMA[1:K, 1:K] <- inverse(TAU[,])
  tau ~ dscaled.gamma(2.5, 3)
  
  for (k in 1:K) {
    mu_b[k] ~ dnorm(0, pow(2.5, -2))
    tau_b[k] <- TAU[k, k]
    sigma_b[k] <- sqrt(1/tau_b[k])
    cv_b[k] <- sigma_b[k] / abs(mu_b[k])
  }
}

data {
  
  for (t in 1:Nyear) {
    for(j in 1:Nsp) {
      scl_x_t[t, j] <- x_t[t, j] / sd(x_t[, j])
      x_t[t, j] <- sum(x[t, ])
    }
  }
  
  for (i in 1:Nsample) {
    x[Year[i], Species[i]] <- X[i]
  }
}