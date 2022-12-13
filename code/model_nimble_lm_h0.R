m0 <- nimbleCode({
  
  sigma0 <- 5
  scale0 <- 2.5
  df0 <- 3
  
  # liklihood ---------------------------------------------------------------
  for (i in 1:Nsample) {
    Y[i] ~ dnorm(mu[Year[i], Species[i]], sd = sigma)
  }
  
  for (j in 1:Nsp) {
    for (t in 1:Nyear) {
      mu[t, j] <- beta[j] - mu_beta[2] * scl_x_t[t]
    }
  }
  
  # prior -------------------------------------------------------------------
  sigma ~ T(dt(mu = 0, sigma = scale0, df = df0), 0, )
  
  for (j in 1:Nsp) {
    beta[j] ~ dnorm(mu_beta[1], sd = sigma_beta)
  }
  
  for (k in 1:K) {
    mu_beta[k] ~ dnorm(0, sd = sigma0)
  }
  sigma_beta ~ T(dt(mu = 0, sigma = scale0, df = df0), 0,)
  
  # data --------------------------------------------------------------------
  
  for (t in 1:Nyear) {
    scl_x_t[t] <- (x_t[t] - mean(x_t[1:Nyear])) / sd(x_t[1:Nyear])
    x_t[t] <- sum(x[t, 1:Nsp])
    for(j in 1:Nsp) {
      x1[t, j] <- x[t, j]
      x2[t, j] <- sum(x[t, 1:Nsp]) - x[t, j]
    }
  }
  
  for (i in 1:Nsample) {
    x[Year[i], Species[i]] <- X[i]
  }
})
