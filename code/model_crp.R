model {
  
  scale0 <- 2.5
  df0 <- 1
  
  # likelihood --------------------------------------------------------------
  for (i in 1:Nsample) {
    Y[i] ~ dnorm(mu[Year[i], Species[i]], tau)
  }

  for (t in 1:Nyear) {
    for (j in 1:Nsp) {
      mu[t, j] <- beta[j, 1] - beta[j, 2] * (x[t, j] + p[j] * x0[t, j])
    }
  }
  
  # prior -------------------------------------------------------------------
  alpha ~ dunif(0, 1)
  phi ~ dgamma(1, 1)
  mu_p <- 1 / (1 + phi)
  b0 ~ dnorm(0, pow(2.5, -2))T(0, )
  
  for (j in 1:Nsp) {
    
    beta[j, 1] ~ dnorm(0, pow(2.5, -2))
    beta[j, 2] ~ dnorm(0, pow(2.5, -2))T(0, )
    p[j] ~ dbeta(1, phi)
    #beta[j, 3] ~ dnorm(0, pow(2.5, -2))T(, beta[j, 2])#<- alpha * beta[j, 2]
    
  }
  
  tau ~ dscaled.gamma(scale0, df0)
  sigma <- sqrt(1 / tau)
}

data {
  for (t in 1:Nyear) {
    for (j in 1:Nsp) {
      x0[t, j] <- sum(x[t, ]) - x[t, j]
    }
  }
  
  for (i in 1:Nsample) {
    x[Year[i], Species[i]] <- X[i]
  }
}