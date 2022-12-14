model {
  
  scale0 <- 2.5
  df0 <- 6
  
  # liklihood ---------------------------------------------------------------
  for (i in 1:Nsample) {
    Y[i] ~ dnorm(mu[Year[i], Species[i]], tau)
  }
  
  for (t in 1:Nyear) {
    for (j in 1:Nsp) {
      mu[t, j] <- beta[j, 1] - (beta[j, 2] * x1[t, j] + beta[j, 3] * x2[t, j])
    }
  }
  
  # prior -------------------------------------------------------------------
  s <- rep(2.5, K)
  alpha ~ dbeta(1, 1)
  p ~ dbeta(1, 1)
  
  for (j in 1:Nsp) {
    
    beta[j, 1] ~ dnorm(0, pow(2.5, -2))
    beta[j, 2] ~ dnorm(0, pow(2.5, -2))T(0, )
    beta[j, 3] <- alpha * beta[j, 2]
    
    for (k in 1:Nsp) {
      z[j, k] ~ dbern(p)
    }
        
    for (t in 1:Nyear) {
      x1[t, j] <- sum(z[j, ] * x[t, ])
      x2[t, j] <- sum((1 - z[j, ]) * x[t, ])
    }
    
  }
  
  tau ~ dscaled.gamma(scale0, df0)
  sigma <- sqrt(1 / tau)
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