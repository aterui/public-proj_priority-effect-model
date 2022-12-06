model {
  
  scale0 <- 2.5
  df0 <- 6
  
  # liklihood ---------------------------------------------------------------
  for (i in 1:Nsample) {
    Y[i] ~ dnorm(mu[Year[i], Species[i]], tau)
  }
  
  for (t in 1:Nyear) {
    for (j in 1:Nsp) {
      mu[t, j] <- b0[j] - inprod(b[j, ], x[t, ])
    }
  }
  
  
  # prior -------------------------------------------------------------------
  
  for (j in 1:Nsp) {
    b0[j] ~ dnorm(0, 1)
    q[j] ~ dnorm(0, 1)T(0, )
    for (k in 1:Nsp) {
      b[j, k] <- beta[j, k] * q[j]
      beta[j, k] <- W[j, k] + (1 - W[j, k]) * beta_prime[j, k]
      beta_prime[j, k] <- z[j, k] * beta1[j, k] + (1 - z[j, k]) * beta0[j, k]

      beta0[j, k] ~ dnorm(0, 1)T(0, )
      beta1[j, k] <- 1
      z[j, k] <- W[j, k] * p[1] + (1 - W[j, k]) * p[2]
    }
  }
  p[1] <- 1
  p[2] ~ dbeta(1, 1)
  tau ~ dscaled.gamma(scale0, df0)  
}

data {
  for (i in 1:Nsample) {
    x[Year[i], Species[i]] <- X[i]
  }
}