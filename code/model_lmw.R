model {
  
  scale0 <- 2.5
  df0 <- 6
  
  # likelihood --------------------------------------------------------------
  for (i in 1:Nsample) {
    Y[i] ~ dnorm(mu[Year[i], Species[i]], tau)
  }
  
  for (t in 1:Nyear) {
    for (j in 1:Nsp) {
      mu[t, j] <- b0[j] - b1[j] * inprod(w[], x[t, ])
    }
  }
  
  
  # prior -------------------------------------------------------------------
  
  for (j in 1:Nsp) {
    b0[j] ~ dnorm(0, 1)
    b1[j] ~ dnorm(0, 1)T(0, )
    for (k in 1:Nsp) {
      b[j, k] <- b1[j] * w[k] * Nsp
    }
  }
  
  w[1:Nsp] ~ ddirch(alpha[])
  alpha[1:Nsp] <- rep(1, Nsp)
  
  tau ~ dscaled.gamma(scale0, df0)  
}

data {
  for (i in 1:Nsample) {
    x[Year[i], Species[i]] <- X[i]
  }
}