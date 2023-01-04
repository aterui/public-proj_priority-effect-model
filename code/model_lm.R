model {
  
  tau0 <- pow(scale0, -2)
  scale0 <- 5
  df0 <- 4
  
  
  # likelihood --------------------------------------------------------------
  
  for (j in 1:Nsp) {
    for (k in 1:Nsp) {
      for (t in 1:Nyear) {
        y[t, j, k] ~ dnorm(mu[t, j, k], tau[j, k])
        mu[t, j, k] <- r[j, k] + b[j] * (x[t, j, j] + z[j, k] * x[t, j, k])
      }
    }
  }
  
  for (j in 1:Nsp) {
    for (k in (1 + j):Nsp) {
      z[j, k] ~ dbern(0.5)
      z[k, j] <- z[j, k]
    }
  }
  
  # prior -------------------------------------------------------------------
  
  for (j in 1:Nsp) {
    b[j] ~ dnorm(0, tau0)
    
    for (k in 1:Nsp) {
      r[j, k] ~ dnorm(0, tau0)
      tau[j, k] ~ dscaled.gamma(scale0, df0)
    }
  }
  
}

data {
  
  for (i in 1:Nsample) {
    y[Year[i], Sp_i[i], Sp_j[i]] <- Y[i]
    x[Year[i], Sp_i[i], Sp_j[i]] <- X[i]
  }
  
}