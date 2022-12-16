model {
  sd0 <- 5
  tau0 <- pow(sd0, -2)
  df0 <- 3
  v_s0 <- rep(sd0, K)  
  
  # likelihood --------------------------------------------------------------
  for (i in 1:Nsample) {
    Y[i] ~ dpois(lambda[Year[i], Species[i]])
  }
  
  for (t in 1:Nt) {
    for (j in 1:Nsp) {
      lambda[t, j] <- n[t, j]
      log(n[t, j]) <- x[t, j]
    }
    sum_n[t] <- sum(n[t, ])
  }
  
  for (j in 1:Nsp) {
    for (t in 1:(Nt - 1)) {
      x[t + 1, j] ~ dnorm(x_prime[t, j], tau_r[j])
      x_prime[t, j] <- x[t, j] + b[j, 1] + b[j, 2] * sum_n[t]# + b[j, 3] * sum_n[t, j]
    }
  }
  
  
  # prior -------------------------------------------------------------------
  for (k in 1:K) {
    mu_b[k] ~ dnorm(0, tau0)
    for (j in 1:Nsp) {
      b0[j, k] <- sd(sum_n[]) * b[j, k]
    }
  }
  
  for (j in 1:Nsp) {
    b[j, 1:K] ~ dmnorm(mu_b[], TAU[,])
    tau_r[j] ~ dscaled.gamma(sd0, df0)
    x[1, j] ~ dnorm(0, tau0)
  }  
  
  TAU[1:K, 1:K] ~ dscaled.wishart(v_s0[], K)
  SIGMA[1:K, 1:K] <- inverse(TAU[,])  
}
