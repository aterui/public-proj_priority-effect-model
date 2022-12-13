m1 <- nimbleCode({
  
  sigma0 <- 5
  scale0 <- 2.5
  df0 <- 3
  
  # liklihood ---------------------------------------------------------------
  for (i in 1:Nsample) {
    Y[i] ~ dnorm(mu[Year[i], Species[i]], sd = sigma)
  }
  
  for (j in 1:Nsp) {
    for (t in 1:Nyear) {
      mu[t, j] <- beta[j, 1] - beta[j, 2] * scl_x_t[t]
    }
  }
  
  # prior -------------------------------------------------------------------
  sigma ~ T(dt(mu = 0, sigma = scale0, df = df0), 0, )
  
  for (j in 1:Nsp) {
    beta[j, 1:K] ~ dmnorm(mu_beta[1:K], cholesky = m_u[1:K, 1:K], prec_param = 0)
  }
  m_u_star[1:K, 1:K] ~ dlkj_corr_cholesky(2, K)
  m_u[1:K, 1:K] <- uppertri_mult_diag(m_u_star[1:K, 1:K], VSD[1:K])
  
  for (k in 1:K) {
    mu_beta[k] ~ dnorm(0, sd = sigma0)
  }
  
  
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


uppertri_mult_diag <- nimbleFunction(
  run = function(mat = double(2), vec = double(1)) {
    returnType(double(2))
    p <- length(vec)
    out <- matrix(nrow = p, ncol = p, init = FALSE)
    for(i in 1:p)
      out[ , i] <- mat[ , i] * vec[i]
    return(out)
  })