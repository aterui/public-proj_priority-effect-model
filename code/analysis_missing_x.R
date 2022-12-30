
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))

n_sample <- 100
n_var <- 3
mu <- 0 #runif(1, 0, 10)
sigma <- 1 #runif(1, 0, 10)

df0 <- foreach(i = 1:100, .combine = bind_rows) %do% {
  
  # test data ---------------------------------------------------------------
  X <- sapply(1:n_var,
              function(x) rnorm(n = n_sample,
                                mean = 0,
                                sd = 1))
  
  b <- runif(n_var + 1, 0.2, 1)
  b[3] <- -0.2
  
  y <- matrix(c(rep(1, n_sample), X), ncol = n_var + 1) %*% b %>% 
    rnorm(n = n_sample, mean = ., sd = 0.1)
  
  
  # lm test -----------------------------------------------------------------
  
  fit1 <- lm(y ~ X[, 1])
  fit2 <- lm(y ~ X[, 1] + X[, 2] + X[, 3])
  
  return(tibble(b1 = coef(fit1)[1],
                b2 = coef(fit2)[1],
                b_true = b[1]))  
}


g1 <- df0 %>% 
  ggplot(aes(x = b1, y = b_true)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

g2 <- df0 %>% 
  ggplot(aes(x = b2, y = b_true)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

g1 + g2