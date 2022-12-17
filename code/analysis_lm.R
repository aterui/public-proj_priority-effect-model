
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))


# parameters --------------------------------------------------------------

#set.seed(1)
nsp <- 20
k <- 100
v_r <- rnorm(nsp, 1, 1)
v_r[v_r < 0.5] <- 0.5
A <- matrix(0,#rexp(nsp * nsp, rate = 1/0.5),
            nsp,
            nsp)

A[1, ] <- 0.5
A[2, ] <- 0.1

diag(A) <- 1


# simulate ----------------------------------------------------------------

df_est <- foreach(i = 1:100, .combine = c) %do% {
  
  print(i)
  
  list_dyn <- cdyns::cdynsim(n_species = nsp,
                             n_timestep = 20,
                             r_type = "constant",
                             r = v_r,
                             int_type = "manual",
                             alpha = A,
                             k = k, 
                             sd_env = 0.1,
                             model = "ricker",
                             immigration = 0)
  
  df0 <- list_dyn$df_dyn %>% 
    mutate(count = rpois(nrow(.), lambda = density))
  
  df_jags <- df0 %>% 
    group_by(species) %>%
    summarize(index = all(count > 0)) %>%
    right_join(df0, by = "species") %>%
    ungroup() %>%
    filter(index == TRUE) %>%
    group_by(species) %>%
    mutate(count0 = lag(count),
           log_r = log(count) - log(count0)) %>%
    drop_na(log_r) %>% 
    mutate(t = timestep - 1) %>% 
    ungroup() %>% 
    group_by(timestep) %>% 
    mutate(x = sum(count0)) %>% 
    ungroup()
  
  #A <- A[unique(df_jags$species), unique(df_jags$species)]
  
  # analysis ----------------------------------------------------------------
  
  df_jags <- df_jags %>% 
    mutate(species = as.numeric(factor(species)),
           spf = factor(species)) %>% 
    group_by(species) %>% 
    mutate(scl_x = scale(x),
           scl_x0 = scale(count0)) %>% 
    ungroup()
  
  fit1 <- lm(log_r ~ spf * x, df_jags)
  fit2 <- lm(log_r ~ spf + x, df_jags)
  
  bf <- AIC(fit1) - AIC(fit2)
  
  return(bf)
}

plot(density(df_est))
