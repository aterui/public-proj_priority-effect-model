
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))


# jags --------------------------------------------------------------------

#set.seed(1)
nsp <- 10
r_min <- 0.5
r_max <- 2.5
k <- 100
A <- matrix(runif(nsp * nsp, 0, 0.5),
            nsp,
            nsp)

diag(A) <- 1

# data --------------------------------------------------------------------

#set.seed(1)
list_dyn <- cdyns::cdynsim(n_species = nsp,
                           n_timestep = 30,
                           r_type = "constant",
                           r = runif(nsp, r_min, r_max),
                           int_type = "manual",
                           alpha = A,
                           k = k, 
                           sd_env = 0.1,
                           model = "ricker")

df0 <- list_dyn$df_dyn %>% 
  mutate(count = rpois(nrow(.), lambda = density + 5))

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
  group_by(t) %>% 
  mutate(n_t = sum(count0)) %>% 
  ungroup()

A <- A[unique(df_jags$species), unique(df_jags$species)]

df_jags <- df_jags %>% 
  mutate(species = as.numeric(factor(species)),
         f = factor(species))

# common setup ------------------------------------------------------------

## mcmc setup ####
n_ad <- 100
n_iter <- 5.0E+4
n_thin <- max(3, ceiling(n_iter / 1000))
n_burn <- ceiling(max(10, n_iter/2))
n_chain <- 4
n_sample <- ceiling(n_iter / n_thin)

inits <- replicate(n_chain,
                   list(.RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA),
                   simplify = FALSE)

for (j in 1:n_chain) inits[[j]]$.RNG.seed <- (j - 1) * 10 + 1

## parameters ####
para <- c("mu_b0",
          "mu_b1",
          "b0",
          "b1",
          "sigma",
          "sigma_b0",
          "sigma_b1",
          "lh")

# jags --------------------------------------------------------------------

list_waic <- foreach(Z = c(0, 1)) %do% {
  
  d_jags <- list(Y = df_jags$log_r,
                 X = df_jags$count0,
                 Year = df_jags$t,
                 Species = df_jags$species,
                 Nsample = nrow(df_jags),
                 Nsp = n_distinct(df_jags$species),
                 Nyear = n_distinct(df_jags$t),
                 K = 2,
                 Z = Z)
  
  m <- read.jagsfile("code/model_lm.R")
  post <- suppressMessages(run.jags(m$model,
                                    monitor = para,
                                    data = d_jags,
                                    n.chains = n_chain,
                                    inits = inits,
                                    method = "parallel",
                                    burnin = n_burn,
                                    sample = n_sample,
                                    adapt = n_ad,
                                    thin = n_thin,
                                    n.sims = 4,
                                    module = "glm"))
  
  rhat <- max(MCMCvis::MCMCsummary(post$mcmc)[,"Rhat"],
              na.rm = T)
  print(rhat)
  
  # while(rhat >= 1.1) {
  #   post <- extend.jags(post,
  #                       combine = TRUE,
  #                       sample = n_sample,
  #                       thin = n_thin)
  #   
  #   rhat <- max(MCMCvis::MCMCsummary(post$mcmc)[,"Rhat"],
  #               na.rm = T)
  #   print(rhat)
  # }
  
  MCMCvis::MCMCchains(post$mcmc) %>% 
    as_tibble() %>% 
    dplyr::select(starts_with("lh")) %>% 
    data.matrix() %>% 
    loo::waic() %>% 
    return()
  
}

loo_compare(list_waic[[1]], list_waic[[2]])

d_jags$Nsp