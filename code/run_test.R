
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))

# common setup ------------------------------------------------------------

## mcmc setup ####
n_ad <- 100
n_iter <- 5.0E+3
n_thin <- max(3, ceiling(n_iter / 250))
n_burn <- ceiling(max(10, n_iter/2))
n_chain <- 4
n_sample <- ceiling(n_iter / n_thin)

inits <- replicate(n_chain,
                   list(.RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA),
                   simplify = FALSE)

for (j in 1:n_chain) inits[[j]]$.RNG.seed <- (j - 1) * 10 + 1

## model file ####
m <- read.jagsfile("code/model_sparse_int.R")

## parameters ####
para <- c("loglik",
          "p0",
          "log_r",
          "sigma_time",
          "sigma_obs",
          "sigma",
          "alpha",
          "rho")

# jags --------------------------------------------------------------------

## psi = 0: niche-structured
## psi = 1: neutral

psi <- c(0, 1)
list_re <- foreach(j = 1:2) %do% {
  
  # data --------------------------------------------------------------------
  set.seed(1)
  list_dyn <- cdyns::cdynsim(n_species = 10,
                             n_timestep = 50,
                             r_type = "constant",
                             r = log(5),
                             int_type = "constant",
                             alpha = 1,
                             k = 100, 
                             sd_env = 0.1,
                             model = "bh")
  
  df0 <- list_dyn$df_dyn %>% 
    mutate(density_obs = density * exp(rnorm(n = nrow(.), mean = 0, sd = 0.01)),
           count = rpois(nrow(.), lambda = density_obs))
  
  ## data for jags ####
  d_jags <- list(N = df0$count,
                 Year = df0$timestep,
                 Species = df0$species,
                 Nsample = nrow(df0),
                 Nyr1 = 1,
                 Nyr = max(df0$timestep),
                 Nsp = n_distinct(df0$species),
                 Nf = 2,
                 W = diag(n_distinct(df0$species)),
                 Q = 1,
                 Psi = psi[j])
  
  ## run jags ####
  post <- run.jags(m$model,
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
                   module = "glm")
  
  mcmc_summary <- MCMCvis::MCMCsummary(post$mcmc)
  print(max(mcmc_summary$Rhat, na.rm = T))

  
  # model evaluation --------------------------------------------------------
  waic_bar <- MCMCvis::MCMCchains(post$mcmc) %>% 
    as_tibble() %>% 
    dplyr::select(starts_with("loglik")) %>% 
    data.matrix() %>% 
    waic()
  
  return(waic_bar)
}


loo::loo_compare(list_re[[1]], list_re[[2]])
