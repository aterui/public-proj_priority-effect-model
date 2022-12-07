
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))

# common setup ------------------------------------------------------------

## mcmc setup ####
n_ad <- 100
n_iter <- 1.0E+4
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
m <- read.jagsfile("code/model_br.R")

## parameters ####
para <- c("log_r",
          "sigma_obs",
          "sigma",
          "mu_log_br",
          "b")

# jags --------------------------------------------------------------------

set.seed(1)

nsp <- 20
r0 <- 1.5
k <- 100

A <- matrix(runif(nsp * nsp, 0.8, 1),
            nsp,
            nsp)

diag(A) <- 1

df_a <- as_tibble(c(A)) %>% 
  bind_cols(which(!is.na(A), arr.ind = T)) %>% 
  rename(alpha_prime = value) %>% 
  mutate(alpha = alpha_prime * r0 / k)

# data --------------------------------------------------------------------

set.seed(1)
list_dyn <- cdyns::cdynsim(n_species = nsp,
                           n_timestep = 20,
                           r_type = "constant",
                           r = r0,
                           int_type = "manual",
                           alpha = A,
                           k = k, 
                           sd_env = 0.1,
                           model = "ricker")

df0 <- list_dyn$df_dyn %>% 
  mutate(count = rpois(nrow(.), lambda = density))

df0 <- df0 %>% 
  group_by(species) %>% 
  summarize(index = all(count > 0)) %>% 
  right_join(df0, by = "species") %>% 
  ungroup() %>% 
  filter(index == TRUE) %>% 
  mutate(species = as.numeric(factor(species))) %>% 
  group_by(species) %>% 
  mutate(count0 = lag(count),
         log_r = log(count) - log(count0)) %>% 
  drop_na(log_r) %>% 
  mutate(t = timestep - 1)

df0 %>% 
  ggplot(aes(x = count0,
             y = log_r)) +
  geom_point() +
  #geom_abline(intercept = r0, slope = -r0 / k) +
  facet_wrap(facets = ~species)

## data for jags ####
d_jags <- list(N = df0$count,
               Year = df0$timestep,
               Species = df0$species,
               Nsample = nrow(df0),
               Nyr1 = 1,
               Nyr = max(df0$timestep),
               Nsp = n_distinct(df0$species),
               Q = 1)

## run jags ####
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

(mcmc_summary <- MCMCvis::MCMCsummary(post$mcmc))
print(max(mcmc_summary$Rhat, na.rm = T))
