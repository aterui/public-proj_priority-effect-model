
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))

# common setup ------------------------------------------------------------

## mcmc setup ####
n_ad <- 100
n_iter <- 2.0E+3
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
m <- read.jagsfile("code/model_var_int.R")

## parameters ####
para <- c("log_r",
          "sigma_r",
          "sigma_alpha",
          "sigma_obs",
          "alpha")

# jags --------------------------------------------------------------------

set.seed(1)
list_dyn <- cdyns::cdynsim(n_species = 10,
                           n_timestep = 20,
                           r_type = "constant",
                           r = log(10),
                           int_type = "constant",
                           alpha = 1.2,
                           k = 100, 
                           sd_env = 0.1,
                           model = "bh")

df0 <- list_dyn$df_dyn %>% 
  mutate(count = rpois(nrow(.), lambda = density))

df0 %>% 
  group_by(species) %>% 
  summarize(n_obs = sum(count > 0))

df0 %>%
  ggplot(aes(x = timestep,
             y = count,
             color = factor(species))) +
  geom_line() +
  geom_point()

## data for jags ####
d_jags <- list(N = df0$count,
               Year = df0$timestep,
               Species = df0$species,
               Nsample = nrow(df0),
               Nyr1 = 1,
               Nyr = max(df0$timestep),
               Nsp = n_distinct(df0$species),
               W = diag(n_distinct(df0$species)),
               Q = 1)

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
