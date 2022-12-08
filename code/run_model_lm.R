
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))


# jags --------------------------------------------------------------------

#set.seed(1)
nsp <- 5
r0 <- runif(nsp, 0.5, 2.5)
k <- 100
A <- matrix(runif(nsp * nsp, 0, 0),
            nsp,
            nsp)

diag(A) <- 1

df_a <- as_tibble(c(A)) %>% 
  bind_cols(which(!is.na(A), arr.ind = T)) %>% 
  rename(alpha_prime = value) %>% 
  mutate(alpha = alpha_prime * r0 / k)


# data --------------------------------------------------------------------

#set.seed(1)
list_dyn <- cdyns::cdynsim(n_species = nsp,
                           n_timestep = 30,
                           r_type = "constant",
                           r = r0,
                           int_type = "manual",
                           alpha = A,
                           k = k, 
                           sd_env = 0.1,
                           model = "ricker")

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
n_iter <- 1.0E+4
n_thin <- max(3, ceiling(n_iter / 250))
n_burn <- ceiling(max(10, n_iter/2))
n_chain <- 4
n_sample <- ceiling(n_iter / n_thin)

inits <- replicate(n_chain,
                   list(.RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA,
                        mu_b = c(1, 1)),
                   simplify = FALSE)

for (j in 1:n_chain) inits[[j]]$.RNG.seed <- (j - 1) * 10 + 1

## parameters ####
para <- c("mu_b",
          "sigma_b",
          "cv_b",
          "tau",
          "SIGMA")

# jags --------------------------------------------------------------------

d_jags <- list(Y = df_jags$log_r,
               X = df_jags$count0,
               Year = df_jags$t,
               Species = df_jags$species,
               Nsample = nrow(df_jags),
               Nsp = n_distinct(df_jags$species),
               Nyear = n_distinct(df_jags$t),
               K = 2)

## run jags ####
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

(mcmc_summary <- MCMCvis::MCMCsummary(post$mcmc))
print(max(mcmc_summary$Rhat, na.rm = T))

# plot --------------------------------------------------------------------

df_jags %>%
  ggplot(aes(x = n_t,
             y = log_r,
             color = factor(species))) +
  geom_point() +
  facet_wrap(facets = ~species) +
  geom_smooth(method = "lm",
              se = F)

print(A)
mcmc_summary %>%
  as_tibble(rownames = "param") %>%
  filter(str_detect(param, "cv"))
