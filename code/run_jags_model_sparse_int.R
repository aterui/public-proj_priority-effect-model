
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))

# data --------------------------------------------------------------------

## parameters
set.seed(1)

nsp <- 25
r0 <- 1.5
k <- 100

A <- matrix(runif(nsp * nsp, 0, 1),
            nsp,
            nsp)

diag(A) <- 1

df_a <- as_tibble(c(A)) %>% 
  bind_cols(which(!is.na(A), arr.ind = T)) %>% 
  rename(alpha0 = value) %>% 
  mutate(alpha = alpha0 * r0 / k)

## simulate
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

## data sort
df0 <- list_dyn$df_dyn %>% 
  mutate(count = rpois(nrow(.), lambda = density))

df0 <- df0 %>% 
  group_by(species) %>% 
  mutate(index = all(count > 0)) %>% 
  filter(index == TRUE) %>% 
  ungroup() %>% 
  mutate(species = as.numeric(factor(species)))


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
m <- read.jagsfile("code/model_sparse_int.R")

## parameters ####
para <- c("p0",
          "log_r",
          "sigma_obs",
          "sigma",
          "alpha",
          "alpha0",
          "z")

# jags --------------------------------------------------------------------

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

df_plot <- mcmc_summary %>% 
  as_tibble(rownames = "param") %>% 
  filter(str_detect(param, "z")) %>% 
  mutate(id = str_extract(param, "\\d{1,},\\d{1,}")) %>% 
  separate(id, into = c("row", "col"),
           convert = T) %>% 
  dplyr::select(row,
                col,
                mean,
                median = `50%`,
                low = `2.5%`,
                high = `97.5%`) %>% 
  left_join(df_a,
            by = c("row", "col"))

df_plot %>% 
  filter(row != col) %>% 
  ggplot(aes(x = alpha0,
             y = mean)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank())

n_distinct(df0$species)