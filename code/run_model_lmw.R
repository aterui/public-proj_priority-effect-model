
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))


# jags --------------------------------------------------------------------

#set.seed(1)

nsp <- 20
A <- matrix(runif(nsp * nsp, 0, 1),
            nsp,
            nsp)

diag(A) <- 1

df_a <- as_tibble(c(A)) %>% 
  bind_cols(which(!is.na(A), arr.ind = T)) %>% 
  rename(alpha_prime = value) %>% 
  mutate(alpha = alpha_prime * 1.5 / 100)


# data --------------------------------------------------------------------

set.seed(1)
list_dyn <- cdyns::cdynsim(n_species = nsp,
                           n_timestep = 30,
                           r_type = "constant",
                           r = 1,
                           int_type = "manual",
                           alpha = A,
                           k = 100, 
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
  mutate(species = as.numeric(factor(species))) %>% 
  group_by(species) %>% 
  mutate(count0 = lag(count),
         log_r = log(count) - log(count0)) %>% 
  drop_na(log_r) %>% 
  mutate(t = timestep - 1)


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

## parameters ####
para <- c("b0",
          "b",
          "p")

# jags --------------------------------------------------------------------

d_jags <- list(Y = df_jags$log_r,
               X = df_jags$count0,
               Year = df_jags$t,
               Species = df_jags$species,
               Nsample = nrow(df_jags),
               Nsp = n_distinct(df_jags$species),
               Nyear = n_distinct(df_jags$t),
               W = diag(nsp))

## run jags ####
m <- read.jagsfile("code/model_lmw.R")
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
  filter(str_detect(param, "b\\[\\d{1,},\\d{1,}\\]")) %>% 
  mutate(id = str_extract(param, "\\d{1,},\\d{1,}")) %>% 
  separate(id, into = c("row", "col"),
           convert = T) %>% 
  dplyr::select(row,
                col,
                median = `50%`,
                low = `2.5%`,
                high = `97.5%`) %>% 
  left_join(df_a,
            by = c("row", "col"))

df_plot %>% 
  ggplot(aes(x = alpha,
             y = median)) +
  geom_point() +
  geom_abline(intercept = 0,
              slope = 1) +
  theme_bw() +
  theme(panel.grid = element_blank())

