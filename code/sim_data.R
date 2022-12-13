
# setup -------------------------------------------------------------------

source(here::here("code/library.R"))


# interaction matrix ------------------------------------------------------
set.seed(1)
nsp <- 10
A <- matrix(runif(nsp * nsp, 0, 0),
            nsp,
            nsp)

diag(A) <- 1

# data --------------------------------------------------------------------
set.seed(1)
r_min <- 1.5
r_max <- 1.5

list_dyn <- cdyns::cdynsim(n_species = nsp,
                           n_timestep = 30,
                           r_type = "constant",
                           r = runif(nsp, r_min, r_max),
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
