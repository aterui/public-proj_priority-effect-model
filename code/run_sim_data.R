
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))


# parameters --------------------------------------------------------------

#set.seed(2)
nt <- 50
nsp <- 2
nn <- 2
r_min <- 0.5
r_max <- 2.5
k <- 100
A <- matrix(runif(nsp^2, 0, 0.25),
            nsp,
            nsp)

A[1:nn, ] <- A[, 1:nn] <- 0.01
A[1:nn, 1:nn] <- 1

diag(A) <- 1


# run simulations ---------------------------------------------------------

#set.seed(1)
list_dyn <- cdyns::cdynsim(n_species = nsp,
                           n_timestep = nt,
                           r_type = "constant",
                           r = c(rep(0.5, nn), runif(nsp - nn, r_min, r_max)),
                           int_type = "manual",
                           alpha = A,
                           k = k, 
                           sd_env = 0,
                           n_warmup = 0,
                           n_burnin = 0,
                           seed = 100,
                           seed_interval = 0,
                           model = "ricker")

df0 <- list_dyn$df_dyn %>% 
  mutate(x = rpois(nrow(.), lambda = density))

df_r <- df0 %>% 
  group_by(species) %>% 
  summarize(index = all(x > 0)) %>% 
  right_join(df0, by = "species") %>% 
  ungroup() %>% 
  filter(index == TRUE) %>% 
  group_by(species) %>% 
  mutate(x0 = lag(x),
         log_r = log(x) - log(x0)) %>% 
  ungroup() %>% 
  drop_na(log_r) %>%
  dplyr::select(species, 
                timestep,
                x,
                x0,
                log_r)

df_jags <- tibble(sp_i = 1:nsp, sp_j = 1:nsp) %>% 
  expand(sp_i, sp_j) %>% 
  left_join(df_r,
            by = c("sp_i" = "species")) %>% 
  left_join(df_r,
            by = c("sp_j" = "species",
                   "timestep"), 
            suffix = c("_i", "_j")) %>% 
  dplyr::select(-log_r_j) %>% 
  rename(log_r = log_r_i,
         species = sp_i)

A <- A[unique(df_jags$species), unique(df_jags$species)]

df_jags <- df_jags %>% 
  mutate(species = as.numeric(factor(species)),
         species = factor(species),
         t = timestep - 1)

df_a <- tibble(value = c(A),
               row = which(!is.na(A), arr.ind = T)[, 1],
               col = which(!is.na(A), arr.ind = T)[, 2])


# plot --------------------------------------------------------------------

df_jags %>% 
  ggplot(aes(x = x0_j,
             y = log_r)) +
  geom_point() +
  facet_grid(rows = vars(species),
             cols = vars(sp_j),
             scales = "free")

list_dyn$df_dyn %>% 
  filter(species %in% c(1, 2)) %>% 
  ggplot(aes(x = timestep,
             y = density,
             color = factor(species))) +
  geom_line()
    