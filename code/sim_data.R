set.seed(1)
list_dyn <- cdyns::cdynsim(n_species = nsp,
                           n_timestep = nyear,
                           r_type = "constant",
                           r = b1,
                           int_type = "manual",
                           alpha = A,
                           k = k, 
                           sd_env = 0.1,
                           model = "ricker",
                           immigration = 5)

## data sort
df0 <- list_dyn$df_dyn %>% 
  mutate(count = rpois(nrow(.), lambda = density))

df0 <- df0 %>% 
  group_by(species) %>% 
  mutate(index = all(count > 0)) %>% 
  filter(index == TRUE) %>% 
  ungroup() %>% 
  mutate(species = as.numeric(factor(species)))

