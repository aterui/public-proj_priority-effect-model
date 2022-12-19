list_dyn <- cdyns::cdynsim(n_species = nsp,
                           n_timestep = 20,
                           r_type = "constant",
                           r = v_r,
                           int_type = "manual",
                           alpha = A,
                           k = k, 
                           sd_env = 0.1,
                           model = "ricker",
                           immigration = 5)

df0 <- list_dyn$df_dyn %>% 
  mutate(count = rpois(nrow(.), lambda = density))

df0 <- df0 %>% 
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
  ungroup()# %>% 
  # group_by(timestep) %>% 
  # mutate(x = mean(count0)) %>% 
  # ungroup()
