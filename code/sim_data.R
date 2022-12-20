list_dyn <- cdyns::cdynsim(n_species = nsp,
                           n_timestep = nt,
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
  summarize(index = sum(count > 0) / nt) %>%
  right_join(df0, by = "species") %>%
  ungroup() %>%
  filter(index == 1) %>%
  group_by(species) %>%
  mutate(count0 = lag(count)) %>%
  drop_na(count0) %>% 
  mutate(t = timestep - 1) %>% 
  ungroup()