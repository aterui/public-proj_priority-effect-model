
list_dyn <- cdyns::cdynsim(n_species = nsp,
                           n_timestep = nt,
                           r_type = "constant",
                           r = v_r,
                           int_type = "manual",
                           alpha = A,
                           k = k,
                           seed = 100,
                           sd_env = 0.1,
                           model = "ricker",
                           immigration = 5)

df1 <- list_dyn$df_dyn %>%
  mutate(x = rpois(nrow(.), lambda = density)) %>%
  group_by(species) %>% 
  mutate(index = all(x > 0)) %>% 
  ungroup() %>% 
  filter(index == TRUE) %>% 
  group_by(species) %>%
  mutate(x0 = lag(x)) %>%
  drop_na(x0) %>%
  ungroup() %>% 
  mutate(t = timestep - 1,
         log_r = log(x) - log(x0),
         spf = as.numeric(factor(species))) %>% 
  dplyr::select(t, log_r, x, x0, spf) %>% 
  arrange(spf)

df2 <- expand(df1, t, sp1 = spf, sp2 = spf) %>%
  left_join(df1, by = c("sp1" = "spf",
                        "t")) %>%
  left_join(df1, by = c("sp2" = "spf",
                        "t"),
            suffix = c("_i", "_j")) %>%
  arrange(sp1) %>%
  rename(log_r = log_r_i)

