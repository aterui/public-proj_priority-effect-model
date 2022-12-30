
# setup -------------------------------------------------------------------

library(tidyverse)


# parameter ---------------------------------------------------------------

nsp <- 100
nt <- 100
v_r <- runif(nsp, 0.5, 2.5)
b <- runif(nsp, 0, 0.002)
k <- v_r / b

A <- matrix(rexp(nsp^2, 2),
            nsp,
            nsp)

diag(A) <- 1

sp <- sample(1:nsp, 2)

# simulate ----------------------------------------------------------------

list_dyn <- cdyns::cdynsim(alpha = A,
                           int_type = "manual",
                           n_timestep = nt, 
                           n_species = nsp,
                           r = v_r,
                           k = k,
                           model = "bh",
                           extinct = 0.01,
                           immigration = 5)

df_p <- list_dyn$df_dyn %>% 
  group_by(species) %>% 
  mutate(sum = sum(density),
         index = sum > 0) %>% 
  dplyr::select(-sum) %>% 
  ungroup()

df_sum <- df_p %>% 
  filter(index == TRUE,
         !(species %in% sp)) %>% 
  group_by(timestep) %>% 
  summarize(n = sum(density)) %>%
  ungroup()

df_sum %>% 
  left_join(df_p, by = "timestep") %>% 
  filter(species %in% sp) %>% 
  ggplot(aes(x = n,
             y = density)) +
  geom_point() +
  facet_wrap(facets = ~species,
             scales = "free")
