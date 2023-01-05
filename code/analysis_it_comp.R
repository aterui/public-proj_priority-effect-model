
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))


# sim data ----------------------------------------------------------------

nsp <- 10
nn <- 1
nt <- 20
v_r <- rep(1.5, nsp)
b <- rep(0.001, nsp)
k <- v_r / b

# test --------------------------------------------------------------------

a <- 0:10 * 0.1

df_b <- foreach(i = 1:length(a),
                .combine = bind_rows) %do% {
                  
                  print(i)
                  
                  df0 <- foreach(j = 1:100,
                                 .combine = bind_rows) %do% {
                                   
                                   A <- matrix(abs(rnorm(nsp^2,
                                                         mean = a[i],
                                                         sd = 0)),
                                               nsp,
                                               nsp)
                                   
                                   
                                   # if(a[i] == 1) {
                                   #   A <- matrix(abs(rnorm(nsp^2,
                                   #                         mean = 1,
                                   #                         sd = 0)),
                                   #               nsp,
                                   #               nsp)
                                   # } else {
                                   #   A <- matrix(abs(rnorm(nsp^2,
                                   #                         mean = a[i],
                                   #                         sd = 0.1)),
                                   #               nsp,
                                   #               nsp)
                                   # }
                                   
                                   diag(A) <- 1
                                   
                                   source(here::here("code/sim_data.R"))
                                   
                                   df_null <- df2 %>%
                                     filter(sp1 == sp2) %>% 
                                     group_by(t) %>%
                                     mutate(xt = mean(x_j),
                                            xt0 = mean(x0_j)) %>% 
                                     ungroup() %>% 
                                     group_by(sp1) %>% 
                                     mutate(scl_x0_i = scale(x0_i) %>% c(),
                                            scl_xt0 = scale(xt0) %>% c()) %>% 
                                     ungroup()
                                   
                                   beta0 <- lm(log_r ~ scl_x0_i,
                                               df_null) %>% 
                                     summary() %>% 
                                     coef()
                                   
                                   beta1 <- lm(log_r ~ scl_xt0,
                                               df_null) %>%
                                     summary() %>% 
                                     coef()
                                   
                                   return(tibble(alpha = a[i],
                                                 b0 = beta0[2, 1],
                                                 b1 = beta1[2, 1],
                                                 se0 = beta0[2, 2],
                                                 se1 = beta1[2, 2]))
                                   
                                 }
                  
                  return(df0)
                }


df_b %>% 
  mutate(mean = b1 / b0) %>% 
  pivot_longer(cols = c(mean, se1)) %>% 
  ggplot(aes(x = factor(alpha),
             y = value)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_jitter(alpha = 0.2) +
  facet_wrap(facets = ~ name,
             scales = "free")
