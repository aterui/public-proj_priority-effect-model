
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))


# sim data ----------------------------------------------------------------

nsp <- 10
nn <- 1
nt <- 20
v_r <- runif(nsp, 0.5, 2.5)#rep(1.5, nsp)
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
                                                         sd = 0.1)),
                                               nsp,
                                               nsp)
                                   
                                   diag(A) <- 1
                                   
                                   source(here::here("code/sim_data.R"))
                                   
                                   df_null <- df2 %>%
                                     filter(sp1 == sp2) %>% 
                                     group_by(t) %>%
                                     mutate(xt = mean(x_j),
                                            xt0 = mean(x0_j),
                                            p_x = x0_i / xt0) %>% 
                                     ungroup() %>% 
                                     group_by(sp1) %>% 
                                     mutate(scl_x0_i = scale(x0_i) %>% c(),
                                            scl_xt0 = scale(xt0) %>% c()) %>% 
                                     ungroup()
                                   
                                   fit <- glmmTMB::glmmTMB(log_r ~ p_x + (1 | sp1),
                                                           df_null)
                                   
                                   beta <- fit %>% summary() %>% coef
                                   
                                   return(tibble(alpha = a[i],
                                                 b0 = beta$cond[2, 1],
                                                 se0 = beta$cond[2, 2]))
                                   
                                 }
                  
                  return(df0)
                }


df_b %>% 
  pivot_longer(cols = c(b0, se0)) %>% 
  ggplot(aes(x = factor(alpha),
             y = abs(value))) +
  geom_violin(draw_quantiles = 0.5) +
  geom_jitter(alpha = 0.2) +
  facet_wrap(facets = ~ name,
             scales = "free")
