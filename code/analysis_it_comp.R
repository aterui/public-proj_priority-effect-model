
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))


# sim data ----------------------------------------------------------------

nsp <- 20
nt <- 20

# test --------------------------------------------------------------------

a <- seq(0, 1, by = 0.25)

df_b <- foreach(i = 1:length(a),
                .combine = bind_rows) %do% {
                  
                  print(i)
                  
                  df0 <- foreach(j = 1:100,
                                 .combine = bind_rows) %do% {
                                   
                                   v_r <- runif(nsp, 0.5, 2.5)
                                   k <- runif(nsp, 500, 500)
                                   
                                   A <- matrix(runif(nsp^2,
                                                     min = max(a[i] - 0.25, 0),
                                                     max = min(a[i] + 0.25, 1)),
                                               nsp,
                                               nsp)
                                   
                                   diag(A) <- 1
                                   
                                   #A[upper.tri(A)] <- 1
                                   
                                   if(a[i] == 1) {
                                     A <- matrix(runif(nsp^2,
                                                       min = max(a[i] - 0.05),
                                                       max = min(a[i] + 0.05)),
                                                 nsp,
                                                 nsp)
                                     
                                     diag(A) <- 1
                                     v_r <- rep(mean(v_r), nsp)
                                     k <- mean(k)
                                   }
                                   
                                   source(here::here("code/sim_data.R"))
                                   
                                   df_null <- df2 %>%
                                     filter(sp1 == sp2) %>% 
                                     group_by(t) %>%
                                     mutate(xt = sum(x_j),
                                            xt0 = sum(x0_j)) %>% 
                                     ungroup() %>% 
                                     mutate(p_x = x_i / xt,
                                            p_x0 = x0_i / xt0,
                                            log_pr = log(p_x) - log(p_x0),
                                            logit_p_x0 = log(p_x0) / (1 - log(p_x0)))
                                   
                                   df_p <- list_dyn$df_species %>% 
                                     dplyr::select(species, mean_density) %>% 
                                     mutate(p = mean_density / sum(mean_density)) %>% 
                                     rename(sp1 = species)
                                   
                                   df_lm <- df_null %>%
                                     group_by(sp1) %>%
                                     do(lm = MASS::rlm(log_r ~ p_x0, .) %>% coef()) %>%
                                     mutate(b0 = lm[1],
                                            b1 = lm[2]) %>%
                                     dplyr::select(-lm) %>%
                                     left_join(df_p,
                                               by = "sp1")
                                    
                                   fit <- lm(log(abs(b1)) ~ log(p), df_lm) %>% 
                                     summary()
                                      
                                   return(tibble(alpha = a[i],
                                                 r = cor(df_lm$b1, df_lm$p,
                                                          method = "spearman"),
                                                 rsq = fit$r.squared,
                                                 z = abs(coef(fit)[2, 1]))
                                          )

                                 }
                  
                  return(df0)
                }

df_b %>%
  pivot_longer(cols = c(r, rsq, z)) %>%
  ggplot(aes(x = factor(alpha),
             y = value)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_jitter(alpha = 0.2) +
  facet_wrap(facets = ~ name,
             scales = "free")
# 
# df_lm %>%
#   ggplot(aes(x = p,
#              y = -b1)) +
#   geom_point() +
#   scale_x_continuous(trans = "log10") +
#   scale_y_continuous(trans = "log10")
# 
