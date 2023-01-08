
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))


# sim data ----------------------------------------------------------------

nsp <- 10
nt <- 20

# test --------------------------------------------------------------------

a <- seq(0, 1, by = 0.25)

df_b <- foreach(i = 1:length(a),
                .combine = bind_rows) %do% {
                  
                  print(i)
                  
                  df0 <- foreach(j = 1:50,
                                 .combine = bind_rows) %do% {
                                   
                                   v_r <- runif(nsp, 1.5, 1.5)#rep(1.5, nsp)
                                   k <- runif(nsp, 1000, 1000)
                                   
                                   A <- matrix(runif(nsp^2,
                                                     min = max(a[i] - 0.1, 0),
                                                     max = min(a[i] + 0.1, 1)),
                                               nsp,
                                               nsp)
                                   
                                   #A[upper.tri(A)] <- 1
                                   diag(A) <- 1
                                   # 
                                   # if(a[i] == 1) {
                                   #   A[,] <- 1
                                   #   v_r <- rep(mean(v_r), nsp)
                                   #   b <- mean(b)
                                   #   k <- v_r / b
                                   #}
                                   
                                   source(here::here("code/sim_data.R"))
                                   
                                   df_null <- df2 %>%
                                     filter(sp1 == sp2) %>% 
                                     group_by(t) %>%
                                     mutate(xt0 = sum(x0_j)) %>% 
                                     ungroup() %>% 
                                     mutate(p_x = x0_i / xt0,
                                            logit_p_x = log(p_x) / (1 - log(p_x)))
                                   
                                   df_p <- list_dyn$df_species %>% 
                                     dplyr::select(species, mean_density) %>% 
                                     mutate(p = mean_density / sum(mean_density)) %>% 
                                     rename(sp1 = species)
                                   
                                   df_lm <- df_null %>%
                                     group_by(sp1) %>%
                                     do(lm = lm(log_r ~ p_x, .) %>% coef()) %>%
                                     mutate(b0 = lm[1],
                                            b1 = lm[2]) %>%
                                     dplyr::select(-lm) %>%
                                     left_join(df_p,
                                               by = "sp1") %>%
                                     mutate(w_b1 = p * b1)
                                    
                                   fit <- lm(log(abs(b1)) ~ log(p), df_lm) %>% 
                                     summary()
                                      
                                   return(tibble(alpha = a[i],
                                                 r1 = cor(df_lm$b1, df_lm$p,
                                                          method = "spearman"),
                                                 r2 = cor(log(abs(df_lm$b1)), log(df_lm$p),
                                                          method = "pearson"),
                                                 rsq1 = fit$r.squared,
                                                 rsq2 = fit$adj.r.squared)
                                          )

                                 }
                  
                  return(df0)
                }

df_b %>%
  pivot_longer(cols = c(r1, rsq1)) %>%
  ggplot(aes(x = factor(alpha),
             y = value)) +
  #geom_violin(draw_quantiles = 0.5) +
  geom_boxplot() +
  #geom_point(alpha = 0.2) +
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
