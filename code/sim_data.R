
#set.seed(1)
nsp <- 10
nn <- 3
nt <- 20
v_r <- c(rep(1.5, nn), runif(nsp - nn, 0.5, 2.5))
b <- c(rep(0.001, nn), runif(nsp - nn, 0.0005, 0.001))
k <- v_r / b

A <- matrix(1,
            nsp,
            nsp)

A[1:nn, 1:nn] <- 1
A[(nn + 1):nsp, (nn + 1):nsp] <- runif(length((nn + 1):nsp)^2, 0, 0.2)
diag(A) <- 1

#set.seed(2)
list_dyn <- cdyns::cdynsim(n_species = nsp,
                           n_timestep = nt,
                           r_type = "constant",
                           r = v_r,
                           int_type = "manual",
                           alpha = A,
                           k = k,
                           seed = 100,
                           sd_env = 0.05,
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

# df2 %>% 
#   ggplot(aes(x = x0_j,
#              y = log_r,
#              color = factor(sp1))) +
#   geom_point(alpha = 0.1) +
#   facet_grid(rows = vars(sp1),
#              cols = vars(sp2),
#              scales = "free") +
#   geom_smooth(method = "lm")


# # test --------------------------------------------------------------------

y <- 3

df_b <- foreach(j = 2:n_distinct(df1$spf),
        .combine = bind_rows) %do% {
          
          df_lm <- df2 %>%
            filter(sp1 == 1,
                   sp2 == j) %>%
            group_by(t) %>%
            mutate(xt = (x0_i + x0_j) / 2)
          
          beta0 <- MASS::rlm(log_r ~ x0_i, df_lm) %>%
            coef()
          
          beta1 <- MASS::rlm(log_r ~ xt, df_lm) %>%
            coef()
          
          return(tibble(j = j,
                        r0 = beta0[1],
                        r1 = beta1[2],
                        b0 = beta0[2],
                        b1 = beta1[2]))
        }

df_b %>% 
  mutate(y = b1/b0) %>% 
  ggplot(aes(x = j,
             y = y)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 1,
             color = "salmon")
