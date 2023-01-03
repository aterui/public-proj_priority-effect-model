
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))


# parameters --------------------------------------------------------------

#set.seed(1)
nsp <- 10
nt <- 20
nn <- 3
v_r <- c(rep(1.5, nn), runif(nsp - nn, 1, 2.5))
b <- c(rep(0.002, nn), runif(nsp - nn, 0, 0.002))
k <- v_r / b

subsp <- c(2, 3)

# h01 <- foreach(i = 1:100, .combine = bind_rows) %do% {

for(i in 1:100) {  
  print(i)
  
  A <- matrix(0.01,
              nsp,
              nsp)
  
  A[1:nn, 1:nn] <- 1
  A[(nn + 1):nsp, (nn + 1):nsp] <- runif(length((nn + 1):nsp)^2, 0, 0.5)
  diag(A) <- 1
  
  # simulate ----------------------------------------------------------------
  
  sp <- 1
  while(sp < 5) {
    source("code/sim_data.R")
    (sp <- n_distinct(df0$species))
  }
  
  # analysis ----------------------------------------------------------------
  
  df1 <- df0 %>%
    filter(species %in% subsp) %>% 
    mutate(x = count,
           x0 = count0,
           spf = case_when(species == subsp[1] ~ 1,
                           species == subsp[2] ~ 2),
           spf = factor(spf),
           species = factor(species)) %>% 
    pivot_wider(id_cols = c(t, x, x0, species),
                names_prefix = "sp",
                names_from = spf,
                values_from = count0)
  
  df1 <- df1 %>% 
    mutate(sp1 = rep(na.omit(df1$sp1), 2),
           sp2 = rep(na.omit(df1$sp2), 2),
           xt = sp1 + sp2)
  
  h0 <- glm(x ~ xt + offset(log(x0)),
            family = "poisson",
            data = df1)
  
  h1 <- glm(x ~ x0 * species + offset(log(x0)),
            family = "poisson",
            data = df1)
  
  dbic <- BIC(h0) - BIC(h1)
  p <- exp(-0.5 * BIC(h0)) + exp(-0.5 * BIC(h1))
  p0 <- exp(-0.5 * BIC(h0)) / p
  p1 <- exp(-0.5 * BIC(h1)) / p
  
  lr <- c(exp(logLik(h0) - logLik(h1)))
  
  if(p0 < 0.05) stop("er")
  return(tibble(dbic, lr, p0))
}

g1 <- h01 %>% 
  ggplot(aes(x = p0,
             y = dbic)) + 
  geom_point() +
  theme_bw()

g2 <- h01 %>% 
  ggplot(aes(x = log(lr))) +
  geom_density(fill = grey(0, 0.1)) +
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  theme_bw()

g3 <- h01 %>% 
  ggplot(aes(x = dbic)) +
  geom_density(fill = grey(0, 0.1)) +
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  theme_bw()

g4 <- h01 %>% 
  ggplot(aes(x = p0)) +
  geom_density(fill = grey(0, 0.1)) +
  geom_vline(xintercept = 0.05,
             linetype = "dashed") +
  theme_bw()

(g1 + g2) / (g3 + g4)
mean(h01$p0 > 0.05)
