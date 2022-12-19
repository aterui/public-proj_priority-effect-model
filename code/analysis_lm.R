
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))


# parameters --------------------------------------------------------------

#set.seed(1)
nsp <- 10
v_r <- c(1.5, 1.5, runif(nsp - 2, 0.5, 2.5))
b <- c(0.02, 0.02, runif(nsp - 2, 0, 0.02))
k <- v_r / b
A <- matrix(0.1,
            nsp,
            nsp)

A[1:2, 1:2] <- 1
diag(A) <- 1

h01 <- foreach(i = 1:100, .combine = c) %do% {
  
  print(i)
  
  # simulate ----------------------------------------------------------------
  
  sp <- 1
  while(sp < nsp) {
    source("code/sim_data.R")
    (sp <- n_distinct(df0$species))
  }
  
  # analysis ----------------------------------------------------------------
  
  df1 <- df0 %>%
    filter(species %in% c(1, 2)) %>% 
    group_by(timestep) %>% 
    mutate(x = mean(count0)) %>% 
    ungroup() %>% 
    mutate(species = as.numeric(factor(species)),
           spf = factor(species)) %>% 
    group_by(species) %>% 
    ungroup()
  
  h0 <- lm(log_r ~ x, df1)
  h1 <- lm(log_r ~ x * spf, df1)
  
  h01 <- AIC(h0) - AIC(h1)
  
  return(h01)
}

mean(h01 > 2)
plot(density(h01))
