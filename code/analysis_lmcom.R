
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))

# parameters --------------------------------------------------------------

#set.seed(1)
nsp <- 10
nt <- 20
v_r <- runif(nsp, 1.5, 1.5)
b <- runif(nsp, 0, 0.002)
k <- v_r / b

A <- matrix(rexp(nsp^2, 1/0.5), nsp, nsp)
diag(A) <- 1

sp <- 1
while(sp < 5) {
  source("code/sim_data.R")
  (sp <- n_distinct(df0$species))
}

print(n_distinct(df0$species))


# neutral index -----------------------------------------------------------

Am <- sqrt(A * t(A))
Am <- Am[unique(df0$species), unique(df0$species)]

# x <- rowSums(A)
# y <- colSums(A)
# m_d <- cbind(v_r, x, y)
# v_mean <- colMeans(m_d)
# m_cov <- cov(m_d)
# 
# m_dist <- dist(m_d) %>% 
#   data.matrix()
# 
# m_dist <- m_dist[unique(df0$species), unique(df0$species)]

# analysis ----------------------------------------------------------------

M <- matrix(NA, n_distinct(df0$species), n_distinct(df0$species))
combo <- combn(n_distinct(df0$species), 2)

for (i in 1:ncol(combo)) {

  subsp <- combo[, i]
  
  df1 <- df0 %>%
    mutate(species = factor(species) %>% as.numeric()) %>% 
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

  M[combo[1, i], combo[2, i]] <- BIC(h0) - BIC(h1)
}

M[M > 0] <- 0
M[M < 0] <- 1


plot(M ~ Am)
