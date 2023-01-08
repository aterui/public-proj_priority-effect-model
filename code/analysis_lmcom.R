
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

A <- matrix(0.01,
            nsp,
            nsp)

A[1:nn, 1:nn] <- 1
A[(nn + 1):nsp, (nn + 1):nsp] <- runif(length((nn + 1):nsp)^2, 0, 0.5)
diag(A) <- 1

sp <- 1
while(sp < 3) {
  source("code/sim_data.R")
  (sp <- n_distinct(df1$spf))
}

print(n_distinct(df1$spf))


# neutral index -----------------------------------------------------------

Am <- sqrt(A * t(A))
Am <- Am[unique(df1$spf), unique(df1$spf)]


# analysis ----------------------------------------------------------------

M <- matrix(NA, n_distinct(df1$spf), n_distinct(df1$spf))
combo <- combn(n_distinct(df1$spf), 2)

for (i in 1:ncol(combo)) {

  subsp <- combo[, i]
  
  df_lm <- df2 %>% 
    filter(sp1 %in% subsp[1],
           sp2 %in% subsp[2]) %>% 
    mutate(xt = (x0_i + x0_j),
           p_x = x0_i / xt)
    
  fit <- lm(log_r ~ p_x,
            data = df_lm)
  
  M[combo[1, i], combo[2, i]] <- coef(fit)[2]
}


plot(M ~ Am)
