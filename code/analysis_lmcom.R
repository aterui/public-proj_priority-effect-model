
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))

# parameters --------------------------------------------------------------

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
    mutate(xt = (x0_i + x0_j) / 2)
    
  h1 <- MASS::rlm(log_r ~ scale(x0_i),
                  data = df_lm,
                  method = "MM")
  
  h0 <- MASS::rlm(log_r ~ scale(xt),
                  data = df_lm,
                  method = "MM")
  
  p0 <- abs(coef(h0)[2] - coef(h1)[2])
  
  # p <- exp(-0.5 * BIC(h0)) + exp(-0.5 * BIC(h1))
  # p0 <- exp(-0.5 * BIC(h0)) / p
  # p1 <- exp(-0.5 * BIC(h1)) / p
  
  M[combo[1, i], combo[2, i]] <- p0
}


plot(M ~ Am)
