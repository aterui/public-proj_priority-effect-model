
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))


# data --------------------------------------------------------------------

source(here::here("code/run_sim_data.R"))


# common setup ------------------------------------------------------------

## mcmc setup ####
n_ad <- 100
n_iter <- 3.0E+4
n_thin <- max(3, ceiling(n_iter / 1000))
n_burn <- ceiling(max(10, n_iter/2))
n_chain <- 4
n_sample <- ceiling(n_iter / n_thin)

inits <- replicate(n_chain,
                   list(.RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA),
                   simplify = FALSE)

for (j in 1:n_chain) inits[[j]]$.RNG.seed <- (j - 1) * 10 + 1

## parameters ####
para <- c("r",
          "b",
          "z")

# jags --------------------------------------------------------------------

z <- matrix(NA, nsp, nsp)
diag(z) <- 0

d_jags <- list(Y = df_jags$log_r,
               X = df_jags$x0_j,
               Year = df_jags$t,
               Sp_i = df_jags$species,
               Sp_j = df_jags$sp_j,
               Nsample = nrow(df_jags),
               Nsp = n_distinct(df_jags$species),
               Nyear = n_distinct(df_jags$t),
               z = z)

m <- read.jagsfile("code/model_lm.R")
post <- suppressMessages(run.jags(m$model,
                                  monitor = para,
                                  data = d_jags,
                                  n.chains = n_chain,
                                  inits = inits,
                                  method = "parallel",
                                  burnin = n_burn,
                                  sample = n_sample,
                                  adapt = n_ad,
                                  thin = n_thin,
                                  n.sims = 4,
                                  module = "glm"))

mcmc_summary <- MCMCvis::MCMCsummary(post$mcmc)

rhat <- max(mcmc_summary[,"Rhat"],
            na.rm = T)
print(rhat)

# while(rhat >= 1.1) {
#   post <- extend.jags(post,
#                       combine = TRUE,
#                       sample = n_sample,
#                       thin = n_thin)
#   
#   rhat <- max(MCMCvis::MCMCsummary(post$mcmc)[,"Rhat"],
#               na.rm = T)
#   print(rhat)
# }

df_pr <- mcmc_summary %>% 
  as_tibble(rownames = "param") %>% 
  filter(str_detect(param, "z")) %>% 
  mutate(index = str_extract(param, pattern = "\\d{1,},\\d{1,}")) %>% 
  separate(index,
           into = c("row", "col"),
           sep = ",",
           convert = T) %>% 
  dplyr::select(row,
                col,
                pr = mean) %>% 
  left_join(df_a,
            by = c("row", "col")) %>% 
  filter(row != col)

df_pr %>% 
  ggplot(aes(x = value,
            y = pr)) +
  geom_point() +
  scale_y_continuous(limits = c(0, 1))
