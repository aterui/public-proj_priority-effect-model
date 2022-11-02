
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))

# common setup ------------------------------------------------------------

## mcmc setup ####
n_ad <- 100
n_iter <- 1.0E+4
n_thin <- max(3, ceiling(n_iter / 250))
n_burn <- ceiling(max(10, n_iter/2))
n_chain <- 4
n_sample <- ceiling(n_iter / n_thin)

inits <- replicate(n_chain,
                   list(.RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA),
                   simplify = FALSE)

for (j in 1:n_chain) inits[[j]]$.RNG.seed <- (j - 1) * 10 + 1

## model file ####
m <- read.jagsfile("code/model_sparse_int.R")

## parameters ####
para <- c("loglik",
          "p0",
          "log_r",
          "sigma_time",
          "sigma_obs",
          "sigma",
          "alpha",
          "rho")

# jags --------------------------------------------------------------------

## psi = 0: niche-structured
## psi = 1: neutral

psi <- c(0, 1)
df_para <- expand.grid(n_species = 10,
                       n_timestep = c(5, 10, 15, 20),
                       r = 1,
                       alpha = c(1, 0.5),
                       k = 100,
                       sd_env = 0.1)

df_p <- foreach(i = seq_len(nrow(df_para)),
                .combine = bind_rows) %do% {
                  
                  x <- df_para[i,]
                  
                  list_re <- foreach(j = 1:2) %do% {
                    
                    # data --------------------------------------------------------------------
                    set.seed(1)
                    list_dyn <- cdyns::cdynsim(n_species = x$n_species,
                                               n_timestep = x$n_timestep,
                                               r_type = "constant",
                                               r = x$r,
                                               int_type = "constant",
                                               alpha = x$alpha,
                                               k = x$k, 
                                               sd_env = x$sd_env,
                                               model = "bh")
                    
                    df0 <- list_dyn$df_dyn %>% 
                      mutate(density_obs = density * exp(rnorm(n = nrow(.), mean = 0, sd = 0.1)),
                             count = rpois(nrow(.), lambda = density_obs))
                    
                    ## data for jags ####
                    d_jags <- list(N = df0$count,
                                   Year = df0$timestep,
                                   Species = df0$species,
                                   Nsample = nrow(df0),
                                   Nyr1 = 1,
                                   Nyr = max(df0$timestep),
                                   Nsp = n_distinct(df0$species),
                                   Nf = 2,
                                   W = diag(n_distinct(df0$species)),
                                   Q = 1,
                                   Psi = psi[j])
                    
                    ## run jags ####
                    post <- run.jags(m$model,
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
                                     module = "glm",
                                     silent.jags = TRUE)
                    
                    mcmc_summary <- MCMCvis::MCMCsummary(post$mcmc)
                    print(max(mcmc_summary$Rhat, na.rm = T))
                    
                    while(max(mcmc_summary$Rhat, na.rm = T) >= 1.15) {
                      post <- runjags::extend.jags(post,
                                                   burnin = 0,
                                                   sample = n_sample,
                                                   adapt = n_ad,
                                                   thin = n_thin,
                                                   n.sims = n_chain,
                                                   combine = TRUE,
                                                   silent.jags = TRUE)
                      
                      mcmc_summary <- MCMCvis::MCMCsummary(post$mcmc)
                      print(max(mcmc_summary$Rhat, na.rm = T))
                    }

                                        
                    # model evaluation --------------------------------------------------------
                    waic_bar <- MCMCvis::MCMCchains(post$mcmc) %>% 
                      as_tibble() %>% 
                      dplyr::select(starts_with("loglik")) %>% 
                      data.matrix() %>% 
                      waic()
                    
                    return(waic_bar)
                  }
                  
                  v_waic <- lapply(list_re, function(x) x$estimates[3,1]) %>%
                    unlist()
                  
                  x %>%
                    mutate(waic0 = v_waic[1],
                           waic1 = v_waic[2],
                           d_waic = waic1 - waic0)
                }


# plot --------------------------------------------------------------------

g_waic <- df_p %>% 
  ggplot(aes(x = n_timestep,
             y = d_waic,
             color = factor(alpha)),
         alpha = 0.8) +
  geom_hline(yintercept = 0,
             color = grey(0.8, 0.5)) +
  geom_point() +
  geom_line() +
  labs(color = "Competition ratio") +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave(g_waic,
       filename = here::here("output/figure_waic.pdf"),
       height = 5,
       width = 7.5)