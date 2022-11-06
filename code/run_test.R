
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))

# common setup ------------------------------------------------------------

## mcmc setup ####
n_ad <- 100
n_iter <- 5.0E+3
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
m <- read.jagsfile("code/model_var_int.R")

## parameters ####
para <- c("p",
          "log_r",
          "sigma",
          "sigma_r",
          "sigma_obs",
          "alpha_prime")

# jags --------------------------------------------------------------------

df_para <- expand.grid(n_species = c(2, 5),
                       n_timestep = c(5, 10, 20, 40),
                       r = log(10),
                       alpha = c(1, 0.1),
                       k = 100,
                       sd_env = 0.1) %>% 
  mutate(i = row_number())

n_rep <- 5

df_out <- foreach(x = iterators::iter(df_para, by = "row"),
                  .combine = bind_rows) %do% {
                    
                    print(paste("i = ", x$i))
                    
                    foreach(j = 1:n_rep,
                            .combine = bind_rows) %do% {
                              
                              print(paste("j = ", j))
                              set.seed(j)
                              
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
                                mutate(count = rpois(nrow(.), lambda = density))
                              
                              df0 %>% 
                                group_by(species) %>% 
                                summarize(n_obs = sum(count > 0))
                              
                              df0 %>%
                                ggplot(aes(x = timestep,
                                           y = count,
                                           color = factor(species))) +
                                geom_line() +
                                geom_point()
                              
                              ## data for jags ####
                              d_jags <- list(N = df0$count,
                                             Year = df0$timestep,
                                             Species = df0$species,
                                             Nsample = nrow(df0),
                                             Nyr1 = 1,
                                             Nyr = max(df0$timestep),
                                             Nsp = n_distinct(df0$species),
                                             W = diag(n_distinct(df0$species)),
                                             Q = 1)
                              
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
                                               module = "glm")
                              
                              mcmc_summary <- MCMCvis::MCMCsummary(post$mcmc)
                              print(max(mcmc_summary$Rhat, na.rm = T))
                              
                              while(max(mcmc_summary$Rhat, na.rm = T) >= 1.1) {
                                post <- suppressMessages(runjags::extend.jags(post,
                                                                              burnin = 0,
                                                                              sample = n_sample,
                                                                              adapt = n_ad,
                                                                              thin = n_thin,
                                                                              n.sims = n_chain,
                                                                              combine = TRUE,
                                                                              silent.jags = TRUE))
                                
                                mcmc_summary <- MCMCvis::MCMCsummary(post$mcmc)
                                print(max(mcmc_summary$Rhat, na.rm = T))
                              }
                              
                              x %>% 
                                mutate(n_rep = j,
                                       p = mcmc_summary[1, "50%"]) %>% 
                                relocate(n_rep) %>% 
                                return()
                              
                            }
                  }


# plot --------------------------------------------------------------------

df_out %>% 
  mutate(log_odd = log(p / (1 - p))) %>% 
  ggplot(aes(x = n_timestep,
             y = p,
             color = factor(n_species))) +
  geom_point() +
  geom_line() +
  facet_wrap(facets = ~alpha,
             labeller = label_both) +
  labs(color = 'Number of species',
       y = "Pr(neutralilty)",
       x = "Timestep") +
  theme_bw()

ggsave(filename = "output/test.pdf",
       height = 5,
       width = 12)
