
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))
source(here::here("code/set_functions.R"))
source(here::here("code/data_fmt_fishdata.R"))


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
m <- read.jagsfile("code/model_multi_ricker_sparse.R")

## parameters ####
para <- c("p0",
          "log_r",
          "sigma_time",
          "sigma_obs",
          "sigma",
          "alpha",
          "rho")

# jags --------------------------------------------------------------------

## data screening
unique_site <- df_complete %>% 
  mutate(p = ifelse(abundance > 0, 1, 0)) %>% 
  group_by(site_id, taxon) %>% 
  summarize(freq = sum(p, na.rm = T),
            site_id = unique(site_id)) %>% 
  filter(freq > 4) %>% 
  group_by(site_id) %>% 
  summarize(n_taxa = n_distinct(taxon)) %>% 
  filter(n_taxa > 1) %>% 
  pull(site_id)

df_est <- foreach(i = seq_len(length(unique_site)),
                  .combine = bind_rows) %do% {
                    
                    ## subset data ####
                    df_subset <- df_complete %>% 
                      dplyr::filter(site_id == unique_site[i]) %>% 
                      mutate(p = ifelse(abundance > 0, 1, 0)) %>% 
                      group_by(taxon) %>% 
                      summarize(freq = sum(p, na.rm = T),
                                site_id = unique(site_id)) %>% 
                      filter(freq > 4) %>% 
                      dplyr::select(taxon, site_id) %>% 
                      left_join(df_complete,
                                by = c("taxon", "site_id")) %>% 
                      mutate(taxon_id = as.numeric(factor(.$taxon)))
                    
                    ## data for jags ####
                    d_jags <- list(N = df_subset$abundance,
                                   Year = df_subset$year - min(df_subset$year) + 1,
                                   Species = as.numeric(factor(df_subset$taxon)),
                                   Area = df_subset$area,
                                   Nsample = nrow(df_subset),
                                   Nyr1 = min(df_subset$year, na.rm = T) - 1999 + 1,
                                   Nyr = max(df_subset$year, na.rm = T) - min(df_subset$year, na.rm = T) + 1,
                                   Nsp = n_distinct(df_subset$taxon),
                                   Nf = 2,
                                   W = diag(n_distinct(df_subset$taxon)),
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
                      post <- extend.jags(post,
                                          burnin = 0,
                                          sample = n_sample,
                                          adapt = n_ad,
                                          thin = n_thin,
                                          n.sims = n_chain,
                                          combine = TRUE)
                      
                      mcmc_summary <- MCMCvis::MCMCsummary(post$mcmc)
                      print(max(mcmc_summary$Rhat, na.rm = T))
                    }
                    
                    ## reformat mcmc_summary ####
                    n_total_mcmc <- (post$sample / n_sample) * n_iter + n_burn
                    
                    mcmc_summary <- mcmc_summary %>% 
                      mutate(param = rownames(.),
                             site = unique_site[i]) %>% 
                      tibble() %>% 
                      rename(median = `50%`,
                             lower = `2.5%`,
                             upper = `97.5%`) %>% 
                      mutate(param_name = str_remove(param,
                                                     pattern = "\\[.{1,}\\]"),
                             x = fn_brrm(param),
                             n_total_mcmc = n_total_mcmc,
                             n_sample = post$sample,
                             n_thin = n_thin,
                             n_burn = n_burn,
                             n_chain = n_chain) %>% 
                      separate(col = x,
                               into = c("x1", "x2"),
                               sep = ",",
                               fill = "right",
                               convert = TRUE) %>% 
                      filter(!str_detect(param, "loglik")) %>% 
                      left_join(distinct(df_subset, taxon, taxon_id),
                                by = c("x1" = "taxon_id")) %>% 
                      left_join(distinct(df_subset, taxon, taxon_id),
                                by = c("x2" = "taxon_id"))
                    
                    ## trace plot output ####
                    MCMCvis::MCMCtrace(post$mcmc,
                                       filename = paste0("output/mcmc_trace_multi_ricker_sparse_",
                                                         unique_site[i]))
                    
                    print(paste(i, "/", length(unique_site)))
                    return(mcmc_summary)
                  }


# export ------------------------------------------------------------------

saveRDS(df_est,
        file = here::here("output/summary_multi_ricker_sparse.rds"))
