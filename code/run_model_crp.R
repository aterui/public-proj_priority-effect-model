
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))


# data --------------------------------------------------------------------

#set.seed(1)
nsp <- 20
k <- 100
v_r <- rnorm(nsp, 1, 0.5)
v_r[v_r < 0.5] <- 0.5
A <- matrix(rexp(nsp * nsp, rate = 1/1),
            nsp,
            nsp)

diag(A) <- 1
# x <- sample(1:nsp, 3)
# A[x, x] <- 1

df_a <- tibble(value = c(A),
               row = which(!is.na(A), arr.ind = T)[, "row"],
               col = which(!is.na(A), arr.ind = T)[, "col"])

#set.seed(1)
list_dyn <- cdyns::cdynsim(n_species = nsp,
                           n_timestep = 20,
                           r_type = "constant",
                           r = v_r,
                           int_type = "manual",
                           alpha = A,
                           k = k, 
                           sd_env = 0.1,
                           model = "ricker",
                           immigration = 0)

df0 <- list_dyn$df_dyn %>% 
  mutate(count = rpois(nrow(.), lambda = density))

df_jags <- df0 %>% 
  group_by(species) %>%
  summarize(index = all(count > 0)) %>%
  right_join(df0, by = "species") %>%
  ungroup() %>%
  filter(index == TRUE) %>%
  group_by(species) %>%
  mutate(count0 = lag(count),
         log_r = log(count) - log(count0)) %>%
  drop_na(log_r) %>% 
  mutate(t = timestep - 1) %>% 
  ungroup()

A <- A[unique(df_jags$species), unique(df_jags$species)]

df_jags <- df_jags %>% 
  mutate(species = as.numeric(factor(species)))

z <- matrix(NA, n_distinct(df_jags$species), n_distinct(df_jags$species))
diag(z) <- 1

# df_jags %>% 
#   ggplot(aes(x = count0,
#              y = log_r)) +
#   geom_point() +
#   facet_wrap(facets = ~species)


# common setup ------------------------------------------------------------

## mcmc setup ####
n_ad <- 100
n_iter <- 2E+4
n_thin <- max(3, ceiling(n_iter / 500))
n_burn <- ceiling(max(10, n_iter/2))
n_chain <- 4
n_sample <- ceiling(n_iter / n_thin)

inits <- replicate(n_chain,
                   list(.RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA),
                   simplify = FALSE)

for (j in 1:n_chain) inits[[j]]$.RNG.seed <- (j - 1) * 10 + 1

## parameters ####
para <- c("z",
          "beta",
          "sigma",
          "alpha",
          "mu_p")

# jags --------------------------------------------------------------------

d_jags <- list(Y = df_jags$log_r,
               X = df_jags$count0,
               Year = df_jags$t,
               Species = df_jags$species,
               Nsample = nrow(df_jags),
               Nsp = n_distinct(df_jags$species),
               Nyear = n_distinct(df_jags$t),
               z = z)

m <- read.jagsfile("code/model_crp.R")
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

rhat <- max(MCMCvis::MCMCsummary(post$mcmc)[,"Rhat"],
            na.rm = T)
print(rhat)

MCMCvis::MCMCsummary(post$mcmc) %>% round(2)

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



# # plot --------------------------------------------------------------------
# 
# df_plot <- MCMCvis::MCMCsummary(post$mcmc) %>% 
#   as_tibble(rownames = "param") %>% 
#   filter(str_detect(param, "z\\[\\d{1,},\\d{1,}\\]")) %>% 
#   dplyr::select(param, mean) %>% 
#   mutate(id = str_extract(param, "\\d{1,},\\d{1,}")) %>% 
#   separate(id,
#            into = c("row", "col"),
#            convert = T) %>% 
#   left_join(df_a,
#             by = c("row", "col"))
# 
# df_plot %>% 
#   filter(row != col) %>% 
#   ggplot(aes(x = value,
#              y = mean)) +
#   geom_point()
