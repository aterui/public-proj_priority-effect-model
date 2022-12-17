
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))

# data --------------------------------------------------------------------
## parameters
set.seed(1)

nyear <- 20
nsp <- 30
k <- 100
alpha <- runif(nsp * nsp, 0, 0.5)

b1 <- runif(nsp, 1, 2.5)
b2 <- - b1 / k
b3 <- b2 * alpha

df_t <- tibble(param = c(rep("b1", nsp),
                         rep("b2", nsp)),
               true = c(b1, b2),
               row = rep(1:nsp, 2))

A <- matrix(alpha, nsp, nsp)
A[1:5, ] <- 0.5
A[1:5, 1:5] <- 1
diag(A) <- 1

## simulate
source("code/sim_data.R")
print(n_distinct(df0$species))

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
m <- read.jagsfile("code/model_ssm.R")

## parameters ####
para <- c("b",
          "SIGMA")

# jags --------------------------------------------------------------------
## data for jags ####
d_jags <- list(Y = df0$count,
               Year = df0$timestep,
               Species = df0$species,
               Nsample = nrow(df0),
               Nt = max(df0$timestep),
               Nsp = n_distinct(df0$species),
               K = 2)

## run jags ####
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

(mcmc_summary <- MCMCvis::MCMCsummary(post$mcmc))
print(max(mcmc_summary$Rhat, na.rm = T))


# plot --------------------------------------------------------------------

df_plot <- mcmc_summary %>%
  as_tibble(rownames = "param") %>%
  filter(str_detect(param, "b")) %>%
  mutate(id = str_extract(param, "\\d{1,},\\d{1,}")) %>%
  separate(id, into = c("row", "col"),
           convert = T) %>%
  dplyr::select(row,
                col,
                mean,
                median = `50%`,
                low = `2.5%`,
                high = `97.5%`) %>%
  mutate(param = case_when(col == 1 ~ "b1",
                           col == 2 ~ "b2"),
         g = case_when(row <= 5 ~ "g1",
                       row > 5 ~ "g2")) %>%
  left_join(df_t, by = c("param",
                         "row"))

ggplot(df_plot,
       aes(x = median,
           color = g,
           fill = g)) +
  geom_histogram() +
  facet_wrap(facets = ~param,
             scales = "free")
