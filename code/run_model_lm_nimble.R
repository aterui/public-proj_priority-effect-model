
# setup -------------------------------------------------------------------

source(here::here("code/library.R"))
source(here::here("code/set_functions.R"))


# data --------------------------------------------------------------------

source(here::here("code/sim_data.R"))

data <- list(Y = df_jags$log_r,
             X = df_jags$count0)

const <- list(Nsample = length(df_jags$log_r),
              Nyear = n_distinct(df_jags$t),
              Nsp = n_distinct(df_jags$species),
              Species = df_jags$species,
              Year = df_jags$t,
              K = 2,
              VSD = rep(2.5, 2))

inits <- list(mu_beta = c(0, 0.1),
              m_u = diag(const$K),
              sigma = 0.1)


# nimble model ------------------------------------------------------------

source(here::here("code/model_nimble_lm_h0.R"))
source(here::here("code/model_nimble_lm_h1.R"))

niter <- 2E+4

ml0 <- run_nimble(m = m0,
                  data = data,
                  constants = const,
                  inits = inits,
                  niter = niter)

ml1 <- run_nimble(m = m1,
                  data = data,
                  constants = const,
                  inits = inits,
                  niter = niter)

