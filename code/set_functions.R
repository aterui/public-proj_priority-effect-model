
run_nimble <- function(m, data, constants, inits, niter,
                       nburnin = NULL,
                       thin = NULL) {
  
  if(is.null(nburnin)) nburnin <- floor(niter / 2)
  if(is.null(thin)) thin <- floor((niter - nburn) / 500)
  
  nm <- nimble::nimbleModel(m,
                            constants = constants,
                            data = data,
                            inits = inits)
  nmconf <- configureMCMC(nm,
                          enableWAIC = TRUE,
                          monitors = nm$getNodeNames(stochOnly = TRUE,
                                                     includeData = FALSE))
  cnm <- nimble::compileNimble(nm)
  
  mcmc_m <- nimble::buildMCMC(nmconf)
  mcmc_cm <- nimble::compileNimble(mcmc_m, project = nm)
  
  mcmc_fit <- nimble::runMCMC(mcmc_cm,
                              niter = niter,
                              nburnin = nburnin,
                              thin = thin,
                              nchain = 3,
                              WAIC = TRUE)
  
  # posterior prob ----------------------------------------------------------
  
  log_ml <- bridgesampling::bridge_sampler(mcmc_cm,
                                           silent = TRUE)
  
  return(list(log_ml = log_ml,
              estimate = mcmc_fit))
}

