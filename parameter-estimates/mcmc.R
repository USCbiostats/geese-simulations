#!/bin/sh
#SBATCH --account=lc_pdt
#SBATCH --partition=thomas
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=12GB
#SBATCH --job-name=parameter_estimates
#SBATCH --mail-type=ALL

library(aphylo)
library(coda)

gen_starts <- function(x, n, f = .1) {
  x <- matrix(x, ncol = length(x), nrow = n, byrow = TRUE)
  
  upper_half <- x > .5
  lower_half <- which(!upper_half)
  upper_half <- which(upper_half)
  
  x[upper_half] <- x[upper_half] - runif(length(upper_half), min = f/3, max = f)
  x[lower_half] <- x[lower_half] + runif(length(lower_half), min = f/3, max = f)
  
  x
}

shrink_towards_half <- function(x, margin=.01) {
  
  x[x < (.5 - margin)] <- x[x < (.5 - margin)] + margin
  x[x > (.5 + margin)] <- x[x > (.5 + margin)] - margin
  
  x
}

# Common parameters
# prior.  <- bprior(c(2,2,9,9,2,2,2,2,2), c(9,9,2,2,9,9,9,9,9))
prior.  <- bprior(c(2,2,9,9,2,2,2), c(9,9,2,2,9,9,9))
warmup. <- 500
freq.   <- 10
ub.     <- c(1 - 1e-5)
lb.     <- c(1e-5)

mcmc. <- list(
  nchains      = 4L,
  multicore    = FALSE, 
  burnin       = 0L,
  nsteps       = 5000L,
  conv_checker = NULL,
  kernel       = fmcmc::kernel_ram(lb = lb., ub = ub., warmup = warmup., freq = freq.),
  thin         = 10L
)

set.seed(1362)

# Fitting partially annotated trees --------------------------------------------

partially_annotated <- readRDS("data/candidate_trees.rds")
partially_annotated <- do.call(c, partially_annotated)
partially_annotated <- unlist(lapply(partially_annotated, function(d) {
  lapply(1:Nann(d), function(i) d[,i])
}), recursive = FALSE)
partially_annotated <- do.call(c, partially_annotated)

# 1.A: No prior
ans_mle_partially_annotated_no_prior <- aphylo_mle(
  partially_annotated ~ psi + mu_d + mu_s + Pi
)

message("Partially annotated: MLE No prior done.")

saveRDS(
  ans_mle_partially_annotated_no_prior,
  "parameter-estimates/mle_partially_annotated_no_prior.rds"
  )

# In this case we don't need that many samples, this converges faster
mcmc.$nsteps <- 5000L

set.seed(173812)
mcmc.$kernel <- fmcmc::kernel_adapt(lb = lb., ub = ub., warmup = warmup., freq = 1L)
ans_mcmc_partially_annotated_no_prior <- aphylo_mcmc(
  partially_annotated ~ psi + mu_d + mu_s + Pi,
  params  = gen_starts(coef(ans_mle_partially_annotated_no_prior), mcmc.$nchains),
  control = mcmc.
)

message("Partially annotated: MCMC No prior done.")

saveRDS(
  ans_mcmc_partially_annotated_no_prior,
  "parameter-estimates/mcmc_partially_annotated_no_prior.rds"
  )
# 1.B: Prior
ans_mle_partially_annotated_prior <- aphylo_mle(
  partially_annotated ~ psi + mu_d + mu_s + Pi, priors = prior.
)

message("Partially annotated: MLE prior done.")

saveRDS(
  ans_mle_partially_annotated_prior,
  "parameter-estimates/mle_partially_annotated_prior.rds"
  )

set.seed(8812831)
mcmc.$kernel <- fmcmc::kernel_adapt(lb = lb., ub = ub., warmup = warmup., freq = 1L)
ans_mcmc_partially_annotated_prior <- aphylo_mcmc(
  partially_annotated ~ psi + mu_d + mu_s + Pi,
  priors = prior.,
  params  = gen_starts(coef(ans_mle_partially_annotated_prior), mcmc.$nchains),
  control = mcmc.
)

message("Partially annotated: MCMC prior done.")

saveRDS(
  ans_mcmc_partially_annotated_prior,
  "parameter-estimates/mcmc_partially_annotated_prior.rds"
  )


