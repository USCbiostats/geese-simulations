#!/bin/sh
#SBATCH --account=dconti_251
#SBATCH --partition=conti
#SBATCH --mail-user=vegayon@usc.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=aphylo2-sim
#SBATCH --mem=20G

library(aphylo2)

n     <- 100
nsims <- 5e3
NJOBS <- 2555550

# Testing
params <- c(
  # Gaining function 0, and 0->1
  2.5, 1.5,
  # Overall changes
  1.5, -2,
  # Root probabilities
  -10, -10
)

names(params) <- c(
  "gain0", "neofun01", "changes_dpl", "changes_sp",
  "root0", "root1"
)

# Preparing data
annotations <- replicate(n * 2 - 1, c(9, 9), simplify = FALSE)

# Random tree
set.seed(31)
tree <- aphylo::sim_tree(n)$edge - 1L

duplication <- runif(n * 2 - 1) > .8

# Reading the data in
amodel <- new_model(
  annotations = annotations,
  geneid = c(tree[, 2], n),
  parent = c(tree[, 1], -1L),
  duplication = duplication
)

# Preparing the model
invisible({

  term_kgains(amodel, 0, 1, TRUE)
  term_neofun_a2b(amodel, 0, 1, TRUE)
  term_overall_changes(amodel, TRUE)
  term_overall_changes(amodel, FALSE)

  init(amodel)

})



likelihood(amodel, params*0)

# Simulating
set_seed(amodel, 111)
ans <- replicate(nsims, {
  sim_aphylo2(amodel, params)
}, simplify = FALSE)

# Finding MLE in each one of them
library(slurmR)
out <- slurmR::Slurm_lapply(ans, function(a) {

    # Building the model
    amodel <- new_model(
      annotations = a[c(tree[, 2], n) + 1],
      geneid      = c(tree[,2], n),
      parent      = c(tree[,1], -1),
      duplication
      )

    invisible({

      term_kgains(amodel, 0, 1, TRUE)
      term_neofun_a2b(amodel, 0, 1, TRUE)
      term_overall_changes(amodel, TRUE)
      term_overall_changes(amodel, FALSE)

      init(amodel)

    })

    # Fitting the model
    ans_mcmc <- tryCatch(aphylo2_mcmc(
      amodel,
      nsteps  = 40000,
      kernel  = fmcmc::kernel_ram(warmup = 2000),
      prior   = function(p) dlogis(p, scale = 2, log = TRUE)
    ), error = function(e) e)


    if (inherits(ans_mcmc, "error"))
      return(ans_mcmc)

    # Thinning
    ans_mcmc <- window(ans_mcmc, start = 30000)

    list(
      estimates = apply(ans_mcmc, 2, quantile, probs = c(.025, .5, 0.975)),
      vcov      = cov(ans_mcmc)
    )

  },
  njobs = NJOBS,
  job_name = "aphylo2-lapply",
  tmp_path = "/scratch/vegayon/",
  sbatch_opt    = list(
    account     = "dconti_251",
    partition   = "conti",
    "mail-user" = "vegayon@usc.edu",
    "mail-type" = "END,FAIL"
    ),
  mc.cores = 1L,
  export = c("tree", "duplication", "n"),
  overwrite = TRUE
)

# Saving the output
saveRDS(out, file = "simulation-study.rds")

