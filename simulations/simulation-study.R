#!/bin/sh
#SBATCH --account=dconti_251
#SBATCH --partition=conti
#SBATCH --mail-user=vegayon@usc.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=aphylo2-sim
#SBATCH --mem=20G

library(aphylo)
library(geese)

n     <- 100
nsims <- 5e3
NJOBS <- 100

# Simulates tree and duplication events
sim_tree <- function(n) {

  annotations <- replicate(n * 2 - 1, c(9, 9), simplify = FALSE)

  # Random tree
  tree <- aphylo::sim_tree(n)$edge - 1L

  duplication <- c(
    rep(FALSE, n),
    sample.int(2, n - 1, replace = TRUE, prob = c(.8, .2)) == 2
  )

  list(
    tree = tree,
    duplication = duplication
  )

}

# Sets the model
set_model <- function(dat, n) {

  ann <- if (!length(dat$ann))
    replicate(n * 2 - 1, c(9, 9), simplify = FALSE)
  else
    dat$ann

  m <- new_geese(
    annotations = ann,
    geneid = c(dat$tree[, 2], n),
    parent = c(dat$tree[, 1], -1L),
    duplication = dat$duplication
  )

  # Persistence to preserve parent state
  term_overall_changes(m, duplication = 0)
  term_overall_changes(m, duplication = 1)

  # Tendency for only one sib gaining function 0
  term_pairwise_neofun_singlefun(m, 0, duplication = 1)

  # Function 1 is gain from function 0
  term_neofun_a2b(m, 0, 1, duplication = 1)

  init_model(m, verb = 0)

  m

}

# Randomly drop data
drop_ann <- function(ann, ids, m) {

  ann[sample(ids, m, replace = FALSE)] <- rep(list(rep(9L, length(ann[[1]]))), m)
  ann

}

# Random tree
set.seed(3151)
tree   <- sim_tree(n)
amodel <- set_model(tree, n)

# View(do.call(rbind, get_support(amodel)))

# Testing
params <- c(
  # Persist the parent state
  -4, -4,
  # Gain function 0
  3,
  # Neofun 0 -> 1
  5,
  # Root probabilities
  -10, -10
)
names(params) <- c(names(amodel), "root0", "root1")

likelihood(amodel, params*0)

# Simulating
set_seed(amodel, 111)
ans <- replicate(nsims, {

  # Setup and simulation
  atree  <- sim_tree(n)
  amodel <- set_model(atree, n)
  tmp    <- sim_geese(amodel, params)

  # Drop data
  list(
    tree        = atree$tree,
    ann         = drop_ann(tmp, which(atree$tree[,2] < n), 80),
    duplication = atree$duplication
  )


}, simplify = FALSE)


# Finding MLE in each one of them
# library(slurmR)
# opts_slurmR$debug_on()
cl <- parallel::makePSOCKcluster(6)
parallel::clusterEvalQ(cl, {
  library(geese)
  library(aphylo)
  })

parallel::clusterExport(cl, c("n", "set_model", "ans"))

out <- parallel::parLapply(cl, seq_along(ans), function(i) {

    fn <- sprintf("simulations/simulation-study-%04i.rds", i)
    if (file.exists(fn))
      return(0)

    a <- ans[[i]]

    # Building the model
    tmpmodel <- set_model(a, n = n)

    # Fitting the model
    ans_mcmc <- tryCatch(geese_mcmc(
      tmpmodel,
      nsteps  = 20000,
      kernel  = fmcmc::kernel_ram(warmup = 2000),
      prior   = function(p) dnorm(p, sd = 2, log = TRUE)
    ), error = function(e) e)


    if (inherits(ans_mcmc, "error"))
      return(ans_mcmc)

    # Thinning
    ans_mcmc <- window(ans_mcmc, start = end(ans_mcmc) * 3/4)

    # Prediction
    preds <- predict_geese(
      tmpmodel,
      par           = colMeans(ans_mcmc),
      leave_one_out = TRUE,
      use_reduced_sequence = TRUE
      )

    score <- prediction_score(
      do.call(rbind, preds),
      do.call(rbind, a$ann)
      )

    # Model and prediction with aphylo
    ap <- aphylo_from_data_frame(
      tree        = as.phylo(a$tree),
      annotations = data.frame(
        id = c(a$tree[, 2], n),
        do.call(rbind, a$ann)
      ),
      types = data.frame(0:(2*n-2), as.integer(!a$duplication))
    )

    ans_aphylo <- aphylo_mcmc(ap ~ mu_d + mu_s + Pi)

    res <- list(
      estimates   = apply(ans_mcmc, 2, quantile, probs = c(.025, .5, 0.975)),
      vcov        = cov(ans_mcmc),
      score_geese = score,
      score_aphylo = prediction_score(ans_aphylo)
    )

    saveRDS(res, fn)

    return(0)

  } # ,
  # njobs = NJOBS,
  # job_name = "aphylo2-lapply",
  # # tmp_path = "/scratch/vegayon/",
  # sbatch_opt    = list(
  #   account     = "dconti_251",
  #   partition   = "conti",
  #   "mail-user" = "vegayon@usc.edu",
  #   "mail-type" = "END,FAIL"
  #   ),
  # mc.cores = 4L #,
  # export = c("tree", "duplication", "n", "set_model"),
  # overwrite = TRUE
)

parallel::stopCluster(cl)

# Saving the output
# saveRDS(out, file = "simulations/simulation-study.rds")

