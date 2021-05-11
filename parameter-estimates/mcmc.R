#!/bin/sh
#SBATCH --account=lc_pdt
#SBATCH --partition=thomas
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=12GB
#SBATCH --job-name=parameter_estimates
#SBATCH --mail-type=ALL

library(aphylo)
library(geese)
library(coda)

shrink_towards_half <- function(x, margin=.01) {

  x[x < (.5 - margin)] <- x[x < (.5 - margin)] + margin
  x[x > (.5 + margin)] <- x[x > (.5 + margin)] - margin

  x
}

# Common parameters
warmup. <- 500
freq.   <- 10

mcmc. <- list(
  nchains      = 4L,
  multicore    = FALSE,
  burnin       = 0L,
  nsteps       = 5000L,
  conv_checker = NULL,
  kernel       = fmcmc::kernel_ram(warmup = warmup., freq = freq.),
  thin         = 10L
)

set.seed(1362)

# Fitting partially annotated trees --------------------------------------------

partially_annotated <- readRDS("data/candidate_trees.rds")

# Parsing the data
treeids <- sort(unique(names(partially_annotated)))
data    <- vector("list", length(treeids))
adata   <- data # aphylo data
names(data) <- treeids
for (tn in treeids) {

  # Preparing annotations
  tmp_trees <- partially_annotated[names(partially_annotated) == tn]
  tmp_ann   <- lapply(tmp_trees, function(a) rbind(a$tip.annotation, a$node.annotation))
  tmp_ann   <- do.call(cbind, tmp_ann)

  # Creating the aphylo version
  
  if (sum(names(partially_annotated) == tn) > 1)
    adata[[tn]] <- do.call(c, tmp_trees)
  else
    adata[[tn]] <- tmp_trees[[1L]]

  # Offspring -> parent, we need to add the root
  tmp_tree  <- rbind(tmp_trees[[1]]$tree$edge[, 2:1], c(Ntip(tmp_trees[[1]]) + 1, -1))

  # Sorting the annotations accordingly
  tmp_ann <- tmp_ann[tmp_tree[,1],,drop=FALSE]
  tmp_ann <- lapply(1:nrow(tmp_ann), function(i) tmp_ann[i ,])

  data[[tn]] <- list(
    tree = tmp_tree, ann = tmp_ann,
    dpl  = with(tmp_trees[[1L]], c(tip.type, node.type))[tmp_tree[,1]] == 0L
    )

}

# Obtaining features to figure out whether we can deal with it
data_features <- sapply(data, function(d) {

  tab <- table(table(d$tree[,2]))
  ans <- c(
    polytomes = max(as.integer(names(tab))),
    functions = length(d$ann[[1]])
  )

  ans <- c(ans, expected_size = 2^((ans[1] + 1) * ans[2]))

  })


# Resorting to start from the easiest
data_features <- data_features[,order(data_features[3,])]

for (current_tree in colnames(data_features)) {

  # Will fit a model only no more than 1/2 billion
  # combinations
  if (data_features[3, current_tree] >= (1e9/2))
    next

  message(paste(rep("#", options("width")), collapse = ""))
  message("Analyzing tree ", current_tree)
  message(paste(rep("#", options("width")), collapse = ""))

  # Checking if tree was already analyzed
  fn <- sprintf("parameter-estimates/mcmc-%s.rds", current_tree)
  if (file.exists(fn)) {
    message("This tree was already analyzed...")
    next
  }

  model2fit <- with(data[[ current_tree ]], new_geese(
    annotations = ann,
    geneid      = tree[,1],
    parent      = tree[,2],
    duplication = dpl
  ))

  nfunctions <- length(data[[current_tree]]$ann[[1]])

  # Checkingout number of polytomies
  parse_polytomies(model2fit)

  # Building the model
  term_overall_changes(model2fit, duplication = TRUE)
  term_overall_changes(model2fit, duplication = FALSE)
  term_genes_changing(model2fit, duplication = TRUE)
  term_genes_changing(model2fit, duplication = FALSE)
  term_gains(model2fit, 0:(nfunctions - 1))
  term_loss(model2fit, 0:(nfunctions - 1))

  rule_limit_changes(model2fit, 0, 0, 4, TRUE)
  rule_limit_changes(model2fit, 1, 0, 4, FALSE)

  # Currently fails b/c some nodes have 7 offspring.
  # 2^((7 + 1) * 4 functions) = 2 ^ 32 = 4,294,967,296 different cases
  # which the computer cannot handle, unless using restrictions.
  init_model(model2fit)
  message("This model's support size is: ", support_size(model2fit))

  set.seed(112)

  loc <- c(0,0,-.5,-1,rep(.5,nfunctions), rep(-1,nfunctions), rep(-1, nfunctions))

  ans_geese_mcmc <- geese_mcmc(
    model2fit,
    prior  = function(p) dlogis(p, location = loc, scale = 2, log = TRUE),
    nsteps = 2e4,
    kernel = fmcmc::kernel_am(
      warmup = 5e3,
      fixed  = c(TRUE,TRUE, rep(FALSE, nterms(model2fit) - 2))
    ))	     

  # Making predictions
  estimates <- colMeans(window(ans_geese_mcmc, start = 1e4))
  pred_loo <- predict_geese(model2fit, estimates, leave_one_out = TRUE)
  pred_loo <- unlist(pred_loo)

  auc_geese <- aphylo::auc(
    pred   = pred_loo,
    labels = unlist(data[[ current_tree ]]$ann)
  )

  # How about the baseline model?
  ans_aphylo <- aphylo_mcmc(
    adata[[ current_tree ]] ~ mu_d + mu_s + psi + Pi,
    priors = bprior(c(2,2,9,5,2,2,5), c(9,9,2,5,9,9,5))
    )

  auc_aphylo <- prediction_score(ans_aphylo, loo = TRUE)

  output <- list(
    geese_mcmc  = ans_geese_mcmc,
    geese_auc   = auc_geese,
    geese_pred  = pred_loo,
    aphylo_mcmc = ans_aphylo,
    aphylo_auc  = auc_aphylo$auc,
    aphylo_pred = auc_aphylo$predicted,
    labels      = unlist(data[[ current_tree ]]$ann),
    tree        = partially_annotated[[ current_tree ]]
  )

  saveRDS(output, file = fn)

}

#
#
# # Looking at pooled models -----------------------------------------------------
#
# model2 <- new_flock()
#
# for (i in 1:4)
#   with(data[[ least_annotated[i] ]], add_geese(
#     model2,
#     annotations = ann,
#     geneid      = tree[,1],
#     parent      = tree[,2],
#     duplication = dpl
#   ))
#
# term_overall_changes(model2, duplication = TRUE)
# term_overall_changes(model2, duplication = FALSE)
# term_genes_changing(model2, duplication = TRUE)
# term_gains(model2, 0)
# term_loss(model2, 0)
#
# init_model(model2)
#
# set.seed(112)
# ans2 <- geese_mcmc(
#   model2,
#   prior  = function(p) dlogis(p, location = c(1,-1,-1,1,-1,-2), scale = 2, log = TRUE),
#   nsteps = 20000,
#   kernel = fmcmc::kernel_ram(warmup = 5000)
# )
# plot(ans2)
#
# ans2 <- geese_mcmc(model2, nsteps = 20000, kernel = fmcmc::last_kernel())
# graphics.off()
# plot(window(ans2, start = 15000))
#
# # Making predictions
# estimates2 <- colMeans(window(ans2, start = 15000))
# pred_loo2 <- predict_flock(model2, estimates2, leave_one_out = TRUE)
# pred_loo2 <- unlist(pred_loo2)
#
# auc_geese <- aphylo::auc(
#   pred   = pred_loo,
#   labels = unlist(data[[ least_annotated[1L] ]]$ann)
# )
#
#
# partially_annotated <- do.call(c, partially_annotated)
# partially_annotated <- unlist(lapply(partially_annotated, function(d) {
#   lapply(1:Nann(d), function(i) d[,i])
# }), recursive = FALSE)
# partially_annotated <- do.call(c, partially_annotated)
#
# # 1.A: No prior
# ans_mle_partially_annotated_no_prior <- aphylo_mle(
#   partially_annotated ~ psi + mu_d + mu_s + Pi
# )
#
# message("Partially annotated: MLE No prior done.")
#
# saveRDS(
#   ans_mle_partially_annotated_no_prior,
#   "parameter-estimates/mle_partially_annotated_no_prior.rds"
#   )
#
# # In this case we don't need that many samples, this converges faster
# mcmc.$nsteps <- 5000L
#
# set.seed(173812)
# mcmc.$kernel <- fmcmc::kernel_adapt(lb = lb., ub = ub., warmup = warmup., freq = 1L)
# ans_mcmc_partially_annotated_no_prior <- aphylo_mcmc(
#   partially_annotated ~ psi + mu_d + mu_s + Pi,
#   params  = gen_starts(coef(ans_mle_partially_annotated_no_prior), mcmc.$nchains),
#   control = mcmc.
# )
#
# message("Partially annotated: MCMC No prior done.")
#
# saveRDS(
#   ans_mcmc_partially_annotated_no_prior,
#   "parameter-estimates/mcmc_partially_annotated_no_prior.rds"
#   )
# # 1.B: Prior
# ans_mle_partially_annotated_prior <- aphylo_mle(
#   partially_annotated ~ psi + mu_d + mu_s + Pi, priors = prior.
# )
#
# message("Partially annotated: MLE prior done.")
#
# saveRDS(
#   ans_mle_partially_annotated_prior,
#   "parameter-estimates/mle_partially_annotated_prior.rds"
#   )
#
# set.seed(8812831)
# mcmc.$kernel <- fmcmc::kernel_adapt(lb = lb., ub = ub., warmup = warmup., freq = 1L)
# ans_mcmc_partially_annotated_prior <- aphylo_mcmc(
#   partially_annotated ~ psi + mu_d + mu_s + Pi,
#   priors = prior.,
#   params  = gen_starts(coef(ans_mle_partially_annotated_prior), mcmc.$nchains),
#   control = mcmc.
# )
#
# message("Partially annotated: MCMC prior done.")
#
# saveRDS(
#   ans_mcmc_partially_annotated_prior,
#   "parameter-estimates/mcmc_partially_annotated_prior.rds"
#   )
#
#
