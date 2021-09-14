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
  fn <- sprintf("parameter-estimates/mcmc-unif-prior-%s.rds", current_tree)
  if (file.exists(fn)) {
    message("This tree was already analyzed...")
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
  term_gains(model2fit, 0:(nfunctions - 1))
  # term_loss(model2fit, 0:(nfunctions - 1))
  term_gains(model2fit, 0:(nfunctions - 1), FALSE)
  term_loss(model2fit, 0:(nfunctions - 1), FALSE)

  rule_limit_changes(model2fit, 0, 0, 4, TRUE)
  rule_limit_changes(model2fit, 1, 0, 4, FALSE)

  # Currently fails b/c some nodes have 7 offspring.
  # 2^((7 + 1) * 4 functions) = 2 ^ 32 = 4,294,967,296 different cases
  # which the computer cannot handle, unless using restrictions.
  init_model(model2fit)
  message("This model's support size is: ", support_size(model2fit))

  set.seed(112)

  loc <- c(
    # Overall changes
    0, 0,
    # Genes changing at duplication
    -1/2,
    # Gains and loss x nfunctions (duplication)
    rep(1/2, nfunctions), # rep(-1/2, nfunctions),
    # Gains and loss x nfunctions (speciation)
    rep(-1/2, nfunctions), rep(-1/2, nfunctions) #,
    # rep(0, nfunctions)
  )

  ans_geese_mcmc <- geese_mcmc(
    model2fit,
    prior  = function(p) 0, # Uniform prior
    nsteps = 20000, burnin = 0,
    kernel = fmcmc::kernel_am(
      warmup = 2e3,
      fixed  = c(TRUE,TRUE, rep(FALSE, nterms(model2fit) - 2)),
      lb     = -10,
      ub     = 10
    ))

  # Making predictions
  estimates <- colMeans(window(ans_geese_mcmc, start = 10000))
  pred_loo <- predict_geese(model2fit, estimates, leave_one_out = TRUE)
  pred_loo <- unlist(pred_loo)

  auc_geese <- aphylo::prediction_score(
    x   = cbind(pred_loo),
    expected = cbind(unlist(data[[ current_tree ]]$ann))
  )

  # How about the baseline model?
  atree <- adata[[ current_tree ]]
  ans_aphylo <- aphylo_mcmc(
    atree ~ mu_s + mu_d + Pi,
    priors = uprior()
  )

  auc_aphylo <- prediction_score(ans_aphylo, loo = TRUE)

  output <- list(
    geese_mcmc  = ans_geese_mcmc,
    geese_auc   = auc_geese,
    aphylo_mcmc = ans_aphylo,
    aphylo_auc  = auc_aphylo,
    tree        = partially_annotated[[ current_tree ]]
  )

  saveRDS(output, file = fn)

}

