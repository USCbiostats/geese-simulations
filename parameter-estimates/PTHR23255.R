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
    functions = length(d$ann[[1]]),
    nzeros    = sum(unlist(d$ann) == 0L),
    nones     = sum(unlist(d$ann) == 1L),
    nann      = sum(unlist(d$ann) != 9L)
  )

  ans <- c(ans, expected_size = 2^((ans[1] + 1) * ans[2]))

})


# Resorting to start from the easiest
data_features <- data_features[,order(data_features[3,])]


current_tree <- "PTHR10788"

# Model 0 ----------------------------------------------------------------------
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
# term_genes_changing(model2fit, duplication = TRUE)
term_gains(model2fit, funs = 0:1, TRUE)
term_loss(model2fit, funs = 0:1, TRUE)

rule_limit_changes(model2fit, 0, 0, 4, TRUE)
rule_limit_changes(model2fit, 1, 0, 4, FALSE)

init_model(model2fit)
message("This model's support size is: ", support_size(model2fit))

set.seed(112)
loc <- c(
  rep(0, nterms(model2fit) - 2),
  rep(-9, 2)
)

model0 <- geese_mcmc(
  model2fit,
  prior  = function(p) dlogis(p, location = loc, scale = 2, log = TRUE),
  nsteps = 1e4,
  kernel = fmcmc::kernel_ram(
    warmup = 1e3,
    lb     = -10,
    ub     = 10
  ))


plot(window(model0, start=8e3))

# Computing MAE
pred0 <- predict_geese(
  model2fit,
  colMeans(window(model0, start=8e3)), leave_one_out = TRUE
)

pred0 <- do.call(rbind, pred0)
obse <- do.call(rbind, data[[current_tree]]$ann)

pscore0 <- prediction_score(
  x = pred0, expected = obse
)

plot(pscore0$auc)

# Model 1 ----------------------------------------------------------------------
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
term_gains(model2fit, funs = 0:1, TRUE)
term_neofun_a2b(model2fit, 0, 1, TRUE)
term_subfun(model2fit, 0, 1)

rule_limit_changes(model2fit, 0, 0, 4, TRUE)
rule_limit_changes(model2fit, 1, 0, 4, FALSE)

init_model(model2fit)
message("This model's support size is: ", support_size(model2fit))

set.seed(112)
loc <- c(
  rep(0, nterms(model2fit) - 2),
  rep(-9, 2)
)

model1 <- geese_mcmc(
  model2fit,
  prior  = function(p) dlogis(p, location = loc, scale = 2, log = TRUE),
  nsteps = 1e4,
  kernel = fmcmc::kernel_ram(
    warmup = 1e3,
    lb     = -10,
    ub     = 10
  ))


plot(window(model1, start=8e3))

# Computing MAE
pred1 <- predict_geese(
  model2fit,
  colMeans(window(model1, start=8e3)), leave_one_out = TRUE
  )

pred1 <- do.call(rbind, pred1)
obse <- do.call(rbind, data[[current_tree]]$ann)

pscore1 <- prediction_score(
  x = pred1, expected = obse
)

plot(pscore1$auc)

# Model 2 ----------------------------------------------------------------------
model2fit <- with(data[[ current_tree ]], new_geese(
  annotations = ann,
  geneid      = tree[,1],
  parent      = tree[,2],
  duplication = dpl
))

# Checkingout number of polytomies
parse_polytomies(model2fit)

# Building the model
term_overall_changes(model2fit, duplication = TRUE)
term_overall_changes(model2fit, duplication = FALSE)
term_genes_changing(model2fit, duplication = TRUE)
term_gains(model2fit, funs = 0:1, TRUE)
term_loss(model2fit, funs = 0:1, TRUE)

rule_limit_changes(model2fit, 0, 0, 4, TRUE)
rule_limit_changes(model2fit, 1, 0, 4, FALSE)

init_model(model2fit)
message("This model's support size is: ", support_size(model2fit))

set.seed(112)
loc <- c(
  rep(0, nterms(model2fit) - 2),
  rep(-9, 2)
)

model2 <- geese_mcmc(
  model2fit,
  prior  = function(p) dlogis(p, location = loc, scale = 2, log = TRUE),
  nsteps = 1e4,
  kernel = fmcmc::kernel_ram(
    warmup = 1e3,
    lb     = -10,
    ub     = 10
  ))


plot(window(model2, start=8e3))

# Computing MAE
pred2 <- predict_geese(
  model2fit,
  colMeans(window(model2, start=8e3)), leave_one_out = TRUE
)

pred2 <- do.call(rbind, pred2)

pscore2 <- prediction_score(
  x = pred2, expected = obse
)

plot(pscore2$auc)
