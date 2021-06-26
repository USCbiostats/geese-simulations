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

# Subsetting trees to obtain parameter estimates
data_to_include <- which(data_features[3,] < 1e9/2)
data_to_include <- colnames(data_features)[data_to_include]

model2fit <- new_flock()
adata2    <- NULL
data_included <- NULL
for (i in data_to_include) {

  if (length(data[[i]]$ann[[1]]) > 1)
    next

  with(data[[i]], add_geese(
    model2fit,
    annotations = ann,
    geneid      = tree[,1],
    parent      = tree[,2],
    duplication = dpl
    )
  )

  if (!length(adata2))
    adata2 <- adata[[i]]
  else
    adata2 <- c(adata2, adata[[i]])

  data_included <- c(data_included, i)
}

data_obs <- data[data_included]

################################################################################
# Model Fit --------------------------------------------------------------------
################################################################################

# Aphylo
ans_aphylo <- aphylo_mcmc(
  adata2 ~ mu_d + mu_s + psi + Pi,
  priors = bprior(c(2,2,9,5,2,2,5), c(9,9,2,5,9,9,5)),
)

auc_aphylo <- prediction_score(ans_aphylo, loo = TRUE)
# stop()

# Building the model
nfunctions <- 1
term_overall_changes(model2fit, duplication = TRUE)
term_overall_changes(model2fit, duplication = FALSE)
term_genes_changing(model2fit, duplication = TRUE)
term_gains(model2fit, 0:(nfunctions - 1))
term_loss(model2fit, 0:(nfunctions - 1))
term_gains(model2fit, 0:(nfunctions - 1), FALSE)
term_loss(model2fit, 0:(nfunctions - 1), FALSE)

rule_limit_changes(model2fit, 0, 0, 4, TRUE)
rule_limit_changes(model2fit, 1, 0, 4, FALSE)


# Currently fails b/c some nodes have 7 offspring.
# 2^((7 + 1) * 4 functions) = 2 ^ 32 = 4,294,967,296 different cases
# which the computer cannot handle, unless using restrictions.
init_model(model2fit)

set.seed(112)

# Prior
loc <- c(
  -1/2, -1/2, -1/2,
  rep(1/2, nfunctions),
  rep(-1/2, nfunctions),
  rep(-1/2, nfunctions),
  rep(-1/2, nfunctions),
  rep(-9, nfunctions)
)

# loc <- c(0,0, -1/2, 1/2, -1/2, -1/2, -1,-9)
# loc <- c(rep(0, nterms(model2fit) - 1), -9)
ans_geese_mcmc <- geese_mcmc(
  model2fit,
  initial = loc,
  prior  = function(p) dnorm(p, mean = loc, sd = 1, log = TRUE),
  nsteps = 2e4,
  kernel = fmcmc::kernel_am(
    warmup = 1e3,
    fixed  = c(TRUE,TRUE, rep(FALSE, nterms(model2fit) - 3), rep(TRUE, nfunctions)),
    lb     = -10,
    ub     = 10
  ))

estimates <- colMeans(window(ans_geese_mcmc, start = 1e4))

pred <- predict_flock(model2fit, estimates)

ans <- list(
  mcmc     = ans_geese_mcmc,
  pred     = pred,
  data     = data_obs,
  included = data_included,
  aphylo   = ans_aphylo
)

saveRDS(ans, "parameter-estimates/mcmc-joint-geese2.rds")





