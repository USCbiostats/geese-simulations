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

# Subsetting trees to obtain parameter estimates
data_to_include <- which(data_features[3,] < 1e9/2)
data_to_include <- colnames(data_features)[data_to_include]

model2fit <- new_flock()
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
}

# Building the model
#term_overall_changes(model2fit, duplication = TRUE)
#term_overall_changes(model2fit, duplication = FALSE)
term_genes_changing(model2fit, duplication = TRUE)
term_genes_changing(model2fit, duplication = FALSE)
term_gains(model2fit, 0, TRUE)
term_loss(model2fit, 0, TRUE)
term_gains(model2fit, 0, FALSE)
term_loss(model2fit, 0, FALSE)

rule_limit_changes(model2fit, 2, 0, 4, TRUE)
rule_limit_changes(model2fit, 3, 0, 4, FALSE)
rule_limit_changes(model2fit, 4, 0, 4, TRUE)
rule_limit_changes(model2fit, 5, 0, 4, FALSE)


# Currently fails b/c some nodes have 7 offspring.
# 2^((7 + 1) * 4 functions) = 2 ^ 32 = 4,294,967,296 different cases
# which the computer cannot handle, unless using restrictions.
init_model(model2fit)

loc <- c(-.5,-1,.5,-1,-1, -1,-5)

set.seed(1231)
ans_geese_mcmc <- geese_mcmc(
  model2fit,
  prior  = function(p) dlogis(p, location = loc, scale = 2, log = TRUE),
  nsteps = 1e3,
  kernel = fmcmc::kernel_am(
    warmup = 500,
    # fixed  = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
    ub     = 10,
    lb     = -10
  )
)

