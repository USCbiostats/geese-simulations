library(aphylo)
library(geese)
library(coda)

# Fitting partially annotated trees --------------------------------------------
model_data <- readRDS("data/model_data.rds")

# Subsetting trees to obtain parameter estimates
data_to_include <- which(
  (sapply(model_data, "[[", "max_size") < 1e9/2) &
    (sapply(model_data, "[[", "nfuns") == 1)
  )
data_to_include <- names(model_data)[data_to_include]

model2fit <- new_flock()
for (i in data_to_include) {

  with(model_data[[i]], add_geese(
    model2fit,
    annotations = ann,
    geneid      = edges[,1],
    parent      = edges[,2],
    duplication = dpl
    )
  )

}

################################################################################
# Model Fit --------------------------------------------------------------------
################################################################################

# Aphylo
adata <- do.call(c, lapply(model_data[data_to_include], "[[", "tree"))
ans_aphylo <- aphylo_mcmc(
  adata ~ mu_d + mu_s + psi + Pi,
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

rule_limit_changes(model2fit, 0, 0, 5, TRUE)
rule_limit_changes(model2fit, 1, 0, 5, FALSE)


# Currently fails b/c some nodes have 7 offspring.
# 2^((7 + 1) * 4 functions) = 2 ^ 32 = 4,294,967,296 different cases
# which the computer cannot handle, unless using restrictions.
init_model(model2fit)

set.seed(112)

# Prior
loc <- c(
  0, 0, -1/2,
  rep(1/2, nfunctions),
  rep(-1/2, nfunctions),
  rep(-1/2, nfunctions),
  rep(-1/2, nfunctions),
  rep(-4, nfunctions)
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
    # fixed  = c(TRUE,TRUE, rep(FALSE, nterms(model2fit) - 3), rep(TRUE, nfunctions)),
    lb     = -10,
    ub     = 10
  ))

estimates <- colMeans(window(ans_geese_mcmc, start = 1e4))

pred <- predict_flock(model2fit, estimates, leave_one_out = TRUE)

ans <- list(
  mcmc     = ans_geese_mcmc,
  pred     = pred,
  data     = model_data[data_to_include],
  aphylo   = ans_aphylo
)

saveRDS(ans, "parameter-estimates/mcmc-joint-geese3.rds")





