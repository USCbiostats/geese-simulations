library(aphylo)
library(geese)
library(coda)

# Fitting partially annotated trees --------------------------------------------
model_data <- readRDS("data/model_data.rds")
NSTEPS     <- 2e4

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

# Building the model
nfunctions <- 1

# Checkingout number of polytomies
parse_polytomies(model2fit)

# Building the model
term_overall_changes(model2fit, duplication = TRUE)
term_overall_changes(model2fit, duplication = FALSE)
term_genes_changing(model2fit, duplication = TRUE)
term_gains(model2fit, 0:(nfunctions - 1))
term_gains(model2fit, 0:(nfunctions - 1), FALSE)
term_loss(model2fit, 0:(nfunctions - 1), FALSE)

rule_limit_changes(model2fit, 0, 0, 4, TRUE)
rule_limit_changes(model2fit, 1, 0, 4, FALSE)

# Currently fails b/c some nodes have 7 offspring.
# 2^((7 + 1) * 4 functions) = 2 ^ 32 = 4,294,967,296 different cases
# which the computer cannot handle, unless using restrictions.
init_model(model2fit)
message("This model's support size is: ", support_size(model2fit))

fn <- "parameter-estimates/mcmc-joint-geese-prior.rds"

set.seed(112)

loc <- c(
  # Overall changes
  0, 0,
  # Genes changing at duplication
  -1/2,
  # Gains and loss x nfunctions (duplication)
  rep(1/2, nfunctions), # rep(-1/2, nfunctions),
  # Gains and loss x nfunctions (speciation)
  rep(-1/2, nfunctions), rep(-1/2, nfunctions),
  rep(0, nfunctions)
)

# loc <- c(0,0, -1/2, 1/2, -1/2, -1/2, -1,-9)
# loc <- c(rep(0, nterms(model2fit) - 1), -9)
if (!file.exists(fn)) {
  ans_geese_mcmc <- geese_mcmc(
    model2fit,
    initial = loc*0,
    prior  = function(p) dnorm(p, mean = loc, sd = 1, log = TRUE),
    nsteps = NSTEPS,
    kernel = fmcmc::kernel_ram(
      warmup = 1e3,
      fixed  = c(TRUE,TRUE, rep(FALSE, nterms(model2fit) - 2)),
      lb     = -9,
      ub     = 9
    ))

  estimates <- colMeans(window(ans_geese_mcmc, start = NSTEPS/2))

  pred <- predict_flock(model2fit, estimates, leave_one_out = TRUE)

  saveRDS(
    list(
      mcmc = ans_geese_mcmc,
      pred = pred,
      data = model_data[data_to_include]
    ),
    file = fn
    )

}

fn <- "parameter-estimates/mcmc-joint-geese-no-prior.rds"

if (!file.exists(fn)) {
  # No prior
  set.seed(1112)
  ans_geese_mcmc_no_prior <- geese_mcmc(
    model2fit,
    initial = loc*0,
    prior  = function(p) 0,
    nsteps = NSTEPS,
    kernel = fmcmc::kernel_ram(
      warmup = 1e3,
      fixed  = c(TRUE,TRUE, rep(FALSE, nterms(model2fit) - 2)),
      lb     = -9,
      ub     = 9
    ))


  estimates_no_prior <- colMeans(window(ans_geese_mcmc_no_prior, start = NSTEPS/2))

  pred_no_prior <- predict_flock(model2fit, estimates_no_prior, leave_one_out = TRUE)
  saveRDS(
    list(
      mcmc = ans_geese_mcmc_no_prior,
      pred = pred_no_prior,
      data = model_data[data_to_include]
    ),
    file = fn
    )

}

# Aphylo
set.seed(212)
adata <- do.call(c, lapply(model_data[data_to_include], "[[", "tree"))
ans_aphylo <- aphylo_mcmc(
  adata ~ mu_d + mu_s + psi + Pi,
  priors = bprior(c(2,2,9,5,2,2,5), c(9,9,2,5,9,9,5)),
)

auc_aphylo <- prediction_score(ans_aphylo, loo = TRUE)

set.seed(212)
ans_aphylo_no_prior <- aphylo_mcmc(
  adata ~ mu_d + mu_s + psi + Pi,
  priors = function(p) 1,
)

auc_aphylo_no_prior <- prediction_score(ans_aphylo_no_prior, loo = TRUE)

ans <- list(
#  mcmc            = ans_geese_mcmc,
#  mcmc_no_prior   = ans_geese_mcmc_no_prior,
#  pred            = pred,
#  pred_no_prior   = pred_no_prior,
  data            = model_data[data_to_include],
  aphylo          = ans_aphylo,
  aphylo_no_prior = ans_aphylo_no_prior
)

saveRDS(ans, "parameter-estimates/mcmc-joint-aphylo.rds")





