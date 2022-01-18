library(aphylo)
library(geese)
library(coda)

# Fitting partially annotated trees --------------------------------------------
model_data <- readRDS("data/model_data.rds")
NSTEPS     <- 2e4

# init <- c(
#   `Overall changes at duplication` = 0,
#   `Overall changes at speciation` = 0,
#   `Only one gene changes at duplication` = -.5,
#   `Only one gene changes at speciation` = -.5,
#   `Gains 0 at duplication` = 0.0768957326737147,
#   `Gains 0 at speciation` = -2.15269078801566,
#   `Loss 0 at speciation` = -3.42075856213656,
#   `Loss 0 at duplication` = -3.42075856213656,
#   `Root 1` = 1.64267185610882
# ) * 0

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

term_k_genes_changing(model2fit, 1, TRUE)
term_k_genes_changing(model2fit, 1, FALSE)

term_gains(model2fit, 0:(nfunctions - 1), TRUE)
term_gains(model2fit, 0:(nfunctions - 1), FALSE)

term_loss(model2fit, 0:(nfunctions - 1), TRUE)
term_loss(model2fit, 0:(nfunctions - 1), FALSE)

rule_limit_changes(model2fit, 0, 0, 4, TRUE)
rule_limit_changes(model2fit, 1, 0, 4, FALSE)

# Currently fails b/c some nodes have 7 offspring.
# 2^((7 + 1) * 4 functions) = 2 ^ 32 = 4,294,967,296 different cases
# which the computer cannot handle, unless using restrictions.
init_model(model2fit)

print(model2fit)

set.seed(112)

loc <- c(
  # Overall changes
  0, 0,
  # One gene changing
  c(-1/2, -1/2),
  # Gains and loss x nfunctions (duplication)
  rep(1/2, nfunctions), # rep(-1/2, nfunctions),
  # Gains and loss x nfunctions (speciation)
  rep(-1/2, nfunctions), rep(-1/2, nfunctions),
  rep(0, nfunctions)
)

fn <- "parameter-estimates/mcmc-joint-geese.rds"

set.seed(1112)
ans_geese_mcmc <- geese_mcmc(
  model2fit,
  prior  = function(p) 0, # Uniform prior
  nsteps = 20000, burnin = 0, thin = 1,
  kernel = fmcmc::kernel_am(
    warmup = 2e3,
    fixed  = c(TRUE,TRUE, rep(FALSE, nterms(model2fit) - 2)),
    lb     = -10,
    ub     = 10
  ))

estimates_no_prior <- colMeans(
  window(ans_geese_mcmc, start = NSTEPS * (4/5))
  )

pred_no_prior <- predict_flock(
  p   = model2fit,
  par = estimates_no_prior,
  leave_one_out = TRUE
  )

saveRDS(
  list(
    mcmc = ans_geese_mcmc,
    pred = pred_no_prior,
    data = model_data[data_to_include]
  ),
  file = fn
  )






