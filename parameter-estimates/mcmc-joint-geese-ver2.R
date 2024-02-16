library(aphylo)
library(geese)
library(coda)
source("fig/plot_functions.R")

# Fitting partially annotated trees --------------------------------------------
model_data <- readRDS("data/model_data.rds")
NSTEPS     <- 4e4/2

# Subsetting trees to obtain parameter estimates
data_to_include <- which(
  (sapply(model_data, "[[", "max_size") < 1e8) &
  # (sapply(model_data, "[[", "polytomies") <= 25) &
    (sapply(model_data, "[[", "nfuns") == 1)
  )
length(data_to_include)

cumsum(table(sapply(model_data, "[[", "polytomies")))

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
polytoms <- parse_polytomies(model2fit)

# Building the model
term_overall_changes(model2fit, duplication = TRUE)  # Just constrain support
term_overall_changes(model2fit, duplication = FALSE) # Just constrain support

term_gains(model2fit, 0:(nfunctions - 1))
term_gains(model2fit, 0:(nfunctions - 1), FALSE)
term_loss(model2fit, 0:(nfunctions - 1))
term_loss(model2fit, 0:(nfunctions - 1), FALSE)

# Indicator variable (this makes the difference)
term_k_genes_changing(model2fit, 1, TRUE)
term_k_genes_changing(model2fit, 1, FALSE)

rule_limit_changes(model2fit, 0, 0, 4, TRUE)
rule_limit_changes(model2fit, 1, 0, 4, FALSE)

loc <- c(
  # Overall changes
  0, 0,
  # Genes changing at duplication
  # -1/2,
  # Gains and loss x nfunctions (duplication)
  rep(1/2, nfunctions), rep(-1/2, nfunctions),
  # Gains and loss x nfunctions (speciation)
  rep(-1/2, nfunctions), rep(-1/2, nfunctions),
  c(1/2, 1/2)
)


time0 <- proc.time()

# Currently fails b/c some nodes have 7 offspring.
# 2^((7 + 1) * 4 functions) = 2 ^ 32 = 4,294,967,296 different cases
# which the computer cannot handle, unless using restrictions.
init_model(model2fit)

View(support <- do.call(rbind, get_support(model2fit)))
cor(support[,-c(1:3)])


print(model2fit)

fn <- "parameter-estimates/mcmc-joint-geese.rds"

set.seed(1112)
ans_geese_mcmc <- geese_mcmc(
  model2fit,
  prior  = function(p) dnorm(p, sd = 2, log = TRUE), # Uniform prior
  nsteps = 20000, burnin = 0, thin = 1,
  kernel = fmcmc::kernel_am(
    warmup = 2e3,
    fixed  = c(TRUE,TRUE, rep(FALSE, nterms(model2fit) - 2)),
    lb     = -10,
    ub     = 10
))
time1 <- proc.time()
time1 - time0

traceplots(ans_geese_mcmc[,]) #, hlines = c(1,1,-1,-1,-1,-1,0))

window(ans_geese_mcmc[,-c(1,2)], start = end(ans_geese_mcmc)*4/5) |>
  apply(2, quantile, probs = c(.025, .5, .975)) |>
  t()

estimates_no_prior <- colMeans(
  window(ans_geese_mcmc, start = NSTEPS * (4/5))
  )

pred_no_prior <- predict_flock(
  p   = model2fit,
  par = estimates_no_prior,
  leave_one_out = TRUE, only_annotated = TRUE
  )

(pscores <- Map(
  \(a,b) prediction_score(a,b),
  a = lapply(pred_no_prior, \(x) do.call(rbind, x)),
  b = lapply(model_data[data_to_include], \(x) do.call(rbind, x$ann) )
  ))

saveRDS(
  list(
    mcmc = ans_geese_mcmc,
    pred = pred_no_prior,
    data = model_data[data_to_include]
  ),
  file = fn
  )

