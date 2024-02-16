library(aphylo)
library(geese)
library(coda)

# Fitting partially annotated trees --------------------------------------------
model_data <- readRDS("data/model_data.rds")
NSTEPS     <- 2e4

# Subsetting trees to obtain parameter estimates
data_to_include <- which(
  (sapply(model_data, "[[", "max_size") < 1e8) &
  # (sapply(model_data, "[[", "polytomies") <= 25) &
    (sapply(model_data, "[[", "nfuns") == 1)
  )
  
data_to_include <- names(model_data)[data_to_include]

# Aphylo
set.seed(212)
adata <- do.call(c, lapply(model_data[data_to_include], "[[", "tree"))
ans_aphylo <- aphylo_mcmc(
  adata ~ mu_d + mu_s + Pi,
  priors = uprior(),
)

auc_aphylo <- prediction_score(ans_aphylo, loo = TRUE)

saveRDS(
  list(
    mcmc = ans_aphylo,
    pred = auc_aphylo,
    data = model_data[data_to_include]
  ),
  file = "parameter-estimates/mcmc-joint-aphylo.rds"
)




