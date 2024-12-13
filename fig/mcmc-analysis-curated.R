library(aphylo)
library(geese)
source("fig/plot_functions.R")

# Reading the data -----------------------------------------------------------------------
fn <- list.files("parameter-estimates/", full.names = TRUE, pattern = "mcmc-curated-PTHR.+\\.rds")
dat <- lapply(fn, readRDS)

names(dat) <- gsub(".+(PTHR[0-9]+).+", "\\1", fn)

# Correcting AUCs aphylo
dat <- lapply(dat, function(d) {

  if (!inherits(d$aphylo_auc, "aphylo_prediction_score")) {

    d$aphylo_auc <- aphylo::prediction_score(
      x        = do.call(rbind, lapply(d$aphylo_auc, "[[", "predicted")),
      expected = do.call(rbind, lapply(d$aphylo_auc, "[[", "expected"))
    )

  }

  d

})

length(unlist(lapply(dat, function(d) as.vector(d$aphylo_auc$expected)[
  as.vector(d$aphylo_auc$expected) != 9
])))

# Extracting AUCs for the uniform prior
aucs <- lapply(dat, function(d) {
  auc_geese  <- d$geese_auc$auc$auc
  auc_aphylo <- d$aphylo_auc$auc$auc
  c(geese = auc_geese, aphylo = auc_aphylo)
})

aucs <- do.call(rbind, aucs)

# Extracting MAEs
maes <- lapply(dat, function(d) {
  auc_geese  <- 1 - d$geese_auc$obs
  auc_aphylo <- 1 - d$aphylo_auc$obs
  c(geese = auc_geese, aphylo = auc_aphylo)
})

maes <- do.call(rbind, maes)

test_methods <- function(auc., mae.) {

  # Testing the AUCs
  auc <- t.test(auc.[,1], auc.[,2], paired = TRUE, alternative = "two.sided")
  mae <- t.test(mae.[,1], mae.[,2], paired = TRUE, alternative = "two.sided")

  data.frame(
    auc = c(auc$estimate, auc$conf.int, auc$p.value),
    mae = c(mae$estimate, mae$conf.int, mae$p.value),
    row.names = c("Estimate", "lower", "upper", "p.value")
  )
}


# Extracting the AUCs for the other priors
aucs_prior <- lapply(dat, function(d) {
  auc_geese  <- d$geese_auc_prior$auc$auc
  auc_aphylo <- d$aphylo_auc_beta$auc$auc
  c(geese = auc_geese, aphylo = auc_aphylo)
})

aucs_prior <- do.call(rbind, aucs_prior)

maes_prior <- lapply(dat, function(d) {
  auc_geese  <- 1 - d$geese_auc_prior$obs
  auc_aphylo <- 1 - d$aphylo_auc_beta$obs
  c(geese = auc_geese, aphylo = auc_aphylo)
})

maes_prior <- do.call(rbind, maes_prior)

# Running the tests
test_methods(aucs, maes) |> knitr::kable()
# |         |       auc|        mae|
# |:--------|---------:|----------:|
# |Estimate | 0.1977466| -0.2050273|
# |lower    | 0.1187514| -0.2676845|
# |upper    | 0.2767417| -0.1423700|
# |p.value  | 0.0000119|  0.0000001|
test_methods(aucs_prior, maes_prior) |> knitr::kable()
# |         |        auc|        mae|
# |:--------|----------:|----------:|
# |Estimate |  0.0000940| -0.0485871|
# |lower    | -0.1038563| -0.0794822|
# |upper    |  0.1040442| -0.0176919|
# |p.value  |  0.9985474|  0.0029496|



# MCMC analysis ----------------------------------------------------------------
graphics.off()
pdf("fig/mcmc-analysis-curated-traceplots.pdf")
for (n in names(dat)) {

  traceplots(
    dat[[n]]$geese_mcmc[,-c(1,2)],
    col = adjustcolor("black", alpha.f = .5),
    smooth = TRUE
    )  

  title(n)

}
dev.off()

estimates <- lapply(dat, \(x) colMeans(window(x$geese_mcmc, start = 15000)))

# Single function
estimates_1 <- do.call(rbind, estimates[sapply(estimates, length) == 9])[,-c(1,2)]

head(estimates_1[estimates_1[,1] > estimates_1[,2],], 50)

window(dat$PTHR11575$geese_mcmc, start = 15000)[,-c(1,2)] |>
  apply(2, quantile, probs = c(.025, .975)) |>
  t()

# Two functions
estimates_2 <- do.call(rbind, estimates[sapply(estimates, length) == 14])[,-c(1,2)]

head(estimates_2[estimates_2[,1] > estimates_2[,3],], 50)

window(dat$PTHR19443$geese_mcmc, start = 15000)[,-c(1,2)] |>
  apply(2, quantile, probs = c(.025, .5, .975)) |>
  t()

traceplots(
  dat[["PTHR19443"]]$geese_mcmc[,-c(1,2)], col = adjustcolor("black", alpha.f = .5),
  smooth = TRUE
)

# Three functions
estimates_3 <- do.call(rbind, estimates[sapply(estimates, length) == 19])[,-c(1,2),drop=FALSE]

head(estimates_3, 50)

window(dat$PTHR10024$geese_mcmc, start = 15000)[,-c(1,2)] |>
  apply(2, quantile, probs = c(.025, .5, .975)) |>
  t()

traceplots(
  dat[["PTHR10024"]]$geese_mcmc[,-c(1,2)], col = adjustcolor("black", alpha.f = .5),
  smooth = TRUE
)


# Overall MAEs -----------------------------------------------------------------
geese_mae <- lapply(dat, function(d) {
  with(d$geese_auc, cbind(predicted, expected))
})

geese_mae <- do.call(rbind, geese_mae)
(geese_mae <- aphylo::prediction_score(
  x = geese_mae[,1,drop=FALSE], expected = geese_mae[,2,drop=FALSE]
  ))

aphylo_mae <- lapply(dat, function(d) {
  with(d$aphylo_auc, cbind(predicted, expected))

})

aphylo_mae <- do.call(rbind, aphylo_mae)
(aphylo_mae <- aphylo::prediction_score(
  x = aphylo_mae[,1,drop=FALSE], expected = aphylo_mae[,2,drop=FALSE]
))


# Plotting ---------------------------------------------------------------------
plot_mae(x = maes[,1], y = maes[,2], fn = "fig/mcmc-analysis-curated-mae.svg")
plot_auc(
  x = geese_mae, y = aphylo_mae,
  fn = "fig/mcmc-analysis-curated-auc.svg", width = 6, height = 4,
  title_args = list(
    sub = paste("Phylogenetic models predicting", nrow(geese_mae$predicted), "GO annotations.")
  ))

saveRDS(
  list(
    aucs   = aucs,
    maes   = maes,
    geese  = geese_mae,
    aphylo = aphylo_mae,
    raw    = dat
  ),
  file = "fig/mcmc-analysis-curated.rds"
)

