library(aphylo)
library(geese)
source("fig/plot_functions.R")

# Reading the data -----------------------------------------------------------------------
fn <- list.files("parameter-estimates/", full.names = TRUE, pattern = "mcmc-unif-prior-curated-PTHR.+\\.rds")
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

# Extracting AUCs
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

tryCatch(prop.test(table(maes[,1] <= maes[,2])), error = function(e) e)
tryCatch(prop.test(table(aucs[,1] >= aucs[,2])), error = function(e) e)


nfunctions <- aphylo::Nann(dat[[2]]$tree)
loc <- c(
  # Overall changes
  0, 0,
  # Genes changing at duplication
  -1/2,
  # Gains and loss x nfunctions
  rep(c(1/2, -1/2), nfunctions * 2),
  rep(0, nfunctions)
  )

# MCMC analysis ----------------------------------------------------------------
graphics.off()
pdf("fig/mcmc-analysis-curated-traceplots.pdf")
for (n in names(dat)) {
  traceplots(
    dat[[n]]$geese_mcmc[,-c(1,2)], col = adjustcolor("black", alpha.f = .5),
    smooth = TRUE
    )
  title(n)
}
dev.off()

estimates <- lapply(dat, \(x) colMeans(window(x$geese_mcmc, start = 15000)))

# Single function
estimates_1 <- do.call(rbind, estimates[sapply(estimates, length) == 9])[,-c(1,2)]
View(estimates_1[estimates_1[,1] > estimates_1[,2],])

window(dat$PTHR11575$geese_mcmc, start = 15000)[,-c(1,2)] |>
  apply(2, quantile, probs = c(.025, .975)) |>
  t()

# Two functions
estimates_2 <- do.call(rbind, estimates[sapply(estimates, length) == 14])[,-c(1,2)]
View(estimates_2[estimates_2[,1] > estimates_2[,3],])

window(dat$PTHR19443$geese_mcmc, start = 15000)[,-c(1,2)] |>
  apply(2, quantile, probs = c(.025, .5, .975)) |>
  t()

traceplots(
  dat[["PTHR19443"]]$geese_mcmc[,-c(1,2)], col = adjustcolor("black", alpha.f = .5),
  smooth = TRUE
)

# Three functions
estimates_3 <- do.call(rbind, estimates[sapply(estimates, length) == 19])[,-c(1,2),drop=FALSE]
View(estimates_3)

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
plot_mae(x = maes[,1], y = maes[,2], fn = "fig/mcmc-analysis-unif-prior-curated-mae.svg")
plot_auc(
  x = geese_mae, y = aphylo_mae,
  fn = "fig/mcmc-analysis-unif-prior-curated-auc.svg", width = 6, height = 4,
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

