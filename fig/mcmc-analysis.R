library(aphylo)
library(geese)
source("fig/plot_functions.R")

# Reading the data -----------------------------------------------------------------------
fn <- list.files("parameter-estimates/", full.names = TRUE, pattern = "mcmc-PTHR.+\\.rds")
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

prop.test(table(maes[,1] <= maes[,2]))
prop.test(table(aucs[,1] >= aucs[,2]))


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
plot_mae(x = maes[,1], y = maes[,2], fn = "fig/mcmc-analysis-mae.svg")
plot_auc(x = geese_mae, y = aphylo_mae, fn = "fig/mcmc-analysis-auc.svg")
