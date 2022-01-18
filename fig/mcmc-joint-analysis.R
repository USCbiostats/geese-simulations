library(aphylo)
library(geese)
source("fig/plot_functions.R")

# Reading the parameter estimates
# dat <- readRDS("parameter-estimates/mcmc-joint.rds")

dat_geese    <- readRDS("parameter-estimates/mcmc-joint-geese.rds")
dat_aphylo   <- readRDS("parameter-estimates/mcmc-joint-aphylo.rds")

# Geese scores -----------------------------------------------------------------

# With prior
pred_geese <- dat_geese$pred
pred_geese <- lapply(pred_geese, do.call, what=rbind)

dat_labs   <- lapply(dat_geese$data, "[[", "ann")
dat_labs   <- lapply(dat_labs, do.call, what=rbind)

overall_auc_geese <- prediction_score(
  x        = do.call(rbind, pred_geese),
  expected = do.call(rbind, dat_labs)
  )

mae_geese <- Map(function(p,d) prediction_score(x = cbind(p), expected = cbind(d)),
    p = pred_geese, d = dat_labs)
mae_geese <- sapply(mae_geese, function(x) 1 - x$obs)


# Aphylo scores -----------------------------------------------------------------

# With prior
mae_aphylo <- prediction_score(
  window(dat_aphylo$mcmc, start = 10000/2), loo = TRUE)

overall_auc_aphylo <- prediction_score(
  x        = do.call(rbind, lapply(mae_aphylo, "[[", "predicted")),
  expected = do.call(rbind, lapply(mae_aphylo, "[[", "expected"))
)

mae_aphylo <- sapply(mae_aphylo, function(x) 1 - x$obs)

# MCMC Trace GEESE -------------------------------------------------------------
graphics.off()
svg("fig/mcmc-joint-analysis-geese-trace.svg")
op <- par(mfrow = c(4,2), mar = c(.75,1,.5,.5)*2.5, oma = c(4, 2, 0, 0))
for (i in 3:ncol(dat_geese$mcmc)) {

  # ylim (to give space to the names)
  r    <- c(-8,8) #range(dat_aphylo$mcmc$hist[,i])
  r[2] <- r[2] + diff(r)/7

  coda::traceplot(
    window(dat_geese$mcmc[,i], start = 15000),
    col = "darkgray", ylim = r
    )
  # abline(h = loc[i], col = "blue", lwd = 2, lty = "dashed")
  legend(
    "topleft",
    legend = gsub(
      "\\s+[0-9]",
      "",
      colnames(dat_geese$mcmc)[i]), bty = "n")
}

# Adding a legend
plot.new()
# plot.window(xlim =c(0,1), ylim = c(0,1))
par(op)
title(xlab = "Step", y = "Parameter value")
dev.off()

# MCMC Trace aphylo ------------------------------------------------------------
graphics.off()
svg("fig/mcmc-joint-analysis-aphylo-trace.svg")
op <- par(mfrow = c(3,2), mar = c(.75,1,.5,.5)*2.5, oma = c(4, 2, 0, 0))
for (i in 1:ncol(dat_aphylo$mcmc$hist[[1]])) {

  # ylim (to give space to the names)
  r    <- c(0,1) #range(dat_aphylo$mcmc$hist[,i])
  r[2] <- r[2] + diff(r)/7

  coda::traceplot(dat_aphylo$mcmc$hist[,i], col = "darkgray", ylim = r)
  # abline(h = loc[i], col = "blue", lwd = 2, lty = "dashed")
  legend(
    "topleft",
    legend = gsub(
      "\\s+[0-9]",
      "",
      colnames(dat_aphylo$mcmc$hist[[1]])[i]), bty = "n")
}

# Adding a legend
plot.new()
# plot.window(xlim =c(0,1), ylim = c(0,1))
par(op)
title(xlab = "Step", y = "Parameter value")
dev.off()

# Parameter estimates ----------------------------------------------------------

# With prior
dat <- dat_geese
table_ <- data.frame(
  `Mean (Sd)` =
    sprintf(
      "%.2f (%.2f)",
      colMeans(window(dat$mcmc[,-c(1:2)], start=15000)),
      sqrt(diag(cov(window(dat$mcmc[,-c(1:2)], start=15000))))
    ), check.names = FALSE
  )

table_$`Credible Interval` <- apply(window(dat$mcmc[,-c(1:2)], start=15000), 2, function(x) {

  q <- quantile(x, probs = c(.025, .975))
  sprintf("[%.02f, %.02f]", q[1], q[2])

})

rownames(table_) <- gsub("\\s+[0-9]", "", colnames(dat$mcmc)[-c(1, 2)])

writeLines(
  knitr::kable(table_, format="latex", digits = 2, booktabs = TRUE),
  con = "fig/mcmc-joint-analysis-geese-estimates.tex"
  )

# Figures ----------------------------------------------------------------------
plot_mae(x = mae_geese, y = mae_aphylo, fn = "fig/mcmc-joint-analysis-mae.svg")
plot_auc(x = overall_auc_geese, y = overall_auc_aphylo, fn = "fig/mcmc-joint-analysis-auc.svg")

saveRDS(
  list(auc_geese = overall_auc_geese, auc_aphylo_joint = overall_auc_aphylo),
  file = "fig/mcmc-joint-analysis.rds"
)

################################################################################
# Comparing with non-pooled models
################################################################################
dat_unpooled <- readRDS("fig/mcmc-analysis-curated.rds")

used_trees <- names(dat_geese$data)
dat_unpooled <- dat_unpooled$raw[used_trees]

# Computing joint AUCs in the subset
dat_unpooled_geese <- lapply(dat_unpooled, function(d) {
  with(d$geese_auc, cbind(predicted, expected))
})

dat_unpooled_aphylo <- lapply(dat_unpooled, function(d) {
  with(d$aphylo_auc, cbind(predicted, expected))
})

dat_unpooled_geese <- do.call(rbind, dat_unpooled_geese)
prediction_score(dat_unpooled_geese[,1,drop=FALSE],  dat_unpooled_geese[,2,drop=FALSE])

dat_unpooled_aphylo <- do.call(rbind, dat_unpooled_aphylo)
prediction_score(dat_unpooled_aphylo[,1,drop=FALSE],  dat_unpooled_aphylo[,2,drop=FALSE])

# Preparing table

