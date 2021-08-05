library(aphylo)
library(geese)
source("fig/plot_functions.R")

# Reading the parameter estimates
# dat <- readRDS("parameter-estimates/mcmc-joint.rds")

dat_geese_prior    <- readRDS("parameter-estimates/mcmc-joint-geese-prior.rds")
dat_geese_no_prior <- readRDS("parameter-estimates/mcmc-joint-geese-no-prior.rds")
dat_aphylo         <- readRDS("parameter-estimates/mcmc-joint-aphylo.rds")

dat <- list(
  mcmc            = dat_geese_prior$mcmc,
  mcmc_no_prior   = dat_geese_no_prior$mcmc,
  pred            = dat_geese_prior$pred,
  pred_no_prior   = dat_geese_no_prior$pred,
  data            = dat_geese_prior$data,
  aphylo          = dat_aphylo$aphylo,
  aphylo_no_prior = dat_aphylo$aphylo_no_prior
)

# Geese scores -----------------------------------------------------------------

# With prior
pred_geese <- dat$pred
pred_geese <- lapply(pred_geese, do.call, what=rbind)

dat_labs   <- lapply(dat$data, "[[", "ann")
dat_labs   <- lapply(dat_labs, do.call, what=rbind)

overall_auc_geese <- prediction_score(
  x        = do.call(rbind, pred_geese),
  expected = do.call(rbind, dat_labs)
  )

mae_geese <- Map(function(p,d) prediction_score(x = cbind(p), expected = cbind(d)),
    p = pred_geese, d = dat_labs)
mae_geese <- sapply(mae_geese, function(x) 1 - x$obs)

# Without prior
pred_geese_no_prior <- dat$pred_no_prior
pred_geese_no_prior <- lapply(pred_geese_no_prior, do.call, what=rbind)

overall_auc_geese_no_prior <- prediction_score(
  x        = do.call(rbind, pred_geese_no_prior),
  expected = do.call(rbind, dat_labs)
)

mae_geese_no_prior <- Map(function(p,d) prediction_score(x = cbind(p), expected = cbind(d)),
                 p = pred_geese_no_prior, d = dat_labs)
mae_geese_no_prior <- sapply(mae_geese_no_prior, function(x) 1 - x$obs)

# Aphylo scores -----------------------------------------------------------------

# With prior
mae_aphylo <- prediction_score(dat$aphylo, loo = TRUE)

overall_auc_aphylo <- prediction_score(
  x = do.call(rbind, lapply(mae_aphylo, "[[", "predicted")),
  expected = do.call(rbind, lapply(mae_aphylo, "[[", "expected"))
)

mae_aphylo <- sapply(mae_aphylo, function(x) 1 - x$obs)

# Without prior
mae_aphylo_no_prior <- prediction_score(dat$aphylo_no_prior, loo = TRUE)

overall_auc_aphylo_no_prior <- prediction_score(
  x = do.call(rbind, lapply(mae_aphylo_no_prior, "[[", "predicted")),
  expected = do.call(rbind, lapply(mae_aphylo_no_prior, "[[", "expected"))
)

mae_aphylo_no_prior <- sapply(mae_aphylo_no_prior, function(x) 1 - x$obs)

# Overal MAE and AUC
nfunctions <- 1
loc <- c(
  # Overall changes
  0, 0,
  # Genes changing at duplication
  -1/2,
  # Gains and loss x nfunctions (duplication)
  rep(1/2, nfunctions), rep(-1/2, nfunctions),
  # Gains and loss x nfunctions (speciation)
  rep(-1/2, nfunctions), rep(-1/2, nfunctions),
  rep(0, nfunctions)
)

# MCMC Trace -------------------------------------------------------------------
graphics.off()
svg("fig/mcmc-joint-analysis-trace.svg")
op <- par(mfrow = c(3,2), mar = c(.75,1,.5,.5)*2.5, oma = c(4, 2, 0, 0))
for (i in 3:ncol(dat$mcmc)) {

  # ylim (to give space to the names)
  r <- range(dat$mcmc[,i])
  r[2] <- r[2] + diff(r)/7

  coda::traceplot(dat$mcmc[,i], col = "darkgray", ylim = r)
  abline(h = loc[i], col = "blue", lwd = 2, lty = "dashed")
  legend(
    "topleft",
    legend = gsub("\\s+[0-9]", "", colnames(dat$mcmc)[i]), bty = "n")
}

# Adding a legend
plot.new()
plot.window(xlim =c(0,1), ylim = c(0,1))
legend(
  "center",
  legend = c("Prior value", "MCMC Trace"), col = c("blue", "darkgray"),
  lty = c(2, 1), lwd = 2, cex = 1.5, bty = "n"
  )
par(op)
title(xlab = "Step", y = "Parameter value")
dev.off()

# Parameter estimates ----------------------------------------------------------

# With prior
table_ <- data.frame(
  `Mean (Sd)` =
    sprintf(
      "%.2f (%.2f)",
      colMeans(window(dat$mcmc[,-c(1:2)], start=10000)),
      sqrt(diag(cov(window(dat$mcmc[,-c(1:2)], start=10000))))
    ), check.names = FALSE
  )

table_$`Credible Interval` <- apply(window(dat$mcmc[,-c(1:2)], start=10000), 2, function(x) {

  q <- quantile(x, probs = c(.025, .975))
  sprintf("[%.02f, %.02f]", q[1], q[2])

})

rownames(table_) <- gsub("\\s+[0-9]", "", colnames(dat$mcmc_no_prior)[-c(1, 2)])

writeLines(
  knitr::kable(table_, format="latex", digits = 2, booktabs = TRUE),
  con = "fig/mcmc-joint-analysis-estimates.tex"
  )

# Without prior
table_ <- data.frame(
  `Mean (Sd)` =
    sprintf(
      "%.2f (%.2f)",
      colMeans(window(dat$mcmc_no_prior[,-c(1:2)], start=10000)),
      sqrt(diag(cov(window(dat$mcmc_no_prior[,-c(1:2)], start=10000))))
    ),
  check.names = FALSE
)

table_$`Credible Interval` <- apply(window(dat$mcmc_no_prior[,-c(1:2)], start=10000), 2, function(x) {

  q <- quantile(x, probs = c(.025, .975))
  sprintf("[%.02f, %.02f]", q[1], q[2])

})

rownames(table_) <- gsub("\\s+[0-9]", "", colnames(dat$mcmc_no_prior)[-c(1, 2)])

writeLines(
  knitr::kable(table_, format="latex", digits=2, booktabs = TRUE),
  con = "fig/mcmc-joint-analysis-estimates-no-prior.tex"
  )

# Figures ----------------------------------------------------------------------
plot_mae(x = mae_geese, y = mae_aphylo, fn = "fig/mcmc-joint-analysis-mae.svg")
plot_auc(x = overall_auc_geese, y = overall_auc_aphylo, fn = "fig/mcmc-joint-analysis-auc.svg")

saveRDS(
  list(auc_geese = overall_auc_geese, auc_aphylo_joint = overall_auc_aphylo),
  file = "parameter-estimates/mcmc-joint-analysis.rds"
)
