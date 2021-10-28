library(aphylo)

# Reading the parameter estimates
# dat <- readRDS("parameter-estimates/mcmc-joint.rds")

dat_geese_no_prior <- readRDS("parameter-estimates/mcmc-joint-geese-second-run.rds")
dat_aphylo <- readRDS("parameter-estimates/mcmc-joint-aphylo.rds")

dat <- list(
  mcmc = dat_geese_no_prior$mcmc,
  # mcmc_no_prior = dat_geese_no_prior$mcmc,
  pred = dat_geese_no_prior$pred,
  # pred_no_prior = dat_geese_no_prior$pred,
  data = dat_geese_no_prior$data,
  aphylo = dat_aphylo$aphylo_no_prior #,
  # aphylo_no_prior = dat_aphylo$aphylo_no_prior
)

# Geese scores -----------------------------------------------------------------

# Without prior
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

# Aphylo scores -----------------------------------------------------------------

# With no prior
mae_aphylo <- prediction_score(dat$aphylo, loo = TRUE)

overall_auc_aphylo <- prediction_score(
  x = do.call(rbind, lapply(mae_aphylo, "[[", "predicted")),
  expected = do.call(rbind, lapply(mae_aphylo, "[[", "expected"))
)

mae_aphylo <- sapply(mae_aphylo, function(x) 1 - x$obs)



plot(1 - mae_geese, 1 - mae_aphylo, main = "1 - MAE")
abline(a = 0, b = 1)


prop.test(table(mae_geese < mae_aphylo))
# 1-sample proportions test with continuity correction
#
# data:  table(mae_geese < mae_aphylo), null probability 0.5
# X-squared = 17.5, df = 1, p-value = 2.873e-05
# alternative hypothesis: true p is not equal to 0.5
# 95 percent confidence interval:
#   0.1516821 0.3626224
# sample estimates:
#   p
# # 0.2428571

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

op <- par(mfrow = c(3,2), mar = c(1,1,1,.5)*2.5, oma = c(2, 2, 0, 0))
for (i in 3:ncol(dat$mcmc)) {
  coda::traceplot(dat$mcmc[,i], col = "darkgray")
  # lines(lowess(dat$mcmc_no_prior[,i], f = 1/16), col = "blue", lwd = 2)
  abline(h = loc[i], col = "blue", lwd = 2, lty = "dashed")
  legend("topleft", legend = colnames(dat$mcmc)[i], bty = "n")
}
par(op)

table_ <- data.frame(
  Mean = colMeans(window(dat$mcmc[,-c(1:2)], start=10000)),
  Var=diag(cov(window(dat$mcmc[,-c(1:2)], start=10000)))
  )

knitr::kable(table_, format="latex", digits=2)

graphics.off()
svg("parameter-estimates/mcmc-joint-analysis.svg")
colors_curves <- palette.colors(n = 5, palette = "ggplot2", alpha=.9)[-1]
lwd <- 3

op <- par(bg = "white")
plot(overall_auc_geese_no_prior$auc, lwd = lwd, lty = 1, col = colors_curves[1])
with(overall_auc_aphylo_no_prior$auc, lines(x = fpr, y = tpr, lwd=lwd, lty=2, col = colors_curves[2]))
par(op)

legend(
  "bottomright",
  legend = c("MAE-AUC Method", sprintf(
    "%.2f - %.2f %s",
    c(
      1 - overall_auc_geese_no_prior$obs,
      1 - overall_auc_aphylo_no_prior$obs #,
      # 1 - dat_joint$auc_geese$obs,
      # 1 - dat_joint$auc_aphylo_joint$obs
    ),
    c(
      overall_auc_geese_no_prior$auc$auc,
      overall_auc_aphylo_no_prior$auc$auc #,
      # dat_joint$auc_geese$auc$auc,
      # dat_joint$auc_aphylo_joint$auc$auc
    ),
    c("GEESE", "aphylo") #, "GEESE (joint)", "aphylo (joint)")
  )),
  col    = c("transparent", colors_curves),
  lwd    = lwd,
  lty    = c(NA, 1:4),
  bt     = "n"
)

# legend("topleft", legend = "Including predictions\nfor TBD annotations", bty="n")
dev.off()



saveRDS(
  list(auc_geese = overall_auc_geese, auc_aphylo_joint = overall_auc_aphylo),
  file = "parameter-estimates/mcmc-joint-analysis.rds"
)
