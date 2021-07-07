library(aphylo)

# Reading the parameter estimates
dat <- readRDS("parameter-estimates/mcmc-joint.rds")

# Retrieving predictions
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


mae_aphylo <- prediction_score(dat$aphylo, loo = TRUE)

overall_auc_aphylo <- prediction_score(
  x = do.call(rbind, lapply(mae_aphylo, "[[", "predicted")),
  expected = do.call(rbind, lapply(mae_aphylo, "[[", "expected"))
)

mae_aphylo <- sapply(mae_aphylo, function(x) x$obs)

plot(1 - mae_geese, 1 - mae_aphylo, main = "1 - MAE")
abline(a = 0, b = 1)


prop.test(table(mae_geese < mae_aphylo))
#
# 1-sample proportions test with continuity correction
#
# data:  table(pred_geese < pred_aphylo), null probability 0.5
# X-squared = 1.1571, df = 1, p-value = 0.2821
# alternative hypothesis: true p is not equal to 0.5
# 95 percent confidence interval:
#   0.3128160 0.5522044
# sample estimates:
#   p
# 0.4285714

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

# Preparing for ggplot
ggdat <- lapply(3:ncol(dat$mcmc), function(i)
  data.frame(
    Term  = colnames(dat$mcmc)[i],
    Step  = 1:nrow(dat$mcmc),
    Param = as.vector(dat$mcmc[,i]),
    Prior = loc[i]
  )
)
ggdat <- do.call(rbind, ggdat)
library(ggplot2)
ggplot(data = subset(ggdat, Step >= 5000), mapping = aes(y = Param, x = Step)) +
  geom_line(col = "tomato") +
  facet_wrap(vars(Term)) +
  geom_abline(slope = 0, intercept = 0, lty=2) +
  geom_line(aes(x = Step, y = Prior))

op <- par(mfrow = c(3,2), mar = c(1,1,1,.5)*2.5, oma = c(2, 2, 0, 0))
for (i in 3:ncol(dat$mcmc)) {
  coda::traceplot(dat$mcmc[,i], col = "gray")
  lines(lowess(dat$mcmc[,i], f = 1/16), col = "blue", lwd = 2)
  abline(h = loc[i], col = "red", lwd = 2, lty = "dashed")
  legend("topleft", legend = colnames(dat$mcmc)[i], bty = "n")
}
par(op)

# graphics.off()
# pdf("parameter-estimates/mcmc-joint-analysis-auc.pdf")
# plot(overall_auc_geese$auc, lty=1, lwd = 1.5)
# with(overall_auc_aphylo$auc, lines(x = fpr, y = tpr, col = "red", lwd=1.5, lty=2))
# legend("bottomright", legend = c("geese", "aphylo"), col = c("black", "red"), lty = c(1,2), bty="n", lwd=1.5)
# dev.off()

saveRDS(
  list(auc_geese = overall_auc_geese, auc_aphylo_joint = overall_auc_aphylo),
  file = "parameter-estimates/mcmc-joint-analysis.rds"
)
