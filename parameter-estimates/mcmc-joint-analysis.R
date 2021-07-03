library(aphylo)

# Reading the parameter estimates
dat <- readRDS("parameter-estimates/mcmc-joint-geese3.rds")

# Retrieving predictions
pred_geese <- dat$pred
pred_geese <- lapply(pred_geese, do.call, what=rbind)

dat_labs   <- lapply(dat$data, "[[", "ann")
dat_labs   <- lapply(dat_labs, do.call, what=rbind)

overall_auc_geese <- prediction_score(
  x        = do.call(rbind, pred_geese),
  expected = do.call(rbind, dat_labs)
  )

pred_geese <- Map(function(p,d) prediction_score(x = cbind(p), expected = cbind(d)),
    p = pred_geese, d = dat_labs)
pred_geese <- sapply(pred_geese, function(x) x$obs)


pred_aphylo <- prediction_score(dat$aphylo, loo = TRUE)

overall_auc_aphylo <- prediction_score(
  x = do.call(rbind, lapply(pred_aphylo, "[[", "predicted")),
  expected = do.call(rbind, lapply(pred_aphylo, "[[", "expected"))
)

pred_aphylo <- sapply(pred_aphylo, function(x) x$obs)

plot(1 - pred_geese, 1 - pred_aphylo)
abline(a = 0, b = 1)


prop.test(table(pred_geese < pred_aphylo))
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


plot(window(dat$mcmc, start=5e3))

graphics.off()
pdf("parameter-estimates/mcmc-joint-analysis-auc.pdf")
plot(overall_auc_geese$auc, lty="dashed", lwd = 1.5)
with(overall_auc_aphylo$auc, lines(x = fpr, y = tpr, col = "red", lwd=1.5))
legend("bottomright", legend = c("geese", "aphylo"), col = c("black", "red"), lty = c(1,2), bty="n")
dev.off()

