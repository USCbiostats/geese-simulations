library(data.table)
library(aphylo)

fn <- list.files("simulations/", pattern = "simulation-study-[0-9]+\\.rds", full.names = TRUE)

dat <- lapply(fn, readRDS)

scores <- lapply(dat, \(x) {
  cbind(
    auc_geese  = x$score_geese$auc$auc,
    auc_aphylo = x$score_aphylo$auc$auc,
    mae_geese  = 1.0 - x$score_geese$obs,
    mae_aphylo = 1.0 - x$score_aphylo$obs
  )
}) |> do.call(what = rbind) |> as.data.table()


dataplot <- suppressWarnings(melt(scores))
dataplot[, method := gsub("^[a-z]+_", "", variable)]
dataplot[, stat   := gsub("_[a-z]+$", "", variable)]

library(ggplot2)
ggplot(dataplot, aes(x = value)) +
  facet_wrap(facets = ~stat) +
  geom_histogram(aes(colour = method))

ggplot(scores, aes(x = mae_geese, y = mae_aphylo)) +
  geom_point(colour = "tomato") +
  labs(x = "GEESE", y = "aphylo") +
  geom_abline(slope = 1, intercept = 0)


op <- par(mar = par()$mar * c(1,1,.5, 1))
smoothScatter(
  x = scores$mae_geese, y = scores$mae_aphylo, pch = 20,
  xlab = "GEESE", ylab = "aphylo")
abline(a = -.0025, b = 1, lwd = 4, lty = 2, col = "darkgray")
abline(a = 0, b = 1, lwd = 4, lty = 2, col = "black")
par(op)


scores[, hist(mae_geese - mae_aphylo, breaks = 100, xlim = c(-0.5, 0.5))]

ggplot(scores, aes(x = mae_geese - mae_aphylo)) +
  xlim(-.5, .5) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = 0, lwd = 1, lty = 2) +
  labs(x = expression(MAE[GEESE] - MAE[aphylo]), y = "Frequency") +
  theme(text = element_text(size = 12))
ggsave("simulations/simulation-study-analyzer.pdf", width = 4, height = 2.5)

with(scores, hist(mae_geese- mae_aphylo, breaks = 300))
with(scores, t.test(mae_geese, mae_aphylo, paired = TRUE))
with(scores, quantile(mae_geese - mae_aphylo, probs = c(.025,.975)))
scores[, quantile(mae_geese - mae_aphylo, probs = c(.025, .975))]

scores[, sd(mae_geese - mae_aphylo)]

N <- nrow(scores)
scores[, wilcox.test(mae_geese, mae_aphylo)]


scores[order(auc_geese - auc_aphylo)][,plot(auc_geese-auc_aphylo, type = "l")]
abline(h = 0)

