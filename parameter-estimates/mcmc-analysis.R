library(aphylo)
library(geese)

# Reading the data -----------------------------------------------------------------------
fn <- list.files("parameter-estimates/", full.names = TRUE, pattern = "mcmc-PTHR.+\\.rds")
dat <- lapply(fn, readRDS)

# Joint estimates
dat_joint <- readRDS("parameter-estimates/mcmc-joint-analysis.rds")

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

# plot(aucs$geese, aucs$aphylo)
# abline(a=0, b = 1)
table(maes[,1] <= maes[,2])
table(aucs[,1] >= aucs[,2])

plot.new()
plot.window(xlim = c(0,1), ylim = c(0,1))
axis(1)
axis(2)
polygon(
  x = c(0,.5,.5,0), y = c(0, 0, .5,.5),
  density = 20,
  col = adjustcolor("gray", .9))
points(1 - maes)
abline(a=0,b=1)

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

# op <- par(mfrow = c(4,2), mar = c(1,1,1,1)*2.5)
# # coda::traceplot(window(dat[[2]]$geese_mcmc, start = 6000),ylim=c(-10,10))
# for (i in 1:ncol(dat[[2]]$geese_mcmc)) {
#   coda::traceplot(dat[[2]]$geese_mcmc[,i], col = "gray")
#   lines(lowess(dat[[2]]$geese_mcmc[,i], f = 1/16), col = "blue", lwd = 2)
#   abline(h = loc[i], col = "red", lwd = 2, lty = "dashed")
#   legend("topleft", legend = colnames(dat[[2]]$geese_mcmc)[i], bty = "n")
# }
#
# par(op)

# Overall MAEs -----------------------------------------------------------------
geese_mae <- lapply(dat, function(d) {
  with(d$geese_auc, cbind(predicted, expected))
})

geese_mae <- do.call(rbind, geese_mae)
geese_mae <- aphylo::prediction_score(
  x = geese_mae[,1,drop=FALSE], expected = geese_mae[,2,drop=FALSE]
  )

aphylo_mae <- lapply(dat, function(d) {
  with(d$aphylo_auc, cbind(predicted, expected))

})

aphylo_mae <- do.call(rbind, aphylo_mae)
aphylo_mae <- aphylo::prediction_score(
  x = aphylo_mae[,1,drop=FALSE], expected = aphylo_mae[,2,drop=FALSE]
)

graphics.off()
svg("parameter-estimates/mcmc-analysis.svg")
colors_curves <- palette.colors(n = 5, palette = "ggplot2", alpha=.9)[-1]
lwd <- 3

op <- par(bg = "white")
plot(geese_mae$auc, lwd = lwd, lty = 1, col = colors_curves[1])
with(aphylo_mae$auc, lines(x = fpr, y = tpr, lwd=lwd, lty=2, col = colors_curves[2]))
with(dat_joint$auc_geese$auc, lines(x = fpr, y = tpr, lwd=lwd, lty=3, col = colors_curves[3]))
with(dat_joint$auc_aphylo_joint$auc, lines(x = fpr, y = tpr, lwd=lwd, lty=4, col = colors_curves[4]))
par(op)

legend(
  "bottomright",
  legend = c("MAE-AUC Method", sprintf(
    "%.2f - %.2f %s",
    c(
      1 - geese_mae$obs,
      1 - aphylo_mae$obs,
      1 - dat_joint$auc_geese$obs,
      1 - dat_joint$auc_aphylo_joint$obs
      ),
    c(
      geese_mae$auc$auc,
      aphylo_mae$auc$auc,
      dat_joint$auc_geese$auc$auc,
      dat_joint$auc_aphylo_joint$auc$auc
    ),
    c("GEESE", "aphylo", "GEESE (joint)", "aphylo (joint)")
    )),
  col    = c("transparent", colors_curves),
  lwd    = lwd,
  lty    = c(NA, 1:4),
  bt     = "n"
)

legend("topleft", legend = "Including predictions\nfor TBD annotations", bty="n")
dev.off()

