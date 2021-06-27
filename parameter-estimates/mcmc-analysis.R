library(aphylo)
library(geese)

# Listing trees ----------------------------------------------------------------
trees <- list.files("parameter-estimates/", pattern = ".+PTHR[0-9]+.+")
trees <- gsub(".+(PTHR[0-9]+).+", "\\1", trees)
trees <- unique(trees)

# Reading the data in ----------------------------------------------------------
dat <- vector("list", length(trees))
names(dat) <- trees
for (tree in trees) {

  # Reading tree with prior
  dat_geese           <- readRDS(sprintf("parameter-estimates/mcmc-geese-%s.rds", tree))
  dat_geese_no_prior  <- readRDS(sprintf("parameter-estimates/mcmc-no-prior-%s.rds", tree))
  dat_aphylo          <- readRDS(sprintf("parameter-estimates/mcmc-aphylo-%s.rds", tree))

  # If the AUC is not present for aphylo, then need to compute it
  if (Ntrees(dat_aphylo$mcmc) > 1) {

    dat_aphylo$auc <- prediction_score(dat_aphylo$mcmc, loo = TRUE)
    dat_aphylo$auc <- lapply(
      dat_aphylo$auc, function(y) data.frame(x = unname(y$predicted), expected = unname(y$expected))
      )

    dat_aphylo$auc <- do.call(rbind, dat_aphylo$auc)
    dat_aphylo$auc <- with(dat_aphylo$auc, prediction_score(
      x = cbind(x), expected = cbind(expected))
      )
  } else {
    dat_aphylo$auc <- prediction_score(dat_aphylo$mcmc)
  }

  dat_aphylo$mae    <- 1 - dat_aphylo$auc$obs
  dat_aphylo$labels <- dat_aphylo$auc$expected
  dat_aphylo$pred   <- dat_aphylo$auc$predicted
  dat_aphylo$auc    <- dat_aphylo$auc$auc
#
#   } else {
#     dat_aphylo$labels <- with(dat_aphylo$tree, rbind(tip.annotation, node.annotation))
#   }

  # Starting point for window
  s          <- with(dat_geese, floor((end(mcmc) - start(mcmc))/2) + start(mcmc))
  s_no_prior <- with(dat_geese_no_prior, floor((end(mcmc) - start(mcmc))/2) + start(mcmc))


  # Building data to analyze
  dat[[tree]] <- list(
    mcmc_geese = colMeans(window(dat_geese$mcmc, start = s)),
    mcmc_geese_no_prior =
      colMeans(window(dat_geese_no_prior$mcmc, start = s_no_prior)),
    auc_aphylo = dat_aphylo$auc,
    auc_geese   = dat_geese$auc,
    auc_geese_no_prior  = dat_geese_no_prior$auc,
    mae_aphylo = dat_aphylo$mae,
    mae_geese  = dat_geese$mae,
    mae_geese_no_prior  = dat_geese_no_prior$mae,
    pred_geese  = dat_geese$pred,
    pred_aphylo = dat_aphylo$pred,
    lab_geese   = dat_geese$labels,
    lab_aphylo  = dat_aphylo$labels,
    nfuns       = ifelse(ncol(dat_geese$mcmc)==8, 1,2)
  )

}

# Building figures -------------------------------------------------------------

# Collective AUC
auc_joint_geese <- lapply(dat, function(d) {
  cbind(pred = d$pred_geese, expec = d$lab_geese)
})
auc_joint_geese <- do.call(rbind, auc_joint_geese)
auc_joint_geese <- prediction_score(
  cbind(auc_joint_geese[, 1]),
  cbind(auc_joint_geese[, 2])
  )

auc_joint_aphylo <- lapply(dat, function(d) {
  cbind(pred = d$pred_aphylo, expec = d$lab_aphylo)
})
auc_joint_aphylo <- do.call(rbind, auc_joint_aphylo)
auc_joint_aphylo <- prediction_score(
  cbind(auc_joint_aphylo[, 1]),
  cbind(auc_joint_aphylo[, 2])
)

# MAEs are more insteresting
maes <- lapply(dat, function(d) data.frame(
  GEESE  = d$mae_geese,
  GEESE_no_prior  = d$mae_geese_no_prior,
  aphylo = d$mae_aphylo
))
maes <- do.call(rbind, maes)

graphics.off()
pdf("parameter-estimates/mcmc-analysis-mae.pdf", width=7, height=7)
plot(
  1 - maes[,-2],
  main = "1 - Mean Absolute Error (MAE)",
  sub = "vis-à-vis comparison using individual estimates",
  pch = 19
);abline(a=0,b=1)
legend(
  "topleft",title = "Overall MAE",
  legend = sprintf(
    "GEESE: %.02f (* Sign. diff. random coin)\naphylo  : %.2f",
    1-auc_joint_geese$obs, 1-auc_joint_aphylo$obs),
  bty    = "n"
  )
dev.off()
# par(op)

# Distribution of parameter estimates ------------------------------------------

# Single fun
estimates <- lapply(dat, "[[", "mcmc_geese")
single_fun <- which(sapply(estimates, length) == 8)
estimates <- do.call(rbind, estimates[single_fun])

colnames(estimates)[c(1:3)] <- c(
  "Changes at duplication",
  "Changes (general)",
  "# of genes changing\nat duplication"
)

# Prior
nfunctions <- 1
loc <- c(
  0, 0, -1/2,
  rep(1/2, nfunctions),
  rep(-1, nfunctions),
  rep(1/2, nfunctions),
  rep(-1, nfunctions),
  rep(-9, nfunctions)
)

op <- par(mai = par("mai") * c(2, 1, 1, 1))
boxplot(estimates[,-c(1:2)], las =2)
abline(h = 0, lwd = 1.5, lty = 2)
points(x = 1:(ncol(estimates) - 2), y = loc[-c(1:2)], pch = 20, col = "tomato")
par(op)

# Double fun
estimates <- lapply(dat, "[[", "mcmc_geese")
multi_fun <- which(sapply(estimates, length) == 13)
estimates <- do.call(rbind, estimates[multi_fun])

colnames(estimates)[3] <- "# of genes changing\nat duplication"

# Prior
nfunctions <- 2
loc <- c(
  0, 0, -1/2,
  rep(1/2, nfunctions),
  rep(-1, nfunctions),
  rep(1/2, nfunctions),
  rep(-1, nfunctions),
  rep(-9, nfunctions)
)

op <- par(mai = par("mai") * c(2, 1, 1, 1))
boxplot(estimates[,-c(1:2)], las =2)
points(x = 1:(ncol(estimates) - 2), y = loc[-c(1:2)], pch = 20, col = "tomato")
abline(h = 0, lwd = 1.5, lty = 2)
par(op)

