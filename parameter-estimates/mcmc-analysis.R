# Reading the data
fn <- list.files("parameter-estimates/", full.names = TRUE, pattern = "\\.rds")
dat <- lapply(fn, readRDS)

# Extracting AUCs
aucs <- lapply(dat, function(d) {
  auc_geese  <- d$geese_auc$auc
  auc_aphylo <- d$aphylo_auc$auc
  c(geese = auc_geese, aphylo = auc_aphylo)
})

maes <- lapply(dat, function(d){

  c(geese = d$geese_mae, aphylo = d$aphylo_mae)
})

aucs <- as.data.frame(do.call(rbind, aucs))
maes <- as.data.frame(do.call(rbind, maes))

# plot(aucs$geese, aucs$aphylo)
# abline(a=0, b = 1)
table(maes[,1] <= maes[,2])
table(aucs[,1] >= aucs[,2])
