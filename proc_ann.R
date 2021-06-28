
fns       <- list.files("parameter-estimates/", pattern = "-geese-PTHR", full.names = TRUE)
geese_dat <- lapply(fns, readRDS)

names(geese_dat) <- gsub(".+(PTHR[0-9]+).+", "\\1", fns)

geese_dat <-lapply(geese_dat, function(d) {

  # Capturing gene names
  if (aphylo::Ntrees(d$tree) == 1) {
    gnames <- with(d$tree$tree, c(tip.label, node.label))
    ord    <- c(d$tree$tree$edge[,2], ape::Ntip(d$tree$tree) + 1)
  } else {
    gnames <- with(d$tree[[1]]$tree, c(tip.label, node.label))
    ord    <- c(d$tree[[1]]$tree$edge[,2], ape::Ntip(d$tree[[1]]$tree) + 1)
  }

  gnames <- gnames[ord]

  data.frame(
    gene       = gnames,
    pred_geese = d$pred,
    lab_geese  = d$labels
  )

})

geese_dat <- lapply(geese_dat, function(d) {
  d[d$lab_geese != 9L,,drop=FALSE]
})

geese_dat <- cbind(
  tree = rep(names(geese_dat), sapply(geese_dat, nrow)),
  do.call(rbind, geese_dat)
)

aphylo::prediction_score(
  x = cbind(geese_dat$pred_geese),
  expected = cbind(geese_dat$lab_geese)
)

# PTHR10082.457  PTHR10082                 AN11 0.75780565         1
# PTHR10082.458  PTHR10082                 AN12 0.60182882         1
# PTHR10082.465  PTHR10082                 AN19 0.75806413         1
# PTHR10082.466  PTHR10082                 AN20 0.60537293         1
# PTHR10082.677  PTHR10082     UniProtKB=G3S6H4 0.75780565         0
# PTHR10082.678  PTHR10082     UniProtKB=F7GE04 0.60182882         0
# PTHR10082.687  PTHR10082     UniProtKB=F6TXW5 0.75791001         0
# PTHR10082.688  PTHR10082                AN242 0.60363680         0
# summary(geese_dat$PTHR10082$tree)
