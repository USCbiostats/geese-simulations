fns <- list.files(
  "parameter-estimates/", pattern = "mcmc-aphylo-PTHR[0-9]+\\.rds$",
  full.names = TRUE
  )

for (f in fns) {
  d <- readRDS(f)

  # Retrieving the aphylo elements only
  output <- with(d, {
    list(
      mcmc   = aphylo_mcmc,
      auc    = aphylo_auc,
      mae    = aphylo_mae,
      pred   = aphylo_pred,
      labels = labels,
      tree   = tree
    )

  })

  # fn <- gsub("mcmc-", "mcmc-aphylo-", f)

  out <- tryCatch(saveRDS(output, file = f), error = function(e) e)
  # if (inherits(out, "error")) {
  #   message("File ", f, " failed.")
  # } else {
    # file.remove(f)
    message("File ", f, " OK.")
  # }



}
