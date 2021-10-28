# Reading the MCMC estimates
x <- readRDS("inst/simulation-study/simulation-study.rds")

# Testing
params <- c(
  # Gaining function 0, and 0->1
  2.5, 1.5,
  # Overall changes
  1.5, -2,
  # Root probabilities
  -10, -10
)

names(params) <- c(
  "gain0", "neofun01", "changes_dpl", "changes_sp",
  "root0", "root1"
)


# Checking coverage
estimates <- lapply(x, "[[", "estimates")
coverage <- sapply(estimates, function(e) {
  (params >= e[1,]) & (params <= e[3,])
})

colMeans(t(coverage)) # Almost 99 pcent

# Checking bias
bias <- sapply(estimates, function(e) {
  e[2,]
})

rownames(bias) <- c(
  "gain0", "neofun01", "changes_dpl", "changes_sp",
  "root0", "root1"
)

# Pretty plot
graphics.off()
svg("figures.svg")
boxplot(t(bias), border = "darkgray", col = "steelblue", outline = FALSE,
        ylab = "Parameter Value", xlab = "Parameter",
        ylim = c(-10,6))
# grid()
abline(h = 0, lty = 2, lwd = 2, col = "gray")
points(
  x = 1:length(params), lwd = 2,
  y = params, cex = 2, pch=23, col = "black", bg = "tomato"
)
text(
  x = 1:length(params) + .3,
  y = params + .15,
  labels = sprintf("%.2f", params)
)
dev.off()
