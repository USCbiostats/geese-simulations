#' @param x,y Objects of class `aphylo_prediction_score`.
#' @param fn Name of the file where to save the svg figure
plot_auc <- function(x, y, fn, title_args = list(), ...) {

  graphics.off()
  svg(fn, ...)
  colors_curves <- palette.colors(n = 5, palette = "ggplot2", alpha=.9)[-1]
  lwd <- 3

  op <- par(bg = "white")
  plot(x$auc, lwd = lwd, lty = 1, col = colors_curves[1])
  with(y$auc, lines(x = fpr, y = tpr, lwd=lwd, lty=2, col = colors_curves[2]))
  par(op)

  legend(
    "bottomright",
    legend = c("MAE-AUC Method", sprintf(
      "%.2f - %.2f %s",
      c(
        1 - x$obs,
        1 - y$obs
      ),
      c(
        x$auc$auc,
        y$auc$auc
      ),
      c("GEESE", "aphylo")
    )),
    col    = c("transparent", colors_curves),
    lwd    = lwd,
    lty    = c(NA, 1:2),
    bt     = "n"
  )

  if (length(title_args))
    do.call(graphics::title, title_args)

  dev.off()
}

# Maes one-at-a-time
#' @param x,y Numeric vectors of MAEs.
plot_mae <- function(x, y, fn, xlab. = "GEESE", ylab. = "aphylo", main. = "MAE per tree", size = 7, cex. = 1.75) {
  graphics.off()
  svg(fn, width = size, height = size)
  ranges <- apply(cbind(x, y), 2, range)
  op <- par(bg = "lightgray", cex = cex., cex.lab = cex., cex.axis = cex.,
            cex.main = cex., no.readonly = TRUE)
  plot.new()
  plot.window(xlim = ranges[,1] * 1.2, ylim = ranges[,2] * 1.2)
  rect(-1,-1,2,2, col = "white", border = "transparent")
  points(x = x, y = y, pch=19, col = adjustcolor("black", alpha = .7))
  abline(a=0, b = 1, lwd = 2, lty = "dashed", col = "darkgray")
  axis(1)
  axis(2)
  title(xlab = xlab., ylab = ylab., main = main.)
  par(op)
  dev.off()
}


traceplots <- function(dat, ylim = c(-10,10), hlines = NULL, ...) {
  k <- ncol(dat)

  op <- par(mfrow = c(ceiling(k/3), 3), mar = c(2.5, 2.5, .5, .5))
  for (i in 1:ncol(dat)) {
    coda::traceplot(dat[,i,drop=FALSE], ylim = ylim, main = "", ...)
    abline(h = 0, lwd = 1.5, lty = 2, col = "tomato")
    legend("topleft", legend = colnames(dat)[i], bt = "n")
    if (length(hlines)) {
      abline(h = hlines[i], lty = 3, col = "steelblue")
      text(x = 0, y = hlines[i], labels = "prior")
    }

  }

  par(op)

}
