#' @title
#' Calculate SSR recursively
#'
#' @details
#' Based on Brown, Durbin and Evans (1975).
#'
#' @param y (Tx1)-vector of the dependent variable.
#' @param x (Txk)-vector of the explanatory stochastic regressors.
#' @param beg The start of SSR calculating period.
#' @param end The end of SSR calculating period.
#' @param width Minimum spacing between the breaks.
#'
#' @return The vector of calculated recursive SSR.
#'
#' @import MASS
#' @importFrom zeallot %<-%
ssr.recursive <- function(y, x, beg, end, width = 2) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    if (beg < 1) beg <- 1
    if (nrow(y) < end) end <- nrow(y)

    N <- nrow(y)

    result <- matrix(data = Inf, nrow = N, ncol = 1)

    y.0 <- y[beg:(beg + width - 1), , drop = FALSE]
    x.0 <- x[beg:(beg + width - 1), , drop = FALSE]

    inv.XX.0 <- qr.solve(t(x.0) %*% x.0)
    c(beta.0, resid.0, ., .) %<-% olsqr(y.0, x.0)

    for (step in beg:(beg + width - 2)) result[step, 1] <- 0

    result[beg + width - 1, 1] <- drop(t(resid.0) %*% resid.0)

    for (step in (beg + width):end) {
        if (step > end) break

        step.x <- x[step, , drop = FALSE]

        step.resid <- y[step, , drop = FALSE] - step.x %*% beta.0
        step.resid <- drop(step.resid)

        denom <- 1 + step.x %*% inv.XX.0 %*% t(step.x)
        denom <- drop(denom)

        inv.XX.1 <- inv.XX.0 -
            (inv.XX.0 %*% t(step.x) %*% step.x %*% inv.XX.0) / denom

        beta.1 <- beta.0 + inv.XX.1 %*% t(step.x) * step.resid

        result[step, 1] <- result[step - 1, 1] + step.resid^2 / denom

        inv.XX.0 <- inv.XX.1

        beta.0 <- beta.1
    }

    return(result)
}
