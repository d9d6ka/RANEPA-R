#' Calculate SSR recursively
#'
#' @details
#' Based on Brown, Durbin and Evans (1975).
#'
#' @import MASS
#' @importFrom zeallot %<-%
ssr_recursive <- function(y, x, beg, end, width = 2) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    if (beg < 1) beg <- 1
    if (nrow(y) < end) end <- nrow(y)

    N <- nrow(y)

    result <- matrix(data = Inf, nrow = N, ncol = 1)

    y_0 <- y[beg:(beg + width - 1), , drop = FALSE]
    x_0 <- x[beg:(beg + width - 1), , drop = FALSE]

    inv_xx_0 <- qr.solve(t(x_0) %*% x_0)
    c(beta_0, resid_0, predict_0, t_b_0) %<-% olsqr(y_0, x_0)

    for (r in beg:(beg + width - 2)) result[r, 1] <- 0

    result[beg + width - 1, 1] <- drop(t(resid_0) %*% resid_0)

    for (r in (beg + width):end) {
        if (r > end) break

        x_r <- x[r, , drop = FALSE]

        resid_r <- y[r, , drop = FALSE] - x_r %*% beta_0
        resid_r <- drop(resid_r)

        denom <- 1 + x_r %*% inv_xx_0 %*% t(x_r)
        denom <- drop(denom)

        inv_xx_1 <- inv_xx_0 -
            (inv_xx_0 %*% t(x_r) %*% x_r %*% inv_xx_0) / denom
        beta_1 <- beta_0 + inv_xx_1 %*% t(x_r) * resid_r
        result[r, 1] <- result[r - 1, 1] + resid_r^2 / denom

        inv_xx_0 <- inv_xx_1
        beta_0 <- beta_1
    }

    return(result)
}
