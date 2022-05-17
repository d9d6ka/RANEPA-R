#' @title
#' Calculating long-run variance
#'
#' @description
#' Procedure ALRVR to estimate the long-run variance
#' as in Andrews (1991) and Kurozumi (2002).
#'
#' @param e (Tx1) vector or residuals.
#'
#' @return Long-run variance.
#'
#' @import MASS
alrvr <- function(e) {
    if (!is.matrix(e)) e <- as.matrix(e)

    N <- nrow(e)
    k <- 0.8
    a <- qr.solve(t(e[1:(N - 1)]) %*% e[1:(N - 1)]) %*%
        t(e[1:(N - 1)]) %*% e[2:N]
    l <- min(
        1.1447 * (4 * a^2 * N / ((1 + a)^2 * (1 - a)^2))^(1 / 3),
        1.1447 * (4 * k^2 * N / ((1 + k)^2 * (1 - k)^2))^(1 / 3) 
    )
    l <- trunc(l)
    lrv <- drop(t(e) %*% e) / N
    for (i in 1:l) {
        w <- (1 - i / (l + 1))
        lrv <- lrv + 2 * drop(t(e[1:(N - i)]) %*% e[(1 + i):N]) * w / N
    }
    return(lrv)
}
