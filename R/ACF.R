#' @title
#' Calculating ACF values.
#'
#' @description
#' A simple auxiliary function providing the estimates of autocorrelations.
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param y An input time series of interest.
#'
#' @return
#' The vector of ACF values for s from 0 to \eqn{N - 1}.
ACF <- function(y) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)

    mean.y <- mean(y)

    result <- rep(0, N)

    for (i in 0:(N - 1)) {
        result[i + 1] <-
            t(y[1:(N - i)] - mean.y) %*% (y[(1 + i):N] - mean.y) / N
    }

    return(result)
}
