#' @title
#' Calculating ACF values
#'
#' @description
#' A simple auxiliary function providing the estimates of autocorrelations of
#' orders from \eqn{0} to \eqn{N - 1}.
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param y An input time series of interest.
#'
#' @return The vector of ACF values of orders from 0 to \eqn{N - 1}.
#' Due to the R's way of indexing the 1-st element is the autocorrelation of
#' order 0, and the N-th value is the autcorrelation of order \eqn{N - 1}.
ACF <- function(y) {
    if (!is.matrix(y)) y <- as.matrix(y)

    n.obs <- nrow(y)

    mean.y <- mean(y)

    result <- rep(0, n.obs)

    for (i in 0:(n.obs - 1)) {
        result[i + 1] <-
            t(y[1:(n.obs - i)] - mean.y) %*% (y[(1 + i):n.obs] - mean.y) / n.obs
    }

    return(result)
}
