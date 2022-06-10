#' @title
#' Calculating ACF.
#'
#' @param y The input time series.
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
