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


#' @importFrom zeallot %<-%
h0W <- function(y) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)

    c(beta, ., ., .) %<-% OLS(y[2:N, ], y[1:(N-1), ])
    beta <- drop(beta)

    a <- (4 * beta^2) / ((1 - beta)^4)
    R <- ACF(y)

    lambda <- matrix(0, nrow = N - 1, ncol = 1)
    s <- as.matrix(1:(N - 1))

    m <- 1.3221 * (a * N)^(1 / 5)
    delta <- (6 * pi * s) / (5 * m)

    for (i in 1:(N - 1)) {
        lambda[i] <-
            3 * (sin(delta[i]) / delta[i] - cos(delta[i])) / (delta[i]^2)
    }

    return(
        list(
            h0 = drop(R[1] + 2 * t(lambda) %*% R[2:N]),
            m = m
        )
    )
}
