#' @title
#' reindex - function that makes reindexing.
#'
#' @param u The residuals series for reindexing.
#'
#' @references
#' Cavaliere, Giuseppe, and A. M. Robert Taylor.
#' “Time-Transformed Unit Root Tests for Models with Non-Stationary Volatility.”
#' Journal of Time Series Analysis 29, no. 2 (March 2008): 300–330.
#' https://doi.org/10.1111/j.1467-9892.2007.00557.x.
#'
#' Kurozumi, Eiji, Anton Skrobotov, and Alexey Tsarev.
#' “Time-Transformed Test for the Explosive Bubbles under
#' Non-Stationary Volatility.”
#' arXiv, November 15, 2021. http://arxiv.org/abs/2012.13937.
reindex <- function(u) {
    N <- length(u)
    u.2 <- as.numeric(u^2)

    s <- (0:N) / N

    eta.hat <- rep(0, (N + 1))
    for (i in 2:(N + 1)) {
        sT <- floor(s[i] * N)
        eta.hat[i] <- (sum(u.2[1:sT]) + (s[i] * N - sT) * u.2[sT + 1]) /
            sum(u.2)
    }
    eta.hat[N + 1] <- 1

    eta.hat.inv <- rep(0, (N + 1))
    for (i in 2:(N + 1)) {
        k <- length(eta.hat[eta.hat < s[i]])
        s0 <- (k - 1) / N
        eta.hat.inv[i] <- s0 +
            (s[i] - eta.hat[k]) / (eta.hat[k + 1] - eta.hat[k]) / N
    }
    eta.hat.inv[N + 1] <- 1

    new.index <- floor(eta.hat.inv * N)

    return(
        list(
            u = u,
            s = s,
            eta.hat = eta.hat,
            eta.hat.inv = eta.hat.inv,
            new.index = new.index
        )
    )
}
