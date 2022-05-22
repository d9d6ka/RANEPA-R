#' @title
#' reindex - function that makes reindexing.
#'
#' @details
#' See Cavaliere and Taylor (2007, 2008) and
#' Kurozumi, Skrobotov and Tsaryov (2021) for more information.
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
