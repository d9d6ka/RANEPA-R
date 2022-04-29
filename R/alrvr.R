## Procedure ALRVR to estimate the long-run variance as
## in Andrews (1991) and Kurozumi (2002).
alrvr <- function(e) {
    if (!is.matrix(e)) e <- as.matrix(e)

    t <- nrow(e)
    k <- 0.8
    a <- solve(t(e[1:(t - 1)]) %*% e[1:(t - 1)]) %*% t(e[1:(t - 1)]) %*% e[2:t]
    l <- min(
        1.1447 * (4 * a^2 * t / ((1 + a)^2 * (1 - a)^2))^(1 / 3), # nolint
        1.1447 * (4 * k^2 * t / ((1 + k)^2 * (1 - k)^2))^(1 / 3)  # nolint
    )
    l <- trunc(l)
    lrv <- as.numeric(t(e) %*% e) / t
    for (i in 1:l) {
        w <- (1 - i / (l + 1))
        lrv <- lrv + 2 * as.numeric(t(e[1:(t - i)]) %*% e[(1 + i):t]) * w / t
    }
    return(lrv)
}
