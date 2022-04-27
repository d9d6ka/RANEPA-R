# Carrion-i-Silvestre and Sansó 2006 OBES

determi <- function(model, t, tb) {
    du <- rbind(
        matrix(data = 0, nrow = tb, ncol = 1),
        matrix(data = 1, nrow = t - tb, ncol = 1)
    )
    dt <- rbind(
        matrix(data = 0, nrow = tb, ncol = 1),
        matrix(data = 1:(t - tb), nrow = t - tb, ncol = 1)
    )

    if (model == 1) {
        xt <- cbind(
            matrix(data = 1, nrow = t, ncol = 1),
            du
        )
    }
    else if (model == 2) {
        xt <- cbind(
            matrix(data = 1, nrow = t, ncol = 1),
            du,
            matrix(data = 1:t, nrow = t, ncol = 1)
        )
    }
    else if (model == 3) {
        xt <- cbind(
            matrix(data = 1, nrow = t, ncol = 1),
            matrix(data = 1:t, nrow = t, ncol = 1),
            dt
        )
    }
    else if (model == 4) {
        xt <- cbind(
            matrix(data = 1, nrow = t, ncol = 1),
            du,
            matrix(data = 1:t, nrow = t, ncol = 1),
            dt
        )
    }
    else {
        stop("ERROR: Try to specify the deterministic component again")
    }

    return(xt)
}

alrvr <- function(e) {
    t <- length(e)
    k <- 0.8
    a <- solve(e[1:(t - 1)] %*% e[1:(t - 1)]) %*% e[1:(t - 1)] %*% e[2:t]
    l <- min(
        1.1447 * (4 * a^2 * t / ((1 + a)^2 * (1 - a)^2))^(1 / 3), # nolint
        1.1447 * (4 * k^2 * t / ((1 + k)^2 * (1 - k)^2))^(1 / 3)  # nolint
    )
    l <- trunc(l)
    lrv <- as.numeric((e %*% e) / t)
    i <- 1
    while (i <= l) {
        w <- (1 - i / (l + 1))
        lrv <- lrv + 2 * as.numeric(e[1:(t - i)] %*% e[(1 + i):t]) * w / t
        i <- i + 1
    }
    return(lrv)
}
