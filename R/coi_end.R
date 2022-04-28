# Carrion-i-Silvestre and Sansó 2006 OBES

require(zeallot)

## `[` <- function(...) base::`[`(..., drop = FALSE)

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
    a <- solve(t(e[1:(t - 1)]) %*% e[1:(t - 1)]) %*% t(e[1:(t - 1)]) %*% e[2:t]
    l <- min(
        1.1447 * (4 * a^2 * t / ((1 + a)^2 * (1 - a)^2))^(1 / 3), # nolint
        1.1447 * (4 * k^2 * t / ((1 + k)^2 * (1 - k)^2))^(1 / 3)  # nolint
    )
    l <- trunc(l)
    lrv <- as.numeric((t(e) %*% e) / t)
    for (i in 1:l) {
        w <- (1 - i / (l + 1))
        lrv <- lrv + 2 * as.numeric(t(e[1:(t - 1)]) %*% e[(1 + i):t]) * w / t
    }
    return(lrv)
}

lagn <- function(x, i) {
    if (!is.matrix(x)) x <- as.matrix(x)
    t <- nrow(x)
    k <- ncol(x)
    if (i > 0)
        return(
            rbind(
                matrix(data = NA, nrow = i, ncol = k),
                x[1:(t - i), , drop = FALSE]
            )
        )
    else
        return(
            rbind(
                x[(1 + abs(i)):t, , drop = FALSE],
                matrix(data = NA, nrow = abs(i), ncol = k)
            )
        )
}

DOLS <- function(y, x, model, klags, kleads, tb) {
    t <- length(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    d_streg <- x[2:t, , drop = FALSE] - x[1:(t - 1), , drop = FALSE]

    for (i in 1:klags) {
        d_streg <- cbind(
            d_streg,
            lagn(d_streg, i)
        )
    }

    d_streg_r <- d_streg[(t-1):1, , drop = FALSE]

    for (i in 1:kleads) {
        d_streg_r <- cbind(
            d_streg_r,
            lagn(d_streg_r, i)
        )
    }

    if (klags != 0 & kleads != 0) {
        lags <- d_streg
        leads <- d_streg_r[(t - 1):1, (ncol(streg) + 1):(ncol(d_streg_r)), drop = FALSE]
        ll <- cbind(lags, leads)[(klags + 1):(t - 1 - kleads), , drop = FALSE]
    }
    else if (klags != 0 & kleads == 0) {
        lags <- d_streg
        ll <- lags[(klags + 1):(t - 1), , drop = FALSE]
    }
    else if (klags == 0 & kleads != 0) {
        lags <- d_streg
        leads <- d_streg_r[(t - 1):1, (ncol(streg) + 1):(ncol(d_streg_r)), drop = FALSE]
        ll <- cbind(lags, leads)[1:(t - 1 - kleads), , drop = FALSE]
    }
    else if (klags == 0 & kleads == 0) {
        ll <- d_streg
    }

    if (model == 0) {
        xreg <- cbind(
            streg[(klags + 2):(t - kleads), , drop = FALSE],
            ll
        )
    }
    else if (model >= 1 & model <= 4) {
        deter <- determi(model, t, tb)
        xreg <- cbind(
            deter[(klags + 2):(t - kleads), , drop = FALSE],
            streg[(klags + 2):(t - kleads), , drop = FALSE],
            ll
        )
    }
    else if (model == 5) {
        deter <- determi(1, t, tb)
        xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
        xreg <- cbind(
            deter[(klags + 2):(t - kleads), , drop = FALSE],
            streg[(klags + 2):(t - kleads), , drop = FALSE],
            xdu[(klags + 2):(t - kleads), , drop = FALSE],
            ll
        )
    }
    else if (model == 6) {
        deter <- determi(4, t, tb)
        xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
        xreg <- cbind(
            deter[(klags + 2):(t - kleads), , drop = FALSE],
            streg[(klags + 2):(t - kleads), , drop = FALSE],
            xdu[(klags + 2):(t - kleads), , drop = FALSE],
            ll
        )
    }

    beta <- solve(t(xreg) %*% xreg) %*% t(xreg) %*% y[(klags + 2):(t - kleads), 1, drop = FALSE]
    e <- y[(klags + 2):(t - kleads), 1, drop = FALSE] - xreg %*% beta
    s2 <- (t(e) %*% e) / (nrow(xreg) - ncol(xreg))
    t_beta <- sweep(beta, 1, sqrt(diag(s2 * solve(t(xreg) %*% xreg))))

    bic <- ln(s2) + ncol(xreg) * ln(nrow(xreg)) / nrow(xreg)

    return(
        list(
            beta = beta,
            e = e,
            bic = bic,
            t_beta = t_beta
        )
    )
}
