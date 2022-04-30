#' @import MASS
dols <- function(y, x, model, klags, kleads, tb) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    t <- nrow(y)

    streg <- x
    d_streg_step <- streg[2:t, , drop = FALSE] - streg[1:(t - 1), , drop = FALSE]
    d_streg <- d_streg_step
    d_streg_r <- d_streg_step

    for (i in 1:klags) {
        d_streg <- cbind(
            d_streg,
            lagn(d_streg_step, i)
        )
    }

    for (i in 1:kleads) {
        d_streg_r <- cbind(
            d_streg_r,
            lagn(d_streg_step, -i)
        )
    }

    if (klags != 0 & kleads != 0) {
        lags <- d_streg
        leads <- d_streg_r[, (ncol(streg) + 1):(ncol(d_streg_r)), drop = FALSE]
        ll <- cbind(lags, leads)[(klags + 1):(t - 1 - kleads), , drop = FALSE]
    }
    else if (klags != 0 & kleads == 0) {
        lags <- d_streg
        ll <- lags[(klags + 1):(t - 1), , drop = FALSE]
    }
    else if (klags == 0 & kleads != 0) {
        lags <- d_streg
        leads <- d_streg_r[, (ncol(streg) + 1):(ncol(d_streg_r)), drop = FALSE]
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

    beta <- qr.solve(t(xreg) %*% xreg) %*% t(xreg) %*%
        y[(klags + 2):(t - kleads), 1, drop = FALSE]
    e <- y[(klags + 2):(t - kleads), 1, drop = FALSE] - xreg %*% beta
    s2 <- as.numeric(t(e) %*% e) / (nrow(xreg) - ncol(xreg))
    t_beta <- sweep(beta, 1, sqrt(diag(s2 * qr.solve(t(xreg) %*% xreg))))

    bic <- log(s2) + ncol(xreg) * log(nrow(xreg)) / nrow(xreg)

    return(
        list(
            beta = beta,
            e = e,
            bic = bic,
            t_beta = t_beta
        )
    )
}
