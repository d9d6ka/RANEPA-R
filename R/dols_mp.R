#' @import MASS
dols_mp <- function(y, x, model, break_point, const = FALSE, trend = FALSE, klags, kleads) { # nolint
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y) # nolint

    streg <- x
    d_streg_step <- streg[2:N, , drop = FALSE] - streg[1:(N - 1), , drop = FALSE] # nolint
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
        ll <- cbind(lags, leads)[(klags + 1):(N - 1 - kleads), , drop = FALSE]
    }
    else if (klags != 0 & kleads == 0) {
        lags <- d_streg
        ll <- lags[(klags + 1):(N - 1), , drop = FALSE]
    }
    else if (klags == 0 & kleads != 0) {
        lags <- d_streg
        leads <- d_streg_r[, (ncol(streg) + 1):(ncol(d_streg_r)), drop = FALSE]
        ll <- cbind(lags, leads)[1:(N - 1 - kleads), , drop = FALSE]
    }
    else if (klags == 0 & kleads == 0) {
        ll <- d_streg
    }
    deter <- determi_kpss_mp(model, N, break_point, const, trend)

    xreg <- cbind(
        deter[(klags + 2):(N - kleads), , drop = FALSE],
        streg[(klags + 2):(N - kleads), , drop = FALSE],
        ll
    )
    print(head(xreg))
    beta <- qr.solve(t(xreg) %*% xreg) %*% t(xreg) %*%
        y[(klags + 2):(N - kleads), 1, drop = FALSE]
    resid <- y[(klags + 2):(N - kleads), 1, drop = FALSE] - xreg %*% beta
    s2 <- drop(t(resid) %*% resid) / (nrow(xreg) - ncol(xreg))
    t_beta <- sweep(beta, 1, sqrt(diag(s2 * qr.solve(t(xreg) %*% xreg))))

    bic <- log(s2) + ncol(xreg) * log(nrow(xreg)) / nrow(xreg)

    print(paste(klags, " : ", bic, sep=""))

    return(
        list(
            beta = beta,
            resid = resid,
            bic = bic,
            t_beta = t_beta
        )
    )
}
