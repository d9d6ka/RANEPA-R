#' @import MASS
DOLS <- function(y, x, model, break.point, k.lags, k.leads) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (is.null(x)) stop("ERROR! Explanatory variables needed for DOLS")
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)

    d.x.step <- x[2:N, , drop = FALSE] - x[1:(N - 1), , drop = FALSE]
    d.x.lag <- d.x.step
    d.x.lead <- d.x.step

    for (i in 1:k.lags) {
        d.x.lag <- cbind(
            d.x.lag,
            lagn(d.x.step, i)
        )
    }

    for (i in 1:k.leads) {
        d.x.lead <- cbind(
            d.x.lead,
            lagn(d.x.step, -i)
        )
    }

    if (k.lags != 0 & k.leads != 0) {
        lags <- d.x.lag
        leads <- d.x.lead[, (ncol(x) + 1):(ncol(d.x.lead)), drop = FALSE]
        lags.leads <- cbind(lags, leads)
        lags.leads <- lags.leads[(k.lags + 1):(N - 1 - k.leads), , drop = FALSE]
    }
    else if (k.lags != 0 & k.leads == 0) {
        lags <- d.x.lag
        lags.leads <- lags[(k.lags + 1):(N - 1), , drop = FALSE]
    }
    else if (k.lags == 0 & k.leads != 0) {
        lags <- d.x.lag
        leads <- d.x.lead[, (ncol(x) + 1):(ncol(d.x.lead)), drop = FALSE]
        lags.leads <- cbind(lags, leads)
        lags.leads <- lags.leads[1:(N - 1 - k.leads), , drop = FALSE]
    }
    else if (k.lags == 0 & k.leads == 0) {
        lags.leads <- d.x.lag
    }

    if (model == 0) {
        xreg <- cbind(
            x[(k.lags + 2):(N - k.leads), , drop = FALSE],
            lags.leads
        )
    }
    else if (model >= 1 & model <= 4) {
        deter <- determinants.KPSS.1.break(model, N, break.point)
        xreg <- cbind(
            deter[(k.lags + 2):(N - k.leads), , drop = FALSE],
            x[(k.lags + 2):(N - k.leads), , drop = FALSE],
            lags.leads
        )
    }
    else if (model == 5) {
        deter <- determinants.KPSS.1.break(1, N, break.point)
        xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
        xreg <- cbind(
            deter[(k.lags + 2):(N - k.leads), , drop = FALSE],
            x[(k.lags + 2):(N - k.leads), , drop = FALSE],
            xdu[(k.lags + 2):(N - k.leads), , drop = FALSE],
            lags.leads
        )
    }
    else if (model == 6) {
        deter <- determinants.KPSS.1.break(4, N, break.point)
        xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
        xreg <- cbind(
            deter[(k.lags + 2):(N - k.leads), , drop = FALSE],
            x[(k.lags + 2):(N - k.leads), , drop = FALSE],
            xdu[(k.lags + 2):(N - k.leads), , drop = FALSE],
            lags.leads
        )
    }

    beta <- qr.solve(t(xreg) %*% xreg) %*% t(xreg) %*%
        y[(k.lags + 2):(N - k.leads), 1, drop = FALSE]

    resid <- y[(k.lags + 2):(N - k.leads), 1, drop = FALSE] - xreg %*% beta

    s2 <- drop(t(resid) %*% resid) / (nrow(xreg) - ncol(xreg))

    t.beta <- sweep(beta, 1, sqrt(diag(s2 * qr.solve(t(xreg) %*% xreg))), `/`)

    bic <- log(s2) + ncol(xreg) * log(nrow(xreg)) / nrow(xreg)

    return(
        list(
            beta   = beta,
            resid  = resid,
            bic    = bic,
            t.beta = t.beta
        )
    )
}
