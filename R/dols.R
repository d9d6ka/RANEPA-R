#' @import MASS
dols <- function(y, x, model, break.point, k.lags, k.leads) { # nolint
    if (!is.matrix(y)) y <- as.matrix(y)
    if (is.null(x)) stop("ERROR! Explanatory variables needed for DOLS")
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y) # nolint

    d_x_step <- x[2:N, , drop = FALSE] - x[1:(N - 1), , drop = FALSE] # nolint
    d_x <- d_x_step
    d_x_r <- d_x_step

    for (i in 1:k.lags) {
        d_x <- cbind(
            d_x,
            lagn(d_x_step, i)
        )
    }

    for (i in 1:k.leads) {
        d_x_r <- cbind(
            d_x_r,
            lagn(d_x_step, -i)
        )
    }

    if (k.lags != 0 & k.leads != 0) {
        lags <- d_x
        leads <- d_x_r[, (ncol(x) + 1):(ncol(d_x_r)), drop = FALSE]
        lags_leads <- cbind(lags, leads)
        lags_leads <- lags_leads[(k.lags + 1):(N - 1 - k.leads), , drop = FALSE]
    }
    else if (k.lags != 0 & k.leads == 0) {
        lags <- d_x
        lags_leads <- lags[(k.lags + 1):(N - 1), , drop = FALSE]
    }
    else if (k.lags == 0 & k.leads != 0) {
        lags <- d_x
        leads <- d_x_r[, (ncol(x) + 1):(ncol(d_x_r)), drop = FALSE]
        lags_leads <- cbind(lags, leads)
        lags_leads <- lags_leads[1:(N - 1 - k.leads), , drop = FALSE]
    }
    else if (k.lags == 0 & k.leads == 0) {
        lags_leads <- d_x
    }

    if (model == 0) {
        xreg <- cbind(
            x[(k.lags + 2):(N - k.leads), , drop = FALSE],
            lags_leads
        )
    }
    else if (model >= 1 & model <= 4) {
        deter <- determi_kpss_1p(model, N, break.point)
        xreg <- cbind(
            deter[(k.lags + 2):(N - k.leads), , drop = FALSE],
            x[(k.lags + 2):(N - k.leads), , drop = FALSE],
            lags_leads
        )
    }
    else if (model == 5) {
        deter <- determi_kpss_1p(1, N, break.point)
        xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
        xreg <- cbind(
            deter[(k.lags + 2):(N - k.leads), , drop = FALSE],
            x[(k.lags + 2):(N - k.leads), , drop = FALSE],
            xdu[(k.lags + 2):(N - k.leads), , drop = FALSE],
            lags_leads
        )
    }
    else if (model == 6) {
        deter <- determi_kpss_1p(4, N, break.point)
        xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
        xreg <- cbind(
            deter[(k.lags + 2):(N - k.leads), , drop = FALSE],
            x[(k.lags + 2):(N - k.leads), , drop = FALSE],
            xdu[(k.lags + 2):(N - k.leads), , drop = FALSE],
            lags_leads
        )
    }

    beta <- qr.solve(t(xreg) %*% xreg) %*% t(xreg) %*%
        y[(k.lags + 2):(N - k.leads), 1, drop = FALSE]

    resid <- y[(k.lags + 2):(N - k.leads), 1, drop = FALSE] - xreg %*% beta

    s2 <- drop(t(resid) %*% resid) / (nrow(xreg) - ncol(xreg))

    t_beta <- sweep(beta, 1, sqrt(diag(s2 * qr.solve(t(xreg) %*% xreg))))

    bic <- log(s2) + ncol(xreg) * log(nrow(xreg)) / nrow(xreg)

    return(
        list(
            beta   = beta,
            resid  = resid,
            bic    = bic,
            t_beta = t_beta
        )
    )
}
