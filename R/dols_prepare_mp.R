dols_prepare_mp <- function(y, x,
                            model, break.point,
                            const = FALSE, trend = FALSE,
                            k.lags, k.leads) { # nolint
    if (!is.matrix(y)) y <- as.matrix(y)
    if (is.null(x)) stop("ERROR! Explanatory variables needed for DOLS")
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y) # nolint

    d_x_step <- x[2:N, , drop = FALSE] - x[1:(N - 1), , drop = FALSE]
    d_x_lag <- d_x_step
    d_x_lead <- d_x_step

    for (i in 1:k.lags) {
        d_x <- cbind(
            d_x_lag,
            lagn(d_x_step, i)
        )
    }

    for (i in 1:k.leads) {
        d_x_lead <- cbind(
            d_x_lead,
            lagn(d_x_step, -i)
        )
    }

    if (k.lags != 0 & k.leads != 0) {
        lags <- d_x_lag
        leads <- d_x_lead[, (ncol(x) + 1):(ncol(d_x_lead)), drop = FALSE]
        lags_leads <- cbind(lags, leads)
        lags_leads <- lags_leads[(k.lags + 1):(N - 1 - k.leads), , drop = FALSE]
    }
    else if (k.lags != 0 & k.leads == 0) {
        lags <- d_x_lag
        lags_leads <- lags[(k.lags + 1):(N - 1), , drop = FALSE]
    }
    else if (k.lags == 0 & k.leads != 0) {
        lags <- d_x_lag
        leads <- d_x_lead[, (ncol(x) + 1):(ncol(d_x_lead)), drop = FALSE]
        lags_leads <- cbind(lags, leads)
        lags_leads <- lags_leads[1:(N - 1 - k.leads), , drop = FALSE]
    }
    else if (k.lags == 0 & k.leads == 0) {
        lags_leads <- d_x_lag
    }
    deter <- determi_kpss_mp(model, N, break.point, const, trend)

    xreg <- cbind(
        deter[(k.lags + 2):(N - k.leads), , drop = FALSE],
        x[(k.lags + 2):(N - k.leads), , drop = FALSE],
        lags_leads
    )

    yreg <- y[(k.lags + 2):(N - k.leads), 1, drop = FALSE]

    return(
        list(
            yreg = yreg,
            xreg = xreg
        )
    )
}
