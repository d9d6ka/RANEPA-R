#' @import MASS
DOLS.N.breaks <- function(y,
                          x,
                          model,
                          break.point,
                          const = FALSE,
                          trend = FALSE,
                          k.lags,
                          k.leads) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (is.null(x)) {
        stop("ERROR! Explanatory variables needed for DOLS")
    }
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)

    c(yreg, xreg) %<-%
        DOLS.vars.N.breaks(
            y, x,
            model, break.point,
            const, trend,
            k.lags, k.leads
        )

    c(beta, resid, ., t.beta) %<-% OLS(yreg, xreg)

    s2 <- drop(t(resid) %*% resid) / nrow(xreg)

    criterions <- info.criterion(resid, ncol(xreg))

    return(
        list(
            beta = beta,
            resid = resid,
            criterions = criterions,
            t.beta = t.beta
        )
    )
}


DOLS.vars.N.breaks <- function(y, x,
                               model, break.point,
                               const = FALSE, trend = FALSE,
                               k.lags, k.leads) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (is.null(x)) stop("ERROR! Explanatory variables needed for DOLS")
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)

    d.x.step <- x[2:N, , drop = FALSE] - x[1:(N - 1), , drop = FALSE]
    d.x.lag <- d.x.step
    d.x.lead <- d.x.step

    for (i in 1:k.lags) {
        d_x <- cbind(
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
    } else if (k.lags != 0 & k.leads == 0) {
        lags <- d.x.lag
        lags.leads <- lags[(k.lags + 1):(N - 1), , drop = FALSE]
    } else if (k.lags == 0 & k.leads != 0) {
        lags <- d.x.lag
        leads <- d.x.lead[, (ncol(x) + 1):(ncol(d.x.lead)), drop = FALSE]
        lags.leads <- cbind(lags, leads)
        lags.leads <- lags.leads[1:(N - 1 - k.leads), , drop = FALSE]
    } else if (k.lags == 0 & k.leads == 0) {
        lags.leads <- d.x.lag
    }
    deter <- determinants.KPSS.N.breaks(model, N, break.point, const, trend)

    xreg <- cbind(
        deter[(k.lags + 2):(N - k.leads), , drop = FALSE],
        x[(k.lags + 2):(N - k.leads), , drop = FALSE],
        lags.leads
    )

    yreg <- y[(k.lags + 2):(N - k.leads), 1, drop = FALSE]

    return(
        list(
            yreg = yreg,
            xreg = xreg
        )
    )
}
