#' @title
#' ADF.test - simple implementation of ADF test
#'
#' @details
#' A function for ADF test with the ability to select the number of lags.
#' Lags are selected by informational criterions which can be modified as in
#' Ng and Perron (2001) and Cavaliere et al. (2015).
#'
#' Due to the Frisch-Waugh-Lovell theorem we first detrend `y` and then apply
#' the test to the detrended series.
#'
#' @param y The input time series of interest.
#' @param const Include const to the model if TRUE.
#' @param trend Include trend to the model if TRUE.
#' @param max.lag Maximum lag number
#' @param criterion A criterion used to select number of lags.
#' If lag selection is not needed keep this NULL.
#' @param modified.criterion Whether the unit-root test modificaton is needed.
#' @param rescale.criterion Whether the rescaling informational criterion
#' is needed. Designed to cope with heteroscedasticity in residuals.
#'
#' @return List containing:
#' \itemize{
#' \item y
#' \item const
#' \item trend
#' \item residuals
#' \item coefficient estimates
#' \item t-statistic value
#' \item critical value
#' \item Number of lags
#' \item indicator of stationarity}
#'
#' @references
#' Cavaliere, Giuseppe, Peter C. B. Phillips, Stephan Smeekes,
#' and A. M. Robert Taylor. “Lag Length Selection for Unit Root Tests
#' in the Presence of Nonstationary Volatility.”
#' Econometric Reviews 34, no. 4 (April 21, 2015): 512–36.
#' https://doi.org/10.1080/07474938.2013.808065.
#'
#' Ng, Serena, and Pierre Perron. “Lag Length Selection and the Construction of
#' Unit Root Tests with Good Size and Power.”
#' Econometrica 69, no. 6 (2001): 1519–54.
#' https://doi.org/10.1111/1468-0262.00256.
#'
#' @importFrom zeallot %<-%
#'
#' @export
ADF.test <- function(y,
                     const = TRUE, trend = FALSE,
                     max.lag = 0,
                     criterion = NULL, modified.criterion = FALSE,
                     rescale.criterion = FALSE) {
    if (!is.matrix(y)) y <- as.matrix(y)

    if (!is.null(criterion)) {
        if (!criterion %in% c("bic", "aic", "lwz", "hq")) {
            warning("WARNING! Unknown criterion, none is used")
            criterion <- NULL
        }
    }

    N <- nrow(y)

    deter <- NULL
    if (const) {
        deter <- cbind(deter, rep(1, N))
    }
    if (trend) {
        deter <- cbind(deter, 1:N)
    }

    ## Detrending
    if (!is.null(deter)) {
        c(., yd, ., .) %<-% OLS(y, deter)
    } else {
        yd <- y
    }

    d.y <- as.matrix(c(yd[1], diff(yd)))

    x <- lagn(yd, 1, na = 0)
    if (max.lag > 0) {
        for (l in 1:max.lag) {
            x <- cbind(x, lagn(d.y, l, na = 0))
        }
    }

    if (is.null(criterion)) {
        res.lag <- max.lag
    } else {
        if (rescale.criterion) {
            c(d.yr, xr) %<-% CPST.rescale(d.y, x, deter, 0, max.lag)
        } else {
            d.yr <- d.y
            xr <- x
        }

        c(b, e, ., .) %<-%
            OLS(
                d.yr[(2 + max.lag):N, , drop = FALSE],
                xr[(2 + max.lag):N, 1, drop = FALSE]
            )

        res.ic <- info.criterion(
            e, 0,
            modification = modified.criterion,
            alpha = b[1], y = xr[(2 + max.lag):N, 1, drop = FALSE]
        )[[criterion]]
        res.lag <- 0

        for (l in 1:max.lag) {
            if (max.lag == 0) break

            if (rescale.criterion)
                c(d.yr, xr) %<-% CPST.rescale(d.y, x, deter, l, max.lag)
            else {
                d.yr <- d.y
                xr <- x
            }

            c(b, e, ., .) %<-% OLS(
                d.yr[(2 + max.lag):N, , drop = FALSE],
                xr[(2 + max.lag):N, 1:(1 + l), drop = FALSE]
            )

            tmp.ic <- info.criterion(
                e, l,
                modification = modified.criterion,
                alpha = b[1], y = xr[(2 + max.lag):N, 1, drop = FALSE]
            )[[criterion]]

            if (tmp.ic < res.ic) {
                res.ic <- tmp.ic
                res.lag <- l
            }
        }
    }

    c(res.beta, res.resid, ., res.t.beta) %<-%
        OLS(
            d.y[(2 + max.lag):N, , drop = FALSE],
            x[(2 + max.lag):N, 1:(1 + res.lag), drop = FALSE]
        )

    Z.stat <- (N - max.lag - 1) * drop(res.beta[1] - 1)

    return(
        list(
            y = drop(y),
            yd = drop(yd),
            const = const,
            trend = trend,
            beta = res.beta,
            t.beta = drop(res.t.beta),
            alpha = drop(res.beta[1]),
            t.alpha = drop(res.t.beta[1]),
            Z.stat = Z.stat,
            residuals = res.resid,
            lag = res.lag
        )
    )
}


#' Generating rescaled series as in Cavaliere et al. (2015).
#'
#' @param d.y The series of first differences.
#' @param x The matrix of ADF RHS variables.
#' @param deter The matrix of deterministic variables for detrending.
#' @param k The lag of the corresponding ADF model.
#' @param max.lag The maximum possible lag.
#'
#' @references
#' Cavaliere, Giuseppe, Peter C. B. Phillips, Stephan Smeekes,
#' and A. M. Robert Taylor. “Lag Length Selection for Unit Root Tests
#' in the Presence of Nonstationary Volatility.”
#' Econometric Reviews 34, no. 4 (April 21, 2015): 512–36.
#' https://doi.org/10.1080/07474938.2013.808065.
#'
#' @importFrom zeallot %<-%
CPST.rescale <- function(d.y, x, deter, k, max.lag) {
    c(., e, ., .) %<-% OLS(d.y, x[, 1:(1 + k), drop = FALSE])

    NW.se <- NW.volatility(
        e,
        NW.loocv(e^2, rep(1, nrow(e)))$h
    )$se

    yr <- cumsum(d.y / NW.se)

    if (!is.null(deter))
        c(., yr, ., .) %<-% OLS(yr, deter)

    d.yr <- as.matrix(c(0, diff(yr)))

    xr <- lagn(yr, 1, na = 0)
    if (max.lag > 0) {
        for (l in 1:max.lag) {
            xr <- cbind(xr, lagn(d.yr, l, na = 0))
        }
    }

    return(list(d.y = d.yr, x = xr))
}
