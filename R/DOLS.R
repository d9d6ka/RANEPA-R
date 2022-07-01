#' @title
#' Estimating DOLS regression for multiple known break points
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param y A dependent (LHS) variable.
#' @param x A matrix of explanatory (RHS) variables.
#' @param model See Carrion-i-Silvestre and Sansó (2006)
#' * 1: for model An.
#' * 2: for model A.
#' * 3: for model B.
#' * 4: for model C.
#' * 5: for model D.
#' * 6: for model E.
#' @param break.point A position of the break point.
#' @param k.lags,k.leads A number of lags and leads in DOLS regression.
#'
#' @return A list of:
#' * Estimates of coefficients,
#' * Estimates of residuals,
#' * A value of BIC,
#' * \eqn{t}-statistics for the estimates of coefficients.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “Testing the Null of Cointegration with Structural Breaks.”
#' Oxford Bulletin of Economics and Statistics 68, no. 5 (October 2006): 623–46.
#' https://doi.org/10.1111/j.1468-0084.2006.00180.x.
DOLS <- function(y, x, model, break.point, k.lags, k.leads) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (is.null(x)) {
        stop("ERROR! Explanatory variables needed for DOLS")
    }
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

    if (model == 0) {
        xreg <- cbind(
            x[(k.lags + 2):(N - k.leads), , drop = FALSE],
            lags.leads
        )
    } else if (model >= 1 & model <= 4) {
        deter <- determinants.KPSS.1.break(model, N, break.point)
        xreg <- cbind(
            deter[(k.lags + 2):(N - k.leads), , drop = FALSE],
            x[(k.lags + 2):(N - k.leads), , drop = FALSE],
            lags.leads
        )
    } else if (model == 5) {
        deter <- determinants.KPSS.1.break(1, N, break.point)
        xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
        xreg <- cbind(
            deter[(k.lags + 2):(N - k.leads), , drop = FALSE],
            x[(k.lags + 2):(N - k.leads), , drop = FALSE],
            xdu[(k.lags + 2):(N - k.leads), , drop = FALSE],
            lags.leads
        )
    } else if (model == 6) {
        deter <- determinants.KPSS.1.break(4, N, break.point)
        xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
        xreg <- cbind(
            deter[(k.lags + 2):(N - k.leads), , drop = FALSE],
            x[(k.lags + 2):(N - k.leads), , drop = FALSE],
            xdu[(k.lags + 2):(N - k.leads), , drop = FALSE],
            lags.leads
        )
    }

    res.OLS <- OLS(
        y[(k.lags + 2):(N - k.leads), 1, drop = FALSE],
        xreg
    )

    bic <- log(drop(t(res.OLS$residuals) %*% res.OLS$residuals) / nrow(xreg)) +
        ncol(xreg) * log(nrow(xreg)) / nrow(xreg)

    return(
        list(
            beta   = res.OLS$beta,
            resid  = res.OLS$residuals,
            bic    = bic,
            t.beta = res.OLS$t.beta
        )
    )
}
