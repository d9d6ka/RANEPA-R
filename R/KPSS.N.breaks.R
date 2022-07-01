#' @title
#' KPSS-test with multiple known structural breaks
#'
#' @description
#' Procedure to compute the KPSS test with multiple known structural breaks
#'
#' @param y An input (LHS) time series of interest.
#' @param x A matrix of (RHS) explanatory stochastic regressors.
#' @param model A scalar or vector of
#' * 1: for the break in const,
#' * 2: for the break in trend,
#' * 3: for the break in const and trend.
#' @param break.point Array of structural breaks.
#' @param const,trend Whether a constant or trend should be included.
#' @param weakly.exog Boolean where we specify
#' whether the stochastic regressors are exogenous or not
#' * `TRUE`: if the regressors are weakly exogenous,
#' * `FALSE`: if the regressors are not weakly exogenous
#' (DOLS is used in this case).
#' @param lags.init,leads.init Scalars defininig the initial number of lags and
#' leads for DOLS.
#' @param max.lag scalar, with the maximum order of the parametric correction.
#' The final order of the parametric correction is selected
#' using the BIC information criterion.
#' @param kernel Kernel for calculating long-run variance
#' * `bartlett`: for Bartlett kernel,
#' * `quadratic`: for Quadratic Spectral kernel,
#' * `NULL` for the Kurozumi's proposal, using Bartlett kernel.
#' @param criterion Information criterion for DOLS lags and leads selection:
#' aic, bic, hq, or lwz.
#'
#' @return A list of
#' * `beta`: DOLS estimates of the coefficients regressors,
#' * `tests`: SC test (coinKPSS-test),
#' * `resid`: Residuals of the model,
#' * `t.beta`: \eqn{t}-statistics for `beta`,
#' * `DOLS.lags`: The estimated number of lags and leads in DOLS,
#' * `break_point`: Break points.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “Testing the Null of Cointegration with Structural Breaks.”
#' Oxford Bulletin of Economics and Statistics 68, no. 5 (October 2006): 623–46.
#' https://doi.org/10.1111/j.1468-0084.2006.00180.x.
#'
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “The KPSS Test with Two Structural Breaks.”
#' Spanish Economic Review 9, no. 2 (May 16, 2007): 105–27.
#' https://doi.org/10.1007/s10108-006-9017-8.
#'
#' @export
KPSS.N.breaks <- function(y, x,
                          model, break.point,
                          const = FALSE, trend = FALSE,
                          weakly.exog = TRUE,
                          lags.init, leads.init,
                          max.lag, kernel,
                          criterion = "bic") {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.null(x)) {
        if (!is.matrix(x)) x <- as.matrix(x)
    }

    N <- nrow(y)

    if (weakly.exog) {
        xt <- cbind(
            x,
            determinants.KPSS.N.breaks(model, N, break.point, const, trend)
        )

        tmp.OLS <- OLS(y, xt)
        beta <- tmp.OLS$beta
        resids <- tmp.OLS$residuals
        t.beta <- tmp.OLS$t.beta
        DOLS.lags <- 0
        DOLS.leads <- 0
        rm(tmp.OLS)
    } else {
        info.crit.min <- Inf
        for (i in lags.init:1) {
            for (j in leads.init:1) {
                tmp.DOLS <- DOLS.N.breaks(
                    y, x, model, break.point, const, trend, i, j
                )
                beta <- tmp.DOLS$beta
                resids <- tmp.DOLS$residuals
                t.beta <- tmp.DOLS$t.beta
                info.crit <- tmp.DOLS$criterions
                if (info.crit[[criterion]] < info.crit.min) {
                    info.crit.min <- info.crit[[criterion]]
                    beta.min <- beta
                    t.beta.min <- t.beta
                    resid.min <- resids
                    DOLS.lags <- i
                    DOLS.leads <- j
                }
                rm(tmp.DOLS)
            }
        }
        resids <- resid.min
        beta <- beta.min
        t.beta <- t.beta.min
    }

    if (!is.null(kernel)) {
        test <- KPSS(resids, lr.var.SPC(resids, max.lag, kernel))
    } else {
        test <- KPSS(resids, lr.var.bartlett.AK(resids))
    }

    return(
        list(
            beta = beta,
            test = test,
            residuals = resids,
            t.beta = t.beta,
            DOLS.lags = DOLS.lags,
            DOLS.leads = DOLS.leads,
            break.point = break.point
        )
    )
}
