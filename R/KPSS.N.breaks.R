#' @title
#' KPSS-test with multiple known structural breaks
#'
#' @description
#' Procedure to compute the KPSS test with multiple known structural breaks
#'
#' @param y An input (LHS) time series of interest.
#' @param x A matrix of (RHS) explanatory stochastic regressors.
#' @param model A scalar or vector of
#' \describe{
#' \item{1}{for the break in const.}
#' \item{2}{for the break in trend.}
#' \item{3}{for the break in const and trend.}
#' }
#' @param break.point Array of structural breaks.
#' @param trend Include trend if `TRUE`.
#' @param weakly.exog Boolean where we specify
#' whether the stochastic regressors are exogenous or not
#' \describe{
#' \item{TRUE}{if the regressors are weakly exogenous,}
#' \item{FALSE}{if the regressors are not weakly exogenous
#' (DOLS is used in this case).}
#' }
#' @param lags.init,leads.init Scalars defininig the initial number of lags and
#' leads for DOLS.
#' @param max.lag scalar, with the maximum order of the parametric correction.
#' The final order of the parametric correction is selected
#' using the BIC information criterion.
#' @param kernel Kernel for calculating long-run variance
#' \describe{
#' \item{bartlett}{for Bartlett kernel.}
#' \item{quadratic}{for Quadratic Spectral kernel.}
#' \item{NULL}{for the Kurozumi's proposal, using Bartlett kernel.}
#' }
#' @param criterion Information criterion for DOLS lags and leads selection:
#' aic, bic, hq, or lwz.
#'
#' @return \describe{
#' \item{beta}{DOLS estimates of the coefficients regressors.}
#' \item{tests}{SC test (coinKPSS-test).}
#' \item{resid}{Residuals of the model.}
#' \item{t.beta}{t-statistics for `beta`.}
#' \item{DOLS.lags}{The estimated number of lags and leads in DOLS.}
#' \item{break_point}{Break points.}
#' }
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
#' @importFrom zeallot %<-%
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

        c(beta, resid, ., t.beta) %<-% OLS(y, xt)
        DOLS.lags <- 0
        DOLS.leads <- 0
    } else {
        info.crit.min <- Inf
        for (i in lags.init:1) {
            for (j in leads.init:1) {
                c(beta, resid, info.crit, t.beta) %<-%
                    DOLS.N.breaks(y, x, model, break.point, const, trend, i, j)
                if (info.crit[[criterion]] < info.crit.min) {
                    info.crit.min <- info.crit[[criterion]]
                    beta.min <- beta
                    t.beta.min <- t.beta
                    resid.min <- resid
                    DOLS.lags <- i
                    DOLS.leads <- j
                }
            }
        }
        resid <- resid.min
        beta <- beta.min
        t.beta <- t.beta.min
    }

    if (!is.null(kernel)) {
        test <- KPSS(resid, lr.var.SPC(resid, max.lag, kernel))
    } else {
        test <- KPSS(resid, lr.var.bartlett.AK(resid))
    }

    return(
        list(
            beta = beta,
            test = test,
            resid = resid,
            t.beta = t.beta,
            DOLS.lags = DOLS.lags,
            DOLS.leads = DOLS.leads,
            break.point = break.point
        )
    )
}
