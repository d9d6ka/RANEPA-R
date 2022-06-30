#' @title
#' KPSS-test with 2 known structural breaks
#'
#' @description
#' Procedure to compute the KPSS test with two structural breaks
#'
#' @details
#' The break points are known
#'
#' The code provided is the original GAUSS code ported to R.
#'
#' See Carrion-i-Silvestre and Sansó (2007) for further details.
#'
#' @param y An input (LHS) time series of interest.
#' @param model \describe{
#' \item{1}{for the AA (without trend) model.}
#' \item{2}{for the AA (with trend) model.}
#' \item{3}{for the BB model.}
#' \item{4}{for the CC model.}
#' \item{5}{for the AC-CA model.}
#' }
#' @param tb1,tb2 The first and the second break points.
#' @param max.lag scalar, with the maximum order of the parametric correction.
#' The final order of the parametric correction is selected using the
#' BIC information criterion.
#' @param kernel Kernel for calculating long-run variance
#' \describe{
#' \item{bartlett}{for Bartlett kernel.}
#' \item{quadratic}{for Quadratic Spectral kernel.}
#' \item{NULL}{for the Kurozumi's proposal, using Bartlett kernel.}
#' }
#'
#' @return \describe{
#' \item{beta}{DOLS estimates of the coefficients regressors.}
#' \item{tests}{SC test (coinKPSS-test).}
#' \item{resid}{Residuals of the model.}
#' \item{break_point}{Break points.}
#' }
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “The KPSS Test with Two Structural Breaks.”
#' Spanish Economic Review 9, no. 2 (May 16, 2007): 105–27.
#' https://doi.org/10.1007/s10108-006-9017-8.
#'
#' @importFrom zeallot %<-%
#' @export
KPSS.2.breaks <- function(y, model, break.point, max.lag, kernel) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)

    z <- determinants.KPSS.2.breaks(model, N, break.point)
    c(beta, resid, ., t.beta) %<-% OLS(y, z)

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
            break_point = break.point
        )
    )
}


#' @title
#' KPSS-test with 2 unknown structural breaks
#'
#' @description
#' Procedure to compute the KPSS test with two structural breaks
#'
#' @details
#' The break points are known
#'
#' The code provided is the original GAUSS code ported to R.
#'
#' See Carrion-i-Silvestre and Sansó (2007) for further details.
#'
#' @param y (Tx1)-vector of time series.
#' @param model \describe{
#' \item{1}{for the AA (without trend) model.}
#' \item{2}{for the AA (with trend) model.}
#' \item{3}{for the BB model.}
#' \item{4}{for the CC model.}
#' \item{5}{for the AC-CA model.}
#' }
#' @param max.lag scalar, with the maximum order of the parametric correction.
#' The final order of the parametric correction is selected using
#' the BIC information criterion.
#' @param kernel Kernel for calculating long-run variance
#' \describe{
#' \item{bartlett}{for Bartlett kernel.}
#' \item{quadratic}{for Quadratic Spectral kernel.}
#' \item{NULL}{for the Kurozumi's proposal, using Bartlett kernel.}
#' }
#'
#' @return Value of test statistic.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “Testing the Null of Cointegration with Structural Breaks.”
#' Oxford Bulletin of Economics and Statistics 68, no. 5 (October 2006): 623–46.
#' https://doi.org/10.1111/j.1468-0084.2006.00180.x.
#'
#' @importFrom zeallot %<-%
#' @export
KPSS.2.breaks.unknown <- function(y, model, max.lag = 0, kernel = "bartlett") {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)

    c(resid, tb1, tb2) %<-% segments.OLS.double(y, model)

    if (!is.null(kernel)) {
        test <- KPSS(resid, lr.var.SPC(resid, max.lag, kernel))
    } else {
        test <- KPSS(resid, lr.var.bartlett.AK(resid))
    }

    return(
        list(
            test = test,
            tb1 = tb1,
            tb2 = tb2
        )
    )
}
