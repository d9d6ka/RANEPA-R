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
#' @param model A scalar equal to
#' * 1: for the AA (without trend) model,
#' * 2: for the AA (with trend) model,
#' * 3: for the BB model,
#' * 4: for the CC model,
#' * 5: for the AC-CA model.
#' @param break.point Positions for the first and second structural breaks
#' (respective to the origin which is 1).
#' @param max.lag scalar, with the maximum order of the parametric correction.
#' The final order of the parametric correction is selected using the
#' BIC information criterion.
#' @param kernel Kernel for calculating long-run variance
#' * `bartlett`: for Bartlett kernel,
#' * `quadratic`: for Quadratic Spectral kernel,
#' * `NULL` for the Kurozumi's proposal, using Bartlett kernel.
#'
#' @return A list of:
#' * `beta`: DOLS estimates of the coefficients regressors,
#' * `tests`: SC test (coinKPSS-test),
#' * `resid`: Residuals of the model,
#' * `t.beta`: \eqn{t}-statistics for `beta`,
#' * `break_point`: Break points.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “The KPSS Test with Two Structural Breaks.”
#' Spanish Economic Review 9, no. 2 (May 16, 2007): 105–27.
#' https://doi.org/10.1007/s10108-006-9017-8.
#'
#' @export
KPSS.2.breaks <- function(y, model, break.point, max.lag, kernel) {
    if (!is.matrix(y)) y <- as.matrix(y)

    n.obs <- nrow(y)

    z <- determinants.KPSS.2.breaks(model, n.obs, break.point)

    res.OLS <- OLS(y, z)

    if (!is.null(kernel)) {
        test <- KPSS(
            res.OLS$residuals,
            lr.var.SPC(res.OLS$residuals, max.lag, kernel)
        )
    } else {
        test <- KPSS(
            res.OLS$residuals,
            lr.var.bartlett.AK(res.OLS$residuals)
        )
    }

    return(
        list(
            beta = res.OLS$beta,
            test = test,
            residuals = res.OLS$residuals,
            t.beta = res.OLS$t.beta,
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
#' @param model A scalar equal to
#' * 1: for the AA (without trend) model,
#' * 2: for the AA (with trend) model,
#' * 3: for the BB model,
#' * 4: for the CC model,
#' * 5: for the AC-CA model.
#' @param max.lag scalar, with the maximum order of the parametric correction.
#' The final order of the parametric correction is selected using
#' the BIC information criterion.
#' @param kernel Kernel for calculating long-run variance
#' * `bartlett`: for Bartlett kernel,
#' * `quadratic`: for Quadratic Spectral kernel,
#' * `NULL` for the Kurozumi's proposal, using Bartlett kernel.
#'
#' @return Value of test statistic.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “Testing the Null of Cointegration with Structural Breaks.”
#' Oxford Bulletin of Economics and Statistics 68, no. 5 (October 2006): 623–46.
#' https://doi.org/10.1111/j.1468-0084.2006.00180.x.
#'
#' @export
KPSS.2.breaks.unknown <- function(y, model, max.lag = 0, kernel = "bartlett") {
    if (!is.matrix(y)) y <- as.matrix(y)

    res.segs <- segments.OLS.double(y, model)

    if (!is.null(kernel)) {
        test <- KPSS(
            res.segs$residuals,
            lr.var.SPC(res.segs$residuals, max.lag, kernel)
        )
    } else {
        test <- KPSS(
            res.segs$residuals,
            lr.var.bartlett.AK(res.segs$residuals)
        )
    }

    return(
        list(
            test = test,
            tb1 = res.segs$tb1,
            tb2 = res.segs$tb2
        )
    )
}
