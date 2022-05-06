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
#' See Carrion-i-Silvestre and Sansó (2006) for further details.
#'
#' @param y (Tx1)-vector of time series.
#' @param model A scalar or vector of
#' \describe{
#' \item{1}{for the break in const.}
#' \item{2}{for the break in trend.}
#' \item{3}{for the break in const and trend.}
#' }
#' @param break_point Array of structural breaks.
#' @param kmax scalar, with the maximum order of the parametric correction. The final order of the parametric correction is selected using the BIC information criterion.
#' @param kernel \describe{
#' \item{bartlett}{for Bartlett kernel.}
#' \item{quadratic}{for Quadratic Spectral kernel.}
#' \item{NULL}{for the Kurozumi's proposal, using Bartlett kernel.}
#' }
#' @param trend Include trend if `TRUE`.
#'
#' @return \describe{
#' \item{beta}{DOLS estimates of the coefficients regressors.}
#' \item{tests}{SC test (coinkpss-test).}
#' \item{resid}{Residuals of the model.}
#' \item{break_point}{Break points.}
#' }
#'
#' @importFrom zeallot %<-%
#'
#' @export
kpss_known_mp <- function(y, model, break_point, kmax, kernel, trend = FALSE) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y) # nolint

    z <- determi_kpss_mp(model, N, break_point, trend)
    c(beta, resid, p, t_b) %<-% olsqr(y, z)

    s_t <- apply(resid, 2, cumsum)

    if (!is.null(kernel))
        test <- N^(-2) * drop(t(s_t) %*% s_t) / alrvr_kernel(resid, kmax, kernel)
    else
        test <- N^(-2) * drop(t(s_t) %*% s_t) / alrvr(resid)

    return(
        list(
            beta = beta,
            test = test,
            resid = resid,
            break_point = break_point
        )
    )
}
