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
#' @param y (Tx1)-vector of time series
#' @param model
#' \describe{
#' \item{1}{for the AA (without trend) model}
#' \item{2}{for the AA (with trend) model}
#' \item{3}{for the BB model}
#' \item{4}{for the CC model}
#' \item{5}{for the AC-CA model}
#' }
#' @param tb1 First break point
#' @param tb2 Second break point
#' @param kmax scalar, with the maximum order of the parametric correction. The final order of the parametric correction is selected using the BIC information criterion.
#' @param kernel
#' \describe{
#' \item{bartlett}{for Bartlett kernel}
#' \item{quadratic}{for Quadratic Spectral kernel}
#' \item{NULL}{for the Kurozumi's proposal, using Bartlett kernel}
#' }
#'
#' @return Value of test statistic
#'
#' @importFrom zeallot %<-%
#' @export
kpss_known_2p <- function(y, model, tb1, tb2, kmax, kernel) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)

    z <- determi_kpss_2p(model, N, tb1, tb2)
    c(b, e, p) %<-% olsqr(y, z)

    s_t <- apply(e, 2, cumsum)

    if (!is.null(kernel))
        test <- N^(-2) * drop(t(s_t) %*% s_t) / alrvr_kernel(e, kmax, kernel)
    else
        test <- N^(-2) * drop(t(s_t) %*% s_t) / alrvr(e)

    return(test)
}
