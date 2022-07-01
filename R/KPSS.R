#' @title
#' Auxiliary function returning KPSS statistic value.
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param resids The series of residuals.
#' @param variance The value of the long-run variance.
KPSS <- function(resids, variance) {
    if (!is.matrix(resids)) resids <- as.matrix(resids)
    N <- nrow(resids)
    S.t <- apply(resids, 2, cumsum)
    return(drop(t(S.t) %*% S.t) / (N^2 * variance))
}
