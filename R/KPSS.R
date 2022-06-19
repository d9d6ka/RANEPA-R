#' Auxiliary function returning KPSS statistic value.
#'
#' @param resid The series of residuals.
#' @param variance The value of the long-run variance.
KPSS <- function(resid, variance) {
    if (!is.matrix(resid)) resid <- as.matrix(resid)
    N <- nrow(resid)
    S.t <- apply(resid, 2, cumsum)
    return(drop(t(S.t) %*% S.t) / (N^2 * variance))
}
