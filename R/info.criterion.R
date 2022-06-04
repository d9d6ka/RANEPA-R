#' @title
#' Information criterions
#'
#' @details
#' Calculating the value of the following informational criterions:
#' \itemize{
#' \item Akaike,
#' \item Schwarz (Bayesian),
#' \item Hannan-Quinn,
#' \item Liu et al.}
#'
#' @param resid Input residuals.
#' @param extra Number of coefficients etc.
#' @param criterion Needed information criterion: aic, bic, hq or lwz.
#' @param modification Whether the unit-root test modificaton is needed.
#' @param ... Extra parameters needed for modification. Needed are
#' \itemize{
#' \item The coefficient \eqn{\alpha} of \eqn{y_{t-1}},
#' \item The vector of \eqn{y_{t-1}},
#' }
#' exactly in this order.
#'
#' @return
#' The value of the selected informational criterion.
#'
#' @export
info.criterion <- function(resid, extra,
                           modification = FALSE, ...) {
    if (!is.matrix(resid)) resid <- as.matrix(resid)

    N <- nrow(resid)

    if (modification) {
        s2 <- drop(t(resid) %*% resid) / (N - 1)
        tau <- (...elt(1) ^ 2) * drop(t(...elt(2)) %*% ...elt(2)) / s2
    } else {
        tau <- 0
    }

    log.RSS <- log(drop(t(resid) %*% resid) / N)

    return(
        list(
            aic = log.RSS + 2 * (tau + extra) / N,
            bic = log.RSS + (tau + extra) * log(N) / N,
            hq = log.RSS + 2 * (tau + extra) * log(log(N)) / N,
            lwz = log.RSS + 0.299 * (tau + extra) * (log(N))^2.1
        )
    )
}
