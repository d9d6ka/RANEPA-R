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
                           criterion = "aic",
                           modification = FALSE, ...) {
    if (!is.matrix(resid)) resid <- as.matrix(resid)

    N <- nrow(resid)

    if (modification) {
        s2 <- drop(t(resid) %*% resid) / (N - extra)
        tau <- (...elt(1) ^ 2) * drop(t(...elt(2)) %*% ...elt(2)) / s2
    } else {
        tau <- 0
    }

    if (criterion == "aic") {
        result <- log(drop(t(resid) %*% resid) / (N - extra)) +
            2 * (tau + extra) / (N - extra)
    } else if (criterion == "bic") {
        result <- log(drop(t(resid) %*% resid) / (N - extra)) +
            (tau + extra) * log(N - extra) / (N - extra)
    } else if (criterion == "lwz") {
        result <- log(drop(t(resid) %*% resid) / (N - extra)) +
            0.299 * (tau + extra) * (log(N - extra))^2.1
    } else if (criterion == "hq") {
        result <- log(drop(t(resid) %*% resid) / (N - extra)) +
            2 * (tau + extra) * log(log(N - extra)) / (N - extra)
    }

    return(result)
}
