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
#' @param resid Input residuals needed for estimating the values of
#' information criterions.
#' @param extra Number of extra parameters needed for estimating the punishment
#' term.
#' @param modification Whether the unit-root test modificaton is needed.
#' See Ng and Perron (2001) for further information.
#' @param alpha The coefficient \eqn{\alpha} of \eqn{y_{t-1}} in ADF model.
#' Needed only for criterion modification purposes.
#' @param y The vector of \eqn{y_{t-1}} in ADF model.
#' Needed only for criterion modification purposes.
#'
#' @return
#' The list of information criterions values.
#'
#' @references
#' Ng, Serena, and Pierre Perron. “Lag Length Selection and the Construction of
#' Unit Root Tests with Good Size and Power.”
#' Econometrica 69, no. 6 (2001): 1519–54.
#' https://doi.org/10.1111/1468-0262.00256.
#'
#' @export
info.criterion <- function(resid, extra,
                           modification = FALSE,
                           alpha = 0, y = NULL) {
    if (!is.matrix(resid)) resid <- as.matrix(resid)

    N <- nrow(resid)

    if (modification) {
        s2 <- drop(t(resid) %*% resid) / N
        tau <- (alpha ^ 2) * drop(t(y) %*% y) / s2
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
