#' @title
#' Calculate superior sign-based SADF statistic.
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param y The series of interest.
#' @param trim Trimming parameter to determine the lower and upper bounds.
#' @param generalized Whether to calculate generalized statistic value.
#'
#' @return A list of
#' * `y`,
#' * `trim`,
#' * `C.t`: the cumulative sum of "signs" (1 or -1) of the first difference of
#' `y`,
#' * `SBADF.values`: series of sign-based ADF statistics,
#' * `supSBADF.value`: the maximum of `SBADF.values`.
#'
#' @references
#' Harvey, David I., Stephen J. Leybourne, and Yang Zu.
#' “Sign-Based Unit Root Tests for Explosive Financial Bubbles in the Presence
#' of Deterministically Time-Varying Volatility.”
#' Econometric Theory 36, no. 1 (February 2020): 122–69.
#' https://doi.org/10.1017/S0266466619000057.
supSBADF.statistic <- function(y,
                               trim = 0.01 + 1.8 / sqrt(length(y)),
                               generalized = FALSE) {
    n.obs <- length(y)

    ## Calculate C.t.
    C.t <- cumsum(sign(diff(y)))

    SBADF.values <- c()
    m <- 1

    if (!generalized) {
        for (j in (floor(trim * n.obs)):n.obs) {
            t.beta <- OLS(diff(C.t)[1:j], C.t[1:j])$t.beta
            SBADF.values[m] <- drop(t.beta)
            m <- m + 1
        }
    } else {
        for (i in 1:(n.obs - floor(trim * n.obs) + 1)) {
            for (j in (i + floor(trim * n.obs) - 1):n.obs) {
                t.beta <- OLS(diff(C.t)[i:j], C.t[i:j])$t.beta
                SBADF.values[m] <- drop(t.beta)
                m <- m + 1
            }
        }
    }

    supSBADF.value <- max(SBADF.values)

    return(
        list(
            y = y,
            trim = trim,
            C.t = C.t,
            SBADF.values = SBADF.values,
            supSBADF.value = supSBADF.value
        )
    )
}
