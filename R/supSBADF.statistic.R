#' @title
#' Calculate superior sign-based SADF statistic.
#'
#' @references
#' Harvey, David I., Stephen J. Leybourne, and Yang Zu.
#' “Sign-Based Unit Root Tests for Explosive Financial Bubbles in the Presence
#' of Deterministically Time-Varying Volatility.”
#' Econometric Theory 36, no. 1 (February 2020): 122–69.
#' https://doi.org/10.1017/S0266466619000057.
#'
#' @importFrom zeallot %<-%
#'
#' @export
supSBADF.statistic <- function(y,
                               trim = 0.01 + 1.8 / sqrt(length(y)),
                               generalized = FALSE) {
    N <- length(y)

    ## Calculate C.t.
    C.t <- cumsum(sign(diff(y)))

    SBADF.values <- c()
    m <- 1

    if (!generalized) {
        for (j in (floor(trim * N)):N) {
            c(., ., ., t.beta) %<-% OLS(diff(C.t)[1:j], C.t[1:j])
            SBADF.values[m] <- drop(t.beta)
            m <- m + 1
        }
    } else {
        for (i in 1:(N - floor(trim * N) + 1)) {
            for (j in (i + floor(trim * N) - 1):N) {
                c(., ., ., t.beta) %<-% OLS(diff(C.t)[i:j], C.t[i:j])
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
