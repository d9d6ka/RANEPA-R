#' @title
#' Procedure to minimize the SSR for 2 break points
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param y (Tx1)-vector of time series
#' @param model A scalar equal to
#' * 1: for the AA (without trend) model,
#' * 2: for the AA (with trend) model,
#' * 3: for the BB model,
#' * 4: for the CC model,
#' * 5: for the AC-CA model.
#'
#' @return A list of
#' * resid: (Tx1) vector of estimated OLS residuals,
#' * tb1: The first break point,
#' * tb2: The second break point.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “The KPSS Test with Two Structural Breaks.”
#' Spanish Economic Review 9, no. 2 (May 16, 2007): 105–27.
#' https://doi.org/10.1007/s10108-006-9017-8.
segments.OLS.double <- function(y, model) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)

    res.r <- 0
    cur.ssr <- Inf
    res.tb1 <- 0
    res.tb2 <- 0

    if (1 <= model & model <= 4) {
        for (i in 2:(N - 4)) {
            for (j in (i + 2):(N - 2)) {
                z <- determinants.KPSS.2.breaks(model, N, c(i, j))
                resids <- OLS(y, z)$residuals
                ssr <- drop(t(resids) %*% resids)
                if (ssr < cur.ssr) {
                    res.r <- resids
                    cur.ssr <- ssr
                    res.tb1 <- i
                    res.tb2 <- j
                }
            }
        }
    } else if (5 <= model & model <= 7) {
        for (i in 2:(N - 4)) {
            for (j in (i + 2):(N - 2)) {
                z <- determinants.KPSS.2.breaks(model, N, c(i, j))
                resids <- OLS(y, z)$residuals
                ssr <- drop(t(resids) %*% resids)
                if (ssr < cur.ssr) {
                    res.r <- resids
                    cur.ssr <- ssr
                    res.tb1 <- i
                    res.tb2 <- j
                }
            }
        }
        for (j in 2:(N - 4)) {
            for (i in (j + 2):(N - 2)) {
                z <- determinants.KPSS.2.breaks(model, N, c(i, j))
                resids <- OLS(y, z)$residuals
                ssr <- drop(t(resids) %*% resids)
                if (ssr < cur.ssr) {
                    res.r <- resids
                    cur.ssr <- ssr
                    res.tb1 <- j
                    res.tb2 <- i
                }
            }
        }
    }

    return(
        list(
            residuals = res.r,
            tb1 = res.tb1,
            tb2 = res.tb2
        )
    )
}
