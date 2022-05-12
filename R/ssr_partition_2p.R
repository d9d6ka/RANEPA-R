#' @title
#' Procedure to minimize the SSR for 2 break points
#'
#' @details
#' See Carrion-i-Silvestre and Sansó (2007) for further details.
#'
#' @param y (Tx1)-vector of time series
#' @param model \describe{
#' \item{1}{for the AA (without trend) model.}
#' \item{2}{for the AA (with trend) model.}
#' \item{3}{for the BB model.}
#' \item{4}{for the CC model.}
#' \item{5}{for the AC-CA model.}
#' }
#'
#' @return List containing
#' \describe{
#' \item{resid}{(Tx1) vector of estimated OLS residuals.}
#' \item{tb1}{The first breaking point.}
#' \item{tb2}{The second breaking point.}
#' }
#'
#' @importFrom zeallot %<-%
ssr_partition_2p <- function(y, model) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)

    r_min <- 0
    ssr_min <- 1000000
    tb1_min <- 0
    tb2_min <- 0

    if (1 <= model & model <= 4) {
        for (i in 2:(N - 4)) {
            for (j in (i + 2):(N - 2)) {
                z <- determi_kpss_2p(model, N, c(i, j))
                c(beta, resid, p) %<-% olsqr(y, z)
                ssr <- drop(t(resid) %*% resid)
                if (ssr < ssr_min) {
                    r_min <- resid
                    ssr_min <- ssr
                    tb1_min <- i
                    tb2_min <- j
                }
            }
        }
    }
    else if (5 <= model & model <= 7) {
        for (i in 2:(N - 4)) {
            for (j in (i + 2):(N - 2)) {
                z <- determi_kpss_2p(model, N, c(i, j))
                c(beta, resid, p) %<-% olsqr(y, z)
                ssr <- drop(t(resid) %*% resid)
                if (ssr < ssr_min) {
                    r_min <- resid
                    ssr_min <- ssr
                    tb1_min <- i
                    tb2_min <- j
                }
            }
        }
        for (j in 2:(N - 4)) {
            for (i in (j + 2):(N - 2)) {
                z <- determi_kpss_2p(model, N, c(i, j))
                c(beta, resid, p) %<-% olsqr(y, z)
                ssr <- drop(t(resid) %*% resid)
                if (ssr < ssr_min) {
                    r_min <- resid
                    ssr_min <- ssr
                    tb1_min <- j
                    tb2_min <- i
                }
            }
        }
    }

    return(
        list(
            resid = r_min,
            tb1 = tb1_min,
            tb2 = tb2_min
        )
    )
}
