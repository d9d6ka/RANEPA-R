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
segs.SSR.2.breaks <- function(y, model) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)

    r.min <- 0
    SSR.min <- 1000000
    tb1.min <- 0
    tb2.min <- 0

    if (1 <= model & model <= 4) {
        for (i in 2:(N - 4)) {
            for (j in (i + 2):(N - 2)) {
                z <- determinants.KPSS.2.breaks(model, N, c(i, j))
                c(., resid, p, .) %<-% OLS(y, z)
                SSR <- drop(t(resid) %*% resid)
                if (SSR < SSR.min) {
                    r.min <- resid
                    SSR.min <- SSR
                    tb1.min <- i
                    tb2.min <- j
                }
            }
        }
    } else if (5 <= model & model <= 7) {
        for (i in 2:(N - 4)) {
            for (j in (i + 2):(N - 2)) {
                z <- determinants.KPSS.2.breaks(model, N, c(i, j))
                c(., resid, p, .) %<-% OLS(y, z)
                SSR <- drop(t(resid) %*% resid)
                if (SSR < SSR.min) {
                    r.min <- resid
                    SSR.min <- SSR
                    tb1.min <- i
                    tb2.min <- j
                }
            }
        }
        for (j in 2:(N - 4)) {
            for (i in (j + 2):(N - 2)) {
                z <- determinants.KPSS.2.breaks(model, N, c(i, j))
                c(., resid, p, .) %<-% OLS(y, z)
                SSR <- drop(t(resid) %*% resid)
                if (SSR < SSR.min) {
                    r.min <- resid
                    SSR.min <- SSR
                    tb1.min <- j
                    tb2.min <- i
                }
            }
        }
    }

    return(
        list(
            resid = r.min,
            tb1 = tb1.min,
            tb2 = tb2.min
        )
    )
}
