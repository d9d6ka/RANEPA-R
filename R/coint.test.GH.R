#' Gregory-Hansen test for the absense of cointegration
#'
#' @description
#' Gregory-Hansen test for the absense of cointegration under the assumption
#' of a structural break at the unknown moment of time.
#'
#' @param ... Variables of interest.
#' @param shift Expected break type.
#' @param trim The trimming parameter to calculate break moment bounds.
#' @param max.lag The maximum number of lags for the internal ADF testing.
#'
#' @references
#' Gregory, Allan W., and Bruce E. Hansen.
#' “Residual-Based Tests for Cointegration in Models with Regime Shifts.”
#' Journal of Econometrics 70, no. 1 (January 1, 1996): 99–126.
#' https://doi.org/10.1016/0304-4076(69)41685-7.
#'
#' @importFrom zeallot %<-%
#'
#' @export
coint.test.GH <- function(...,
                          shift = "level",
                          trim = 0.15,
                          max.lag = 10) {
    if (...length() < 2) {
        stop("ERROR! Two or more variables are needed")
    }

    if (!shift %in% c("level", "level-trend", "regime")) {
        stop("ERROR! Unknown model specification")
    }

    y1 <- as.matrix(...elt(1))
    y2 <- NULL
    for (i in 2:...length())
        y2 <- cbind(y2, ...elt(i))

    N <- nrow(y1)
    x.const <- rep(1, N)
    x.trend <- 1:N

    first.break <- trunc(trim * N)
    last.break <- trunc((1 - trim) * N)

    res.Za <- Inf
    res.Zt <- Inf
    res.ADF <- Inf

    for (tb in first.break:last.break) {
        phi <- as.numeric(x.trend > tb)

        if (shift == "level") {
            x <- cbind(
                x.const,
                phi,
                y2
            )
        } else if (shift == "level-trend") {
            x <- cbind(
                x.const,
                phi,
                x.trend,
                y2
            )
        } else if (shift == "regime") {
            x <- cbind(
                x.const,
                phi,
                y2,
                phi * y2
            )
        }

        c(., e, ., .) %<-% OLS(y1, x)

        rho <- sum(e[1:(N - 1), ] * e[2:N, ]) /
            sum(e[1:(N - 1), ]^2)

        nu <- e - rho * lagn(e, 1, na = 0)

        lrv <- lr.var(nu)
        lambda <- (lrv - drop(t(nu) %*% nu) / N) / 2

        rho.star <- sum(e[1:(N - 1), ] * e[2:N, ] - lambda) /
            sum(e[1:(N - 1), ]^2)

        res.Za <- min(
            N * (rho.star - 1),
            res.Za
        )
        res.Zt <- min(
            (rho.star - 1) * sqrt(sum(e[1:(N - 1), ]^2) / lrv),
            res.Zt
        )

        res.ADF <- min(
            ADF.test(
                e,
                const = FALSE, trend = FALSE,
                max.lag = max.lag,
                criterion = "aic",
                modified.criterion = TRUE
            )$t.alpha,
            res.ADF
        )
    }

    return(
        list(
            Za = res.Za,
            Zt = res.Zt,
            ADF = res.ADF
        )
    )
}
