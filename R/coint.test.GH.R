#' Gregory-Hansen test for the absense of cointegration
#'
#' @description
#' Gregory and Hansen (1996) test for the null hypothesis of no cointegration
#' under a possible structural break at the unknown moment of time.
#'
#' The authors proposed ADF- and Z-type tests, slightly modified to allow
#' the presence of a possible regime shift. Three type of shifts are allowed:
#' \itemize{
#' \item a shift in the constant,
#' \item a shift in the constand with the trend included,
#' \item and a shift in the constant and the cointegrating vector.}
#'
#' Critical values are calculated via the adopted MacKinnon procedure of
#' estimating the model for the response surface.
#'
#' @param ... Variables of interest.
#' @param shift Expected break type.
#' @param trim The trimming parameter to calculate break moment bounds.
#' @param max.lag The maximum number of lags for the internal ADF testing.
#' @param criterion The criterion for lag selection.
#' @param add.criticals Whether critical values are to be returned.
#' This argument is needed to suppress the calculation of critical values
#' during the precalculation of tables needed for the p-values estimating.
#'
#' @references
#' MacKinnon, James G.
#' “Critical Values for Cointegration Tests.”
#' Working Paper. Working Paper.
#' Economics Department, Queen’s University, January 2010.
#' https://ideas.repec.org/p/qed/wpaper/1227.html.
#'
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
                          max.lag = 10,
                          criterion = "aic",
                          add.criticals = TRUE) {
    if (...length() < 2) {
        stop("ERROR! Two or more variables are needed")
    }

    if (!shift %in% c("level", "level-trend", "regime")) {
        stop("ERROR! Unknown model specification")
    }

    if (length(unique(sapply(list(...), length))) != 1) {
        stop("ERROR! Series of different length")
    }

    y1 <- as.matrix(...elt(1))
    y2 <- NULL
    for (i in 2:...length()) {
        y2 <- cbind(y2, ...elt(i))
    }

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
                criterion = criterion,
                modified.criterion = TRUE
            )$t.alpha,
            res.ADF
        )
    }

    ## Critical values
    if (add.criticals) {
        m <- min(6, ncol(y2))

        asy.Za <- .cval_coint_gh[[shift]]$Za$b0[m, 3]
        est.Za <- .cval_coint_gh[[shift]]$Za$b0[m, 3] +
            .cval_coint_gh[[shift]]$Za$b1[m, 3] / N

        asy.Zt <- .cval_coint_gh[[shift]]$Zt$b0[m, 3]
        est.Zt <- .cval_coint_gh[[shift]]$Zt$b0[m, 3] +
            .cval_coint_gh[[shift]]$Zt$b1[m, 3] / N

        asy.ADF <- .cval_coint_gh[[shift]]$ADF$b0[m, 3]
        est.ADF <- .cval_coint_gh[[shift]]$ADF$b0[m, 3] +
            .cval_coint_gh[[shift]]$ADF$b1[m, 3] / N

        result <- list(
            shift = shift,
            Za = list(
                stat = res.Za,
                cv = est.Za,
                asy.cv = asy.Za
            ),
            Zt = list(
                stat = res.Zt,
                cv = est.Zt,
                asy.cv = asy.Zt
            ),
            ADF = list(
                stat = res.ADF,
                cv = est.ADF,
                asy.cv = asy.ADF
            )
        )

        class(result) <- "cointGH"

        return(result)
    }

    return(
        list(
            shift = shift,
            Za = res.Za,
            Zt = res.Zt,
            ADF = res.ADF
        )
    )
}
