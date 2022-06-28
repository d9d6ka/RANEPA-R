#' Detrending bootstrap test by Smeekes (2013)
#'
#' @param y The series of interest.
#' @param const Whether the constant is to included.
#' @param trend Whether the trend is to be included.
#' @param c A filtration parameter.
#' @param gamma A detrending type selection parameter. If 0 the OLS detrending
#' is applied, if 1 the GLS detrending is applied, otherwise the autocorrelation
#' coefficient is calculated as \eqn{1 + с^{\gamma} T^{\gamma}}.
#' @param trim The trimming parameter.
#' @param max.lag The maximum lag for inner ADF testing.
#' @param criterion A criterion used to select number of lags.
#' If lag selection is not needed keep this NULL.
#' @param modified.criterion Whether the unit-root test modificaton is needed.
#' @param iter The number of bootstrap steps.
#'
#' @references Smeekes, Stephan.
#' “Detrending Bootstrap Unit Root Tests.”
#' Econometric Reviews 32, no. 8 (July 2013): 869–91.
#' https://doi.org/10.1080/07474938.2012.690693.
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom zeallot %<-%
#'
#' @export
ADF.test.S <- function(y,
                       const = TRUE, trend = FALSE,
                       c = 0, gamma = 0, trim = 0.15,
                       max.lag = 0,
                       criterion = NULL, modified.criterion = FALSE,
                       iter = 999) {
    inner.detrend <- function(y, x, c, g, l) {
        if (is.null(x)) {
            return(y)
        }

        N <- nrow(y)
        beg <- trunc(l * N)
        ct <- (c / N)^g

        yt <- y - (1 - ct) * lagn(y, 1, na = 0)
        xt <- x - (1 - ct) * lagn(x, 1, na = 0)

        c(., yd, ., .) %<-% OLS(
            yt[1:beg, , drop = FALSE],
            xt[1:beg, , drop = FALSE]
        )

        yd <- rbind(
            yd,
            as.matrix(rep(0, N - beg))
        )

        for (lstar in (beg + 1):N) {
            c(., ystar, ., .) %<-% OLS(
                yt[1:lstar, , drop = FALSE],
                xt[1:lstar, , drop = FALSE]
            )
            yd[lstar, ] <- ystar[lstar, ]
        }

        return(yd)
    }

    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)

    x <- NULL
    if (const) {
        x <- cbind(x, rep(1, N))
    }
    if (trend) {
        x <- cbind(x, 1:N)
    }

    yd <- inner.detrend(y, x, c, gamma, trim)

    res.ADF %<-%
        ADF.test(
            yd, const, trend, max.lag,
            criterion, modified.criterion
        )

    res.lag <- res.ADF$lag

    res.ADF %<-%
        ADF.test(
            yd, const, trend, res.lag,
            criterion = NULL
        )

    res.stat <- res.ADF$t.alpha
    res.beta <- res.ADF$beta[-1]

    e <- res.ADF$residuals

    progress.bar <- txtProgressBar(max = iter, style = 3)
    progress <- function(n) setTxtProgressBar(progress.bar, n)

    cores <- detectCores()
    cluster <- makeCluster(max(cores - 1, 1), type = "SOCK")
    registerDoSNOW(cluster)

    tmp.stats <- foreach(
        i = 1:iter,
        .combine = 'c',
        .inorder = FALSE,
        .errorhandling = "remove",
        .packages = c("breaktest"),
        .options.snow = list(progress = progress)
    ) %dopar% {
        u <- rep(0, res.lag + N)
        eps <- sample(e, N, replace = TRUE)

        if (res.lag > 0) {
            for (s in 1:N) {
                u[res.lag + s] <-
                    u[(res.lag + s - 1):s] %*% res.beta + eps[s]
            }
            u <- u[-(1:res.lag)]
        } else {
            for (s in 1:N) {
                u[s] <- eps[s]
            }
        }

        tmp.y <- as.matrix(rep(0, N))

        tmp.y[1] <- u[1]
        for (s in 2:N) {
            tmp.y[s] <- tmp.y[s - 1] + u[s]
        }

        tmp.yd <- inner.detrend(tmp.y, x, c, gamma, trim)

        tmp.res <- ADF.test(
            tmp.yd, const, trend, res.lag,
            criterion = NULL
        )

        tmp.res$t.alpha
    }

    stopCluster(cluster)

    p.value <- sum(tmp.stats < res.stat) / iter

    return(
        list(
            stat = res.stat,
            p.value = p.value
        )
    )
}
