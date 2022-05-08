#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom zeallot %<-%
#'
#' @export
bootstrap_kpss_mp <- function(y, x,
                              model, break.point,
                              const = FALSE, trend = FALSE,
                              weakly.exog = TRUE,
                              ll.init, corr.max, kernel,
                              iter = 10000) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.null(x))
        if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)

    c(., test, u, ., dols_lags, .) %<-%
        kpss_known_mp(
            y, x,
            model, break.point,
            const, trend,
            weakly.exog,
            ll.init, corr.max, kernel
        )

    if (weakly.exog) {
        xreg <- cbind(
            x,
            determi_kpss_mp(model, N, break.point, const, trend)
        )
    }
    else {
        c(., xreg) %<-%
            dols_prepare_mp(
                y, x,
                model, break.point,
                const, trend,
                dols_lags, dols_lags
            )
    }

    cores <- detectCores() # nolint

    progress_bar <- txtProgressBar(max = iter, style = 3)
    progress <- function(n) setTxtProgressBar(progress_bar, n)

    cluster <- makeCluster(max(cores - 1, 1))
    registerDoSNOW(cluster)

    result <- foreach(
        i = 1:iter,
        .combine = rbind,
        .options.snow = list(progress = progress)
    ) %dopar% {
        z <- rnorm(length(u))
        tmp_y <- z * u

        c(beta, resid, ., t_beta) %<-% olsqr(tmp_y, xreg)

        s_t <- apply(resid, 2, cumsum)
        if (!is.null(kernel))
            tmp_test <- N^(-2) * drop(t(s_t) %*% s_t) /
                alrvr_kernel(resid, corr.max, kernel)
        else
            tmp_test <- N^(-2) * drop(t(s_t) %*% s_t) / alrvr(resid)
        tmp_test
    }

    stopCluster(cluster)

    p_value <- (1 / iter) * sum(I(test <= result))

    return(
        list(
            test = test,
            p_value = p_value,
            bootstrapped = result
        )
    )
}
