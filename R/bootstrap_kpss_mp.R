#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom zeallot %<-%
#'
#' @export
bootstrap_kpss_mp <- function(y, x,
                          model, break.point, const = FALSE, trend = FALSE,
                          weakly.exog = TRUE,
                          ll.init, corr.max, kernel, iter = 10000) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.null(x))
        if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y) # nolint

    c(., test, resid, ., .) %<-%
        kpss_known_mp(
            y, x, model, break.point,
            const, trend, weakly.exog,
            ll.init, corr.max, kernel
        )

    cores <- detectCores() # nolint

    progress_bar <- txtProgressBar(max = length(iter), style = 3)
    progress <- function(n) setTxtProgressBar(progress_bar, n)

    cluster <- makeCluster(max(cores - 1, 1))
    registerDoSNOW(cluster)

    result <- foreach(
        i = 1:iter,
        .combine = rbind,
        .options.snow = list(progress = progress)
    ) %dopar% {
        z <- rnorm(length(resid))
        tmp_y <- z * resid
        c(., tmp_test, ., ., .) %<-%
            kpss_known_mp(
                tmp_y, x, model, break.point,
                const, trend, weakly.exog,
                ll.init, corr.max, kernel
            )
        tmp_test
    }
    
    stopCluster(cluster)
    
    p_value <- (1 / N) *  sum(I(test <= result))

    return(
        list(
            test = test,
            p_value = p_value,
            bootstrapped = result
        )
    )
}
