#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom zeallot %<-%
#'
#' @export
bootstrap_kpss_mp <- function(y, ..., iter = 10000) {
    c(., test, resid, ., .) %<-% kpss_known_mp(y, ...)

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
        z <- rnorm(length(resid))
        tmp_y <- z * resid
        c(., tmp_test, ., ., .) %<-% kpss_known_mp(tmp_y, ...)
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
