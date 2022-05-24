#' @title
#' SADF test with wild bootstrap.
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#'
#' @export
GSADF.bootstrap.test <- function(y,
                                 r0 = 0.01 + 1.8 / sqrt(length(y)),
                                 const = TRUE,
                                 alpha = 0.05,
                                 iter = 4 * 200,
                                 seed = round(10^4 * sd(y))) {
    N <- length(y)

    ## Find SADF.value.
    model <- GSADF.test(y, r0, const)
    t.values <- model$t.values
    GSADF.value <- model$GSADF.value

    ## Do parallel.
    cores <- detectCores()

    progress.bar <- txtProgressBar(max = iter, style = 3)
    progress <- function(n) setTxtProgressBar(progress.bar, n)

    cluster <- makeCluster(max(cores - 1, 1))
    clusterExport(cluster, c("GSADF.test"))
    registerDoSNOW(cluster)

    GSADF.bootstsrap.values <- foreach(
        step = 1:iter,
        .combine = c,
        .options.snow = list(progress = progress)
    ) %dopar% {
        y.star <- cumsum(c(0, rnorm(N - 1) * diff(y)))
        model <- GSADF.test(y.star, r0, const)
        model$GSADF.value
    }
    stopCluster(cluster)

    ## Find critical value.
    cr.value <- as.numeric(quantile(GSADF.bootstsrap.values, 1 - alpha))

    p.value <- round(sum(GSADF.bootstsrap.values > GSADF.value) / iter, 4)

    result <- list(
        y = y,
        r0 = r0,
        const = const,
        alpha = alpha,
        iter = iter,
        seed = seed,
        t.values = t.values,
        GSADF.value = GSADF.value,
        GSADF.bootstsrap.values = GSADF.bootstsrap.values,
        cr.value = cr.value,
        p.value = p.value
    )

    class(result) <- "sadf"

    return(result)
}
