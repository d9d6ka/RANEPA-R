#' @rdname SADF.bootstrap.test
#' @order 2
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom stats quantile
#' @importFrom stats rnorm
#' @importFrom stats sd
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @export
GSADF.bootstrap.test <- function(y,
                                 trim = 0.01 + 1.8 / sqrt(length(y)),
                                 const = TRUE,
                                 alpha = 0.05,
                                 iter = 4 * 200,
                                 seed = round(10^4 * sd(y))) {
    n.obs <- length(y)

    ## Find SADF.value.
    model <- GSADF.test(y, trim, const)
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
        y.star <- cumsum(c(0, rnorm(n.obs - 1) * diff(y)))
        model <- GSADF.test(y.star, trim, const)
        model$GSADF.value
    }
    stopCluster(cluster)

    ## Find critical value.
    cr.value <- as.numeric(quantile(GSADF.bootstsrap.values, 1 - alpha))

    p.value <- round(sum(GSADF.bootstsrap.values > GSADF.value) / iter, 4)

    result <- list(
        y = y,
        trim = trim,
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
