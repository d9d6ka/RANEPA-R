#' @title
#' Generalized supremum ADF test with wild bootstrap.
#'
#' @param y The input time series of interest.
#' @param trim Trimming parameter to determine the lower and upper bounds.
#' @param const Whether the constant needs to be included.
#' @param alpha The significance level of interest.
#' @param iter The number of iterations/
#' @param seed The seed parameter for the random number generator.
#'
#' @references
#' Kurozumi, Eiji, Anton Skrobotov, and Alexey Tsarev.
#' “Time-Transformed Test for the Explosive Bubbles under
#' Non-Stationary Volatility.”
#' arXiv, November 15, 2021. http://arxiv.org/abs/2012.13937.
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#'
#' @export
GSADF.bootstrap.test <- function(y,
                                 trim = 0.01 + 1.8 / sqrt(length(y)),
                                 const = TRUE,
                                 alpha = 0.05,
                                 iter = 4 * 200,
                                 seed = round(10^4 * sd(y))) {
    N <- length(y)

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
        y.star <- cumsum(c(0, rnorm(N - 1) * diff(y)))
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
