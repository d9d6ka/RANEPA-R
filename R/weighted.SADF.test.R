#' @title
#' Weighted supremum ADF test.
#'
#' @param y The input time series of interest.
#' @param trim Trimming parameter to determine the lower and upper bounds.
#' @param const Whether the constant needs to be included.
#' @param alpha The significance level of interest.
#' @param iter The number of iterations.
#' @param urs Use `union of rejections` strategy.
#' @param seed The seed parameter for the random number generator.
#'
#' @references
#' Harvey, David I., Stephen J. Leybourne, and Yang Zu.
#' “Testing Explosive Bubbles with Time-Varying Volatility.”
#' Econometric Reviews 38, no. 10 (November 26, 2019): 1131–51.
#' https://doi.org/10.1080/07474938.2018.1536099.
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#'
#' @export
weighted.SADF.test <- function(y,
                               trim = 0.01 + 1.8 / sqrt(length(y)),
                               const = TRUE,
                               alpha = 0.05,
                               iter = 4 * 200,
                               urs = TRUE,
                               seed = round(10^4 * sd(y))) {
    N <- length(y)

    ## Find supBZ.value.
    supBZ.model <- supBZ.statistic(y, trim)
    sigma.sq <- supBZ.model$sigma.sq
    BZ.values <- supBZ.model$BZ.values
    supBZ.value <- supBZ.model$supBZ.value

    ## Do parallel.
    cores <- detectCores()

    progress.bar <- txtProgressBar(max = iter, style = 3)
    progress <- function(n) setTxtProgressBar(progress.bar, n)

    cluster <- makeCluster(max(cores - 1, 1))
    clusterExport(cluster, c(
        "ADF.test",
        "SADF.test",
        "supBZ.statistic",
        ".cval_SADF_without_const",
        ".cval_SADF_with_const"
    ))
    registerDoSNOW(cluster)

    SADF.supBZ.bootstrap.values <- foreach(
        i = 1:iter,
        .combine = rbind,
        .options.snow = list(progress = progress)
    ) %dopar% {
        y.star <- cumsum(c(0, rnorm(N - 1) * diff(y)))
        tmp.SADF.value <- NA
        if (urs) {
            tmp.sadf.model <- SADF.test(y.star, trim, const)
            tmp.SADF.value <- tmp.sadf.model$SADF.value
        }
        tmp.supBZ.model <- supBZ.statistic(y.star, trim, sigma.sq)
        tmp.supBZ.value <- tmp.supBZ.model$supBZ.value
        c(tmp.SADF.value, tmp.supBZ.value)
    }

    stopCluster(cluster)

    ## Get sadf_supBZ_bootstrap.values.
    supBZ.bootstrap.values <- SADF.supBZ.bootstrap.values[, 2]

    ## Find critical value.
    supBZ.cr.value <- as.numeric(quantile(
        supBZ.bootstrap.values,
        1 - alpha
    ))

    ## A union of rejections strategy.
    if (urs == TRUE) {
        ## Find SADF.value.
        sadf.model <- SADF.test(y, trim, const)
        t.values <- sadf.model$t.values
        SADF.value <- sadf.model$SADF.value

        ## Get sadf_supBZ.bootstrap.values.
        SADF.bootstrap.values <- SADF.supBZ.bootstrap.values[, 1]

        ## Find critical value.
        SADF.cr.value <- as.numeric(quantile(SADF.bootstrap.values, 1 - alpha))

        ## Calculate U value.
        U.value <- max(
            SADF.value,
            SADF.cr.value / supBZ.cr.value * supBZ.value
        )

        ## Find U_bootstrap.values.
        U.bootstrap.values <- c()
        for (b in 1:iter) {
            U.bootstrap.values[b] <- max(
                SADF.bootstrap.values[b],
                SADF.cr.value /
                    supBZ.cr.value *
                    supBZ.bootstrap.values[b]
            )
        }

        ## Find critical value.
        U.cr.value <- as.numeric(quantile(U.bootstrap.values, 1 - alpha))

        p.value <- round(sum(U.bootstrap.values > U.value) / iter, 4)

        is.explosive <- ifelse(U.value > U.cr.value, 1, 0)
    } else {
        p.value <- round(sum(supBZ.bootstrap.values > supBZ.value) / iter, 4)

        is.explosive <- ifelse(supBZ.value > supBZ.cr.value, 1, 0)
    }

    result <- c(
        list(
            y = y,
            trim = trim,
            const = const,
            alpha = alpha,
            iter = iter,
            urs = urs,
            seed = seed,
            sigma.sq = sigma.sq,
            BZ.values = BZ.values,
            supBZ.value = supBZ.value,
            supBZ.bootstrap.values = supBZ.bootstrap.values,
            supBZ.cr.value = supBZ.cr.value,
            p.value = p.value,
            is.explosive = is.explosive
        ),
        if (urs) {
            list(
                t.values = t.values,
                SADF.value = SADF.value,
                SADF.bootstrap.values = SADF.bootstrap.values,
                SADF.cr.value = SADF.cr.value,
                U.value = U.value,
                U.bootstrap.values = U.bootstrap.values,
                U.cr.value = U.cr.value
            )
        } else {
            NULL
        }
    )

    class(result) <- "sadf"

    return(result)
}
