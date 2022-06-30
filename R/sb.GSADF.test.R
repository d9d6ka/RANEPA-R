#' @title
#' Sign-based SADF test (HLZ, 2019).
#'
#' @param y A series of interest.
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#'
#' @export
sb.GSADF.test <- function(y,
                          r0 = 0.01 + 1.8 / sqrt(length(y)),
                          const = TRUE,
                          alpha = 0.05,
                          iter = 4 * 200,
                          urs = TRUE,
                          seed = round(10^4 * sd(y))) {
    N <- length(y)

    ## Find supSBADF_value.
    supSBADF.model <- supSBADF.statistic(y, r0)

    ## Do parallel.
    cores <- detectCores()

    progress.bar <- txtProgressBar(max = iter, style = 3)
    progress <- function(n) setTxtProgressBar(progress.bar, n)

    cluster <- makeCluster(max(cores - 1, 1))
    clusterExport(cluster, c(
        "ADF.test",
        "GSADF.test",
        "supSBADF.statistic",
        ".cval_GSADF_without_const",
        ".cval_GSADF_with_const"
    ))
    registerDoSNOW(cluster)

    GSADF.supSBADF.bootstrap.values <- foreach(
        step = 1:iter,
        .combine = rbind,
        .options.snow = list(progress = progress)
    ) %dopar% {
        y.star <- cumsum(c(0, rnorm(N - 1) * diff(y)))
        tmp.GSADF.value <- NA
        supSBADF.value <- NA
        if (urs == TRUE) {
            gsadf.model <- GSADF.test(y.star, r0, const)
            tmp.GSADF.value <- gsadf.model$sadf.value
        }
        supSBADF.model <- supSBADF.statistic(y.star, r0)
        tmp.supSBADF.value <- supSBADF.model$supSBADF.value
        c(tmp.GSADF.value, tmp.supSBADF.value)
    }

    stopCluster(cluster)

    ## Get sadf_supSBADF_bootstrap_values
    supSBADF.bootstrap.values <- GSADF.supSBADF.bootstrap.values[, 2]

    ## Find critical value.
    supSBADF.cr.value <- as.numeric(quantile(
        supSBADF.bootstrap.values,
        1 - alpha
    ))

    ## A union of rejections strategy.
    if (urs == TRUE) {
        ## Find sadf_value.
        gsadf.model <- GSADF.test(y, r0, const)
        t.values <- gsadf.model$t.values
        GSADF.value <- gsadf.model$GSADF.value

        ## Get sadf_supSBADF_bootstrap_values
        GSADF.bootstrap.values <- SGADF.supSBADF.bootstrap.values[, 1]

        ## Find critical value.
        GSADF.cr.value <- as.numeric(quantile(
            GSADF.bootstrap.values,
            1 - alpha
        ))

        ## Calculate U value.
        U.value <- max(
            GSADF.value,
            GSADF.cr.value / supSBADF.cr.value * supSBADF.value
        )

        ## Find U_bootstrap_values.
        U.bootstrap.values <- c()
        for (b in 1:iter) {
            U.bootstrap.values[b] <- max(
                GSADF.bootstrap.values[b],
                GSADF.cr.value / supSBADF.cr.value * supSBADF.bootstrap.values[b]
            )
        }

        ## Find critical value.
        U.cr.value <- as.numeric(quantile(U.bootstrap.values, 1 - alpha))

        p.value <- round(sum(U.bootstrap.values > U.value) / iter, 4)

        is.explosive <- ifelse(U.value > U.cr.value, 1, 0)
    } else {
        p.value <- round(sum(supSBADF.bootstrap.values > supSBADF.value) /
            iter, 4)

        is.explosive <- ifelse(supSBADF.value > supSBADF.cr.value, 1, 0)
    }

    result <- c(
        list(
            y = y,
            r0 = r0,
            const = const,
            alpha = alpha,
            iter = iter,
            urs = urs,
            seed = seed,
            SBADF.values = supSBADF.model$SBADF.values,
            supSBADF.value = supSBADF.model$supSBADF_value,
            supSBADF.bootstrap.values = supSBADF.bootstrap.values,
            supSBADF.cr.value = supSBADF.cr.value,
            p.value = p.value,
            is.explosive = is.explosive
        ),
        if (urs) {
            list(
                t.values = t.values,
                GSADF.value = GSADF.value,
                GSADF.bootstrap.values = SADF.bootstrap.values,
                GSADF.cr.value = SADF.cr.value,
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
