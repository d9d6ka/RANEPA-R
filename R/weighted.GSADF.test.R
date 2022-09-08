#' @rdname weighted.SADF.test
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
weighted.GSADF.test <- function(y,
                                trim = 0.01 + 1.8 / sqrt(length(y)),
                                const = TRUE,
                                alpha = 0.05,
                                iter = 4 * 200,
                                urs = TRUE,
                                seed = round(10^4 * sd(y))) {
    n.obs <- length(y)

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
        "GSADF.test",
        "supBZ.statistic",
        ".cval_GSADF_without_const",
        ".cval_GSADF_with_const"
    ))
    registerDoSNOW(cluster)

    SADF.supBZ.bootstrap.values <- foreach(
        step = 1:iter,
        .combine = rbind,
        .options.snow = list(progress = progress)
    ) %dopar% {
        y.star <- cumsum(c(0, rnorm(n.obs - 1) * diff(y)))
        tmp.GSADF.value <- NA
        if (urs == TRUE) {
            gsadf.model <- GSADF.test(y.star, trim, const)
            tmp.GSADF.value <- gsadf.model$GSADF.value
        }
        supBZ.model <- supBZ.statistic(y.star, trim, sigma.sq)
        tmp.supBZ.value <- supBZ.model$supBZ.value
        c(tmp.GSADF.value, tmp.supBZ.value)
    }

    stopCluster(cluster)

    ## Get sadf_supBZ.bootstsrap.values.
    supBZ.bootstsrap.values <- SADF.supBZ.bootstrap.values[, 2]

    ## Find critical value.
    supBZ.cr.value <- as.numeric(quantile(
        supBZ.bootstsrap.values,
        1 - alpha
    ))

    ## A union of rejections strategy.
    if (urs == TRUE) {
        ## Find sadf.value.
        gsadf.model <- GSADF.test(y, trim, const)
        t.values <- gsadf.model$t.values
        GSADF.value <- gsadf.model$GSADF.value

        ## Get sadf_supBZ.bootstsrap.values.
        GSADF.bootstsrap.values <- SADF.supBZ.bootstrap.values[, 1]

        ## Find critical value.
        GSADF.cr.value <- as.numeric(quantile(
            GSADF.bootstsrap.values,
            1 - alpha
        ))

        ## Calculate U value.
        U.value <- max(
            GSADF.value,
            GSADF.cr.value / supBZ.cr.value * supBZ.value
        )

        ## Find U.bootstsrap.values.
        U.bootstsrap.values <- c()
        for (b in 1:iter) {
            U.bootstsrap.values[b] <- max(
                GSADF.bootstsrap.values[b],
                GSADF.cr.value / supBZ.cr.value * supBZ.bootstsrap.values[b]
            )
        }

        ## Find critical value.
        U.cr.value <- as.numeric(quantile(U.bootstsrap.values, 1 - alpha))

        p.value <- round(sum(U.bootstsrap.values > U.value) / iter, 4)

        is.explosive <- ifelse(U.value > U.cr.value, 1, 0)
    } else {
        p.value <- round(sum(supBZ.bootstsrap.values > supBZ.value) / iter, 4)

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
            supBZ.bootstsrap.values = supBZ.bootstsrap.values,
            supBZ.cr.value = supBZ.cr.value,
            p.value = p.value,
            is.explosive = is.explosive
        ),
        if (urs) {
            list(
                t.values = t.values,
                GSADF.value = GSADF.value,
                GSADF.bootstsrap.values = GSADF.bootstsrap.values,
                GSADF.cr.value = GSADF.cr.value,
                U.value = U.value,
                U.bootstsrap.values = U.bootstsrap.values,
                U.cr.value = U.cr.value
            )
        } else {
            NULL
        }
    )

    class(result) <- "sadf"

    return(result)
}
