#' @title
#' Weighted SADF test (HLZ, 2018).
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#'
#' @export
weighted.SADF.test <- function(y,
                               r0 = 0.01 + 1.8 / sqrt(length(y)),
                               const = TRUE,
                               alpha = 0.05,
                               iter = 4 * 200,
                               urs = TRUE,
                               seed = round(10^4 * sd(y))) {
    N <- length(y)

    # Find supBZ.value.
    supBZ.model <- supBZ.statistic(y, r0)
    sigma.sq <- supBZ.model$sigma.sq
    BZ.values <- supBZ.model$BZ.values
    supBZ.value <- supBZ.model$supBZ.value

    # Auxiliary function for parallel.
    calc_sadf_supBZ.value <- function(cpu, y, r0, const, B, seed, urs, sigma_sq) {
        set.seed(seed + cpu)
        N <- length(y)
        result <- list()
        for (i in 1:B) {
            y_star <- cumsum(c(0, rnorm(N - 1) * diff(y)))
            if (urs == TRUE) {
                sadf.model <- SADF.test(y_star, r0, const)
                result$sadf_value[i] <- sadf.model$sadf_value
            }
            supBZ.model <- supBZ.statistic(y_star, r0, sigma_sq)
            result$supBZ.value[i] <- supBZ.model$supBZ.value
        }
        return(result)
    }

    # Do parallel.
    cores <- detectCores()

    progress.bar <- txtProgressBar(max = iter, style = 3)
    progress <- function(n) setTxtProgressBar(progress.bar, n)

    cluster <- makeCluster(max(cores - 1, 1))
    clusterExport(cluster, c("ADF.test",
                             "SADF.test",
                             "supBZ.statistic",
                             ".cval_SADF_without_const",
                             ".cval_SADF_with_const"))
    registerDoSNOW(cluster)



    cl <- makePSOCKcluster(rep("localhost", cpu))
    clusterExport(cl, c(
        "SADF.test", "ADF.test", "supBZ.statistic",
        "SADF_cr_values_with_const",
        "SADF_cr_values_without_const"
    ))
    sadf_supBZ_bootstrap_values <- parLapply(
        cl, 1:cpu, calc_sadf_supBZ.value, y,
        r0, const, B, seed, urs, sigma.sq
    )
    stopCluster(cl)

    # Get sadf_supBZ_bootstrap_values.
    result$supBZ_bootstrap_values <- c()
    for (i in 1:cpu) {
        result$supBZ_bootstrap_values <- c(
            result$supBZ_bootstrap_values,
            sadf_supBZ_bootstrap_values[[i]]$supBZ.value
        )
    }

    # Find critical value.
    result$supBZ_cr_value <- as.numeric(quantile(
        result$supBZ_bootstrap_values,
        1 - alpha
    ))

    # A union of rejections strategy.
    if (urs == TRUE) {

        # Find sadf_value.
        sadf.model <- SADF.test(y, r0, const)
        result$t.values <- sadf.model$t.values
        result$sadf_value <- sadf.model$sadf_value

        # Get sadf_supBZ_bootstrap_values.
        result$sadf_bootstrap_values <- c()
        for (i in 1:cpu) {
            result$sadf_bootstrap_values <- c(
                result$sadf_bootstrap_values,
                sadf_supBZ_bootstrap_values[[i]]$sadf_value
            )
        }

        # Find critical value.
        result$sadf_cr_value <- as.numeric(quantile(
            result$sadf_bootstrap_values,
            1 - alpha
        ))

        # Calculate U value.
        result$U_value <- max(
            result$sadf_value,
            result$sadf_cr_value / result$supBZ_cr_value *
                result$supBZ.value
        )

        # Find U_bootstrap_values.
        result$U_bootstrap_values <- c()
        for (b in 1:(cpu * B)) {
            result$U_bootstrap_values[b] <- max(
                result$sadf_bootstrap_values[b],
                result$sadf_cr_value /
                    result$supBZ_cr_value *
                    result$supBZ_bootstrap_values[b]
            )
        }

        # Find critical value.
        result$U_cr_value <- as.numeric(quantile(result$U_bootstrap_values, 1 - alpha))

        result$p_value <- round(sum(result$U_bootstrap_values >
            result$U_value) / (cpu * B), 4)

        result$is_explosive <- ifelse(result$U_value > result$U_cr_value, 1, 0)
    } else {
        result$p_value <- round(sum(result$supBZ_bootstrap_values >
            result$supBZ.value) / (cpu * B), 4)

        result$is_explosive <- ifelse(result$supBZ.value > result$supBZ_cr_value, 1, 0)
    }

    return(
        list(
            y = y,
            r0 = r0,
            const = const,
            alpha = alpha,
            iter = iter,
            urs = urs,
            seed = seed,
            sigma.sq = sigma.sq,
            BZ.values = BZ.values,
            supBZ.value = supBZ.value,

        )
    )
}
