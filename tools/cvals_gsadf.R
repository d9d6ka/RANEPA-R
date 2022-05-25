library(parallel)
library(doSNOW)
library(breaktest)

N_sim <- 999
N_obs <- c(
    seq(from = 30, to = 99, by = 10),
    seq(from = 100, to = 199, by = 20),
    seq(from = 200, to = 399, by = 50),
    seq(from = 400, to = 799, by = 100),
    seq(from = 800, to = 1599, by = 200),
    seq(from = 1600, to = 3200, by = 400)
)

## 30-100 шаг 10 наблюдений
## 100-200 шаг 20 наблюдений
## 200-400 шаг 50 наблюдений
## 400-800 шаг 100 наблюдений
## 800-1600 шаг 200 наблюдений
## 1600-3200 шаг 400 (или сразу 800 бахнуть) наблюдений

progress.bar <- txtProgressBar(max = N_sim, style = 3)
progress <- function(n) setTxtProgressBar(progress.bar, n)

cores <- detectCores()
cluster <- makeCluster(max(cores - 1, 1))
clusterExport(cluster, c("STADF.test", "GSTADF.test"))
registerDoSNOW(cluster)

.cval_SADF_without_const <- list()
.cval_SADF_with_const <- list()
.cval_GSADF_without_const <- list()
.cval_GSADF_with_const <- list()

for (step in N_obs) {
    cat("\nCalculating tables for", step, "observations\n")

    tmp <- NULL
    tmp <- foreach(
        i = 1:N_sim,
        .combine = c,
        .options.snow = list(progress = progress)
    ) %dopar% {
        y <- rnorm(step)
        y <- cumsum(y)

        model <- STADF.test(y, const = FALSE, add.p.value = FALSE)

        model$STADF.value
    }
    tmp <- sort(tmp)
    names(tmp) <- NULL
    .cval_SADF_without_const[[as.character(step)]] <- tmp

    tmp <- NULL
    tmp <- foreach(
        i = 1:N_sim,
        .combine = c,
        .options.snow = list(progress = progress)
    ) %dopar% {
        y <- rnorm(step)
        y <- cumsum(y)

        model <- STADF.test(y, const = TRUE, add.p.value = FALSE)

        model$STADF.value
    }
    tmp <- sort(tmp)
    names(tmp) <- NULL
    .cval_SADF_with_const[[as.character(step)]] <- tmp

    tmp <- NULL
    tmp <- foreach(
        i = 1:N_sim,
        .combine = c,
        .options.snow = list(progress = progress)
    ) %dopar% {
        y <- rnorm(step)
        y <- cumsum(y)

        model <- GSTADF.test(y, const = FALSE, add.p.value = FALSE)

        model$GSTADF.value
    }
    tmp <- sort(tmp)
    names(tmp) <- NULL
    .cval_GSADF_without_const[[as.character(step)]] <- tmp

    tmp <- NULL
    tmp <- foreach(
        i = 1:N_sim,
        .combine = c,
        .options.snow = list(progress = progress)
    ) %dopar% {
        y <- rnorm(step)
        y <- cumsum(y)

        model <- GSTADF.test(y, const = TRUE, add.p.value = FALSE)

        model$GSTADF.value
    }
    tmp <- sort(tmp)
    names(tmp) <- NULL
    .cval_GSADF_with_const[[as.character(step)]] <- tmp
}

stopCluster(cluster)
