library(parallel)
library(doSNOW)
library(breaktest)

NN_sim <- 9999

progress.bar <- txtProgressBar(max = NN_sim, style = 3)
progress <- function(n) setTxtProgressBar(progress.bar, n)

cores <- detectCores()
cluster <- makeCluster(max(cores - 1, 1))
clusterExport(cluster, c("STADF.test", "GSTADF.test"))
registerDoSNOW(cluster)

cat(" --------------------- 100 --------------------- ")
# tau.break - Bubble period.
cval_SADF_without_const_100 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(100)
    y <- cumsum(y)

    model <- STADF.test(y, const = FALSE, add.p.value = FALSE)

    model$stadf.value
}
cval_SADF_without_const_100 <- sort(cval_SADF_without_const_100)
names(cval_SADF_without_const_100) <- NULL

cval_SADF_with_const_100 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(100)
    y <- cumsum(y)

    model <- STADF.test(y, const = TRUE, add.p.value = FALSE)

    model$stadf.value
}
cval_SADF_with_const_100 <- sort(cval_SADF_with_const_100)
names(cval_SADF_with_const_100) <- NULL

# tau.break - Bubble period.
cval_GSADF_without_const_100 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(100)
    y <- cumsum(y)

    model <- GSTADF.test(y, const = FALSE, add.p.value = FALSE)

    model$gstadf.value
}
cval_GSADF_without_const_100 <- sort(cval_GSADF_without_const_100)
names(cval_GSADF_without_const_100) <- NULL

cval_GSADF_with_const_100 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(100)
    y <- cumsum(y)

    model <- GSTADF.test(y, const = TRUE, add.p.value = FALSE)

    model$gstadf.value
}
cval_GSADF_with_const_100 <- sort(cval_GSADF_with_const_100)
names(cval_GSADF_with_const_100) <- NULL

cat(" --------------------- 200 --------------------- ")
# tau.break - Bubble period.
cval_SADF_without_const_200 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(200)
    y <- cumsum(y)

    model <- STADF.test(y, const = FALSE, add.p.value = FALSE)

    model$stadf.value
}
cval_SADF_without_const_200 <- sort(cval_SADF_without_const_200)
names(cval_SADF_without_const_200) <- NULL

cval_SADF_with_const_200 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(200)
    y <- cumsum(y)

    model <- STADF.test(y, const = TRUE, add.p.value = FALSE)

    model$stadf.value
}
cval_SADF_with_const_200 <- sort(cval_SADF_with_const_200)
names(cval_SADF_with_const_200) <- NULL

# tau.break - Bubble period.
cval_GSADF_without_const_200 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(200)
    y <- cumsum(y)

    model <- GSTADF.test(y, const = FALSE, add.p.value = FALSE)

    model$gstadf.value
}
cval_GSADF_without_const_200 <- sort(cval_GSADF_without_const_200)
names(cval_GSADF_without_const_200) <- NULL

cval_GSADF_with_const_200 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(200)
    y <- cumsum(y)

    model <- GSTADF.test(y, const = TRUE, add.p.value = FALSE)

    model$gstadf.value
}
cval_GSADF_with_const_200 <- sort(cval_GSADF_with_const_200)
names(cval_GSADF_with_const_200) <- NULL

cat(" --------------------- 400 --------------------- ")
# tau.break - Bubble period.
cval_SADF_without_const_400 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(400)
    y <- cumsum(y)

    model <- STADF.test(y, const = FALSE, add.p.value = FALSE)

    model$stadf.value
}
cval_SADF_without_const_400 <- sort(cval_SADF_without_const_400)
names(cval_SADF_without_const_400) <- NULL

cval_SADF_with_const_400 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(400)
    y <- cumsum(y)

    model <- STADF.test(y, const = TRUE, add.p.value = FALSE)

    model$stadf.value
}
cval_SADF_with_const_400 <- sort(cval_SADF_with_const_400)
names(cval_SADF_with_const_400) <- NULL

# tau.break - Bubble period.
cval_GSADF_without_const_400 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(400)
    y <- cumsum(y)

    model <- GSTADF.test(y, const = FALSE, add.p.value = FALSE)

    model$gstadf.value
}
cval_GSADF_without_const_400 <- sort(cval_GSADF_without_const_400)
names(cval_GSADF_without_const_400) <- NULL

cval_GSADF_with_const_400 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(400)
    y <- cumsum(y)

    model <- GSTADF.test(y, const = TRUE, add.p.value = FALSE)

    model$gstadf.value
}
cval_GSADF_with_const_400 <- sort(cval_GSADF_with_const_400)
names(cval_GSADF_with_const_400) <- NULL

cat(" --------------------- 800 --------------------- ")
# tau.break - Bubble period.
cval_SADF_without_const_800 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(800)
    y <- cumsum(y)

    model <- STADF.test(y, const = FALSE, add.p.value = FALSE)

    model$stadf.value
}
cval_SADF_without_const_800 <- sort(cval_SADF_without_const_800)
names(cval_SADF_without_const_800) <- NULL

cval_SADF_with_const_800 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(800)
    y <- cumsum(y)

    model <- STADF.test(y, const = TRUE, add.p.value = FALSE)

    model$stadf.value
}
cval_SADF_with_const_800 <- sort(cval_SADF_with_const_800)
names(cval_SADF_with_const_800) <- NULL

# tau.break - Bubble period.
cval_GSADF_without_const_800 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(800)
    y <- cumsum(y)

    model <- GSTADF.test(y, const = FALSE, add.p.value = FALSE)

    model$gstadf.value
}
cval_GSADF_without_const_800 <- sort(cval_GSADF_without_const_800)
names(cval_GSADF_without_const_800) <- NULL

cval_GSADF_with_const_800 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(800)
    y <- cumsum(y)

    model <- GSTADF.test(y, const = TRUE, add.p.value = FALSE)

    model$gstadf.value
}
cval_GSADF_with_const_800 <- sort(cval_GSADF_with_const_800)
names(cval_GSADF_with_const_800) <- NULL

cat(" --------------------- 1600 --------------------- ")
# tau.break - Bubble period.
cval_SADF_without_const_1600 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(1600)
    y <- cumsum(y)

    model <- STADF.test(y, const = FALSE, add.p.value = FALSE)

    model$stadf.value
}
cval_SADF_without_const_1600 <- sort(cval_SADF_without_const_1600)
names(cval_SADF_without_const_1600) <- NULL

cval_SADF_with_const_1600 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(1600)
    y <- cumsum(y)

    model <- STADF.test(y, const = TRUE, add.p.value = FALSE)

    model$stadf.value
}
cval_SADF_with_const_1600 <- sort(cval_SADF_with_const_1600)
names(cval_SADF_with_const_1600) <- NULL

# tau.break - Bubble period.
cval_GSADF_without_const_1600 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(1600)
    y <- cumsum(y)

    model <- GSTADF.test(y, const = FALSE, add.p.value = FALSE)

    model$gstadf.value
}
cval_GSADF_without_const_1600 <- sort(cval_GSADF_without_const_1600)
names(cval_GSADF_without_const_1600) <- NULL

cval_GSADF_with_const_1600 <- foreach(
    i = 1:NN_sim,
    .combine = c,
    .options.snow = list(progress = progress)
) %dopar% {
    y <- rnorm(1600)
    y <- cumsum(y)

    model <- GSTADF.test(y, const = TRUE, add.p.value = FALSE)

    model$gstadf.value
}
cval_GSADF_with_const_1600 <- sort(cval_GSADF_with_const_1600)
names(cval_GSADF_with_const_1600) <- NULL

.cval_SADF_with_const <- list(
    cval_SADF_with_const_100,
    cval_SADF_with_const_200,
    cval_SADF_with_const_400,
    cval_SADF_with_const_800,
    cval_SADF_with_const_1600
)

.cval_SADF_without_const <- list(
    cval_SADF_without_const_100,
    cval_SADF_without_const_200,
    cval_SADF_without_const_400,
    cval_SADF_without_const_800,
    cval_SADF_without_const_1600
)

.cval_GSADF_with_const <- list(
    cval_GSADF_with_const_100,
    cval_GSADF_with_const_200,
    cval_GSADF_with_const_400,
    cval_GSADF_with_const_800,
    cval_GSADF_with_const_1600
)

.cval_GSADF_without_const <- list(
    cval_GSADF_without_const_100,
    cval_GSADF_without_const_200,
    cval_GSADF_without_const_400,
    cval_GSADF_without_const_800,
    cval_GSADF_without_const_1600
)

stopCluster(cluster)
