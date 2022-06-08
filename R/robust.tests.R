#'
#' @importFrom zeallot %<-%
#'
#' @export
robust.tests.1.break <- function(y,
                                 const = FALSE, trend = FALSE, season = FALSE,
                                 trim = 0.15) {
    if (!is.matrix(y)) y <- as.matrix(y)

    result <- list(
        const = const,
        trend = trend,
        season = season
    )

    ## Critical values
    cv.DF.GLS.m <- -1.93078
    cv.DF.OLS.m <- -2.85706

    cv.DF.GLS.t <- -2.84317
	cv.MDF.GLS <- -3.84632
	cv.DF.OLS.t <- -3.39735
	cv.MDF.OLS <- -4.23
	cv.MDF.t <- -4.25

    if (const && trend) {
        cv.MDF.OLS <- -4.21172
        cv.MDF.t <- -4.34468
    }

    if (!const && trend)
        cv.HLT <- 2.563
    else
        cv.HLT <- 3.162

	sap.ur <- 1.1009;
	sap.ur1 <- 1.0514;
	sap.ur2 <- 1.0317;
	sap.ur3 <- 1.0364;

	sap.cv.ur.HLT.k0 <- 1.0154
	sap.cv.A.HLT.k0 <- 1.00

	sap.cv.ur.PY.k0 <- 1.0159
	sap.cv.A.PY.k0 <- 1.00

    ## Start ##
    N <- nrow(y)

    if (season) {
        SEAS <- cbind(
            rep(1, N),
            seasonal.dummies(N)
        )
        c(., y, ., .) %<-% OLS(y, SEAS)
    }

    max.lag <- trunc(12 * (N / 100)^(1 / 4))

    first.break <- trunc(trim * N)
    last.break <- trunc((1 - trim) * N)

    tb <- segments.GLS(
        y, const, trend, 1,
        first.break, last.break,
        trim
    )
    result$break.time <- tb

    tau <- tb / N
    cv.MDF.GLS.lib = cv.MDF.GLS ## CV(tau, cv_MDF_GLS_lib_1,trm);
    cv.MDF.OLS.lib = cv.MDF.OLS ## CV(tau, cv_MDF_OLS_lib_1,trm);

    DU <- c(rep(0, tb), rep(1, N - tb))
    DT <- DU * (1:N - tb)

    x <- cbind(
        rep(1, N),
        1:N,
        if (const) DU else NULL,
        if (trend) DT else NULL
    )

    ## OLS/GLS ##
    ## Mean case
    c(beta.OLS.m, resid.OLS.m, ., .) %<-% OLS(y, x[, 1, drop = FALSE])
    DF.OLS.m <- ADF.test(resid.OLS.m,
                         const = FALSE, trend = FALSE,
                         max.lag = max.lag,
                         criterion = "aic",
                         modified.criterion = TRUE)
    k.m <- max(1, DF.OLS.m$lag)

    DF.OLS.m <- ADF.test(resid.OLS.m,
                         const = FALSE, trend = FALSE,
                         max.lag = k.m,
                         criterion = NULL)

    c(., resid.GLS.m, ., .) %<-% GLS(y, x[, 1, drop = FALSE], -7)
    DF.GLS.m <- ADF.test(resid.GLS.m,
                         const = FALSE, trend = FALSE,
                         max.lag = k.m,
                         criterion = NULL)
    denom.m <- 1 - sum(DF.GLS.m$beta) + DF.GLS.m$alpha

    ## Trend case
    c(beta.OLS.t, resid.OLS.t, ., .) %<-% OLS(y, x[, 1:2])
    DF.OLS.t <- ADF.test(resid.OLS.t,
                         const = FALSE, trend = FALSE,
                         max.lag = max.lag,
                         criterion = "aic",
                         modified.criterion = TRUE)
    k.t <- max(1, DF.OLS.t$lag)

    DF.OLS.t <- ADF.test(resid.OLS.t,
                         const = FALSE, trend = FALSE,
                         max.lag = k.t,
                         criterion = NULL)

    c(., resid.GLS.t, ., .) %<-% GLS(y, x[, 1:2], -13.5)
    DF.GLS.t <- ADF.test(resid.GLS.t,
                         const = FALSE, trend = FALSE,
                         max.lag = k.t,
                         criterion = NULL)
    denom.t <- 1 - sum(DF.GLS.t$beta) + DF.GLS.t$alpha

    ## ADF-OLS (lambda) ##
    c(beta.OLS, resid.OLS, ., .) %<-% OLS(y, x)
    DF1 <- ADF.test(resid.OLS,
                    const = FALSE, trend = FALSE,
                    max.lag = max.lag,
                    criterion = "aic",
                    modified.criterion = TRUE)
    k.tb <- max(1, DF1$lag)
    DF1 <- ADF.test(resid.OLS,
                    const = FALSE, trend = FALSE,
                    max.lag = k.tb,
                    criterion = NULL)
    MDF.OLS <- DF1$t.alpha

    ## MDF ##
    ## One break ##
    MDF.GLS <- Inf
    for (tb1 in first.break:last.break) {
        DU1 <- c(rep(0, tb1), rep(1, N - tb1))
        DT1 <- DU1 * (1:N - tb1)

        z <- cbind(
            rep(1, N),
            if (const) DU1 else NULL,
            1:N,
            if (trend) DT1 else NULL
        )

        c(beta.OLS, resid.OLS, ., .) %<-% OLS(y, z)
        DF1.tb <- ADF.test(resid.OLS,
                           const = FALSE, trend = FALSE,
                           max.lag = max.lag,
                           criterion = "aic",
                           modified.criterion = TRUE)
        k.tb <- max(1, DF1.tb$lag)

        c(., resid.GLS, ., .) %<-% GLS(y, z, -17.6)
        DF1.tb <- ADF.test(resid.GLS,
                           const = FALSE, trend = FALSE,
                           max.lag = k.tb,
                           criterion = NULL)

        if (DF1.tb$t.alpha < MDF.GLS) MDF.GLS <- DF1.tb$t.alpha
    }

    MDF.t <- Inf
    for (tb1 in first.break:last.break) {
        DU1 <- c(rep(0, tb1), rep(1, N - tb1))
        DT1 <- DU1 * (1:N - tb1)

        z <- cbind(
            rep(1, N),
            1:N,
            DT1
        )

        c(beta.OLS, resid.OLS, ., .) %<-% OLS(y, z)
        DF2 <- ADF.test(resid.OLS,
                        const = FALSE, trend = FALSE,
                        max.lag = max.lag,
                        criterion = "aic",
                        modified.criterion = TRUE)
        k.tb <- max(1, DF2$lag)
        DF2 <- ADF.test(resid.OLS,
                        const = FALSE, trend = FALSE,
                        max.lag = k.tb,
                        criterion = NULL)
        DF2.tb <- N * DF2$alpha / (1 - sum(DF2$beta) + DF2$alpha)

#########TODO: Проверить!!!
        if (DF2$t.alpha < MDF.t) MDF.t <- DF2$t.alpha
    }


    t.HLT <- KPSS.HLT(y, const, trim)
    c(t.PY, cv.PY) %<-% PY(y, const, trend, "aic", trim, max.lag)

    t.alpha <- beta.OLS[1] / sqrt(drop(t(resid.OLS) %*% resid.OLS) / N)
    t.alpha.id = as.numeric(abs(t.alpha) > 1)

    ## UR-HLT
    t.lambda.id = as.numeric(t.HLT > cv.HLT)
    ur.id.sa <- as.numeric(
    (DF.GLS.t$t.alpha < (sap.ur * sap.cv.ur.HLT.k0 * cv.DF.OLS.t)) ||
    (MDF.GLS < (sap.ur * sap.cv.ur.HLT.k0 * cv.MDF.GLS)) ||
    (DF.OLS.t$t.alpha < (sap.ur * sap.cv.ur.HLT.k0 * cv.DF.OLS.t)) ||
    (MDF.t < (sap.ur * sap.cv.ur.HLT.k0 * cv.MDF.t))
    )
    ur1.id.sa <- as.numeric(
    (DF.OLS.t$t.alpha < (sap.ur1 * sap.cv.ur.HLT.k0 * cv.DF.OLS.t)) ||
    (MDF.t < (sap.ur1 * sap.cv.ur.HLT.k0 * cv.MDF.t))
    )
    ur2.id.sa <- as.numeric(
    (MDF.GLS < (sap.ur2 * sap.cv.ur.HLT.k0 * cv.MDF.GLS.lib)) ||
    (MDF.OLS < (sap.ur2 * sap.cv.ur.HLT.k0 * cv.MDF.OLS.lib))
    )
    ur3.id.sa <- as.numeric(MDF.OLS < sap.cv.ur.HLT.k0 * cv.MDF.OLS.lib)
    UR.HLT <-
        (1 - t.lambda.id) * (1 - t.alpha.id) * ur.id.sa +
        (1 - t.lambda.id) * t.alpha.id * ur1.id.sa +
        t.lambda.id * (1 - t.alpha.id) * ur2.id.sa +
        t.lambda.id * t.alpha.id * ur3.id.sa

    ## UR-PY
    t.lambda.id = as.numeric(t.PY > cv.PY[2])
    ur.id.sa <- as.numeric(
    (DF.GLS.t$t.alpha < (sap.ur * sap.cv.ur.PY.k0 * cv.DF.GLS.t)) ||
    (MDF.GLS < (sap.ur * sap.cv.ur.PY.k0 * cv.MDF.GLS)) ||
    (DF.OLS.t$t.alpha < (sap.ur * sap.cv.ur.PY.k0 * cv.DF.OLS.t)) ||
    (MDF.t < (sap.ur * sap.cv.ur.PY.k0 * cv.MDF.t))
    )
    ur1.id.sa <- as.numeric(
    (DF.OLS.t$t.alpha < (sap.ur1 * sap.cv.ur.PY.k0 * cv.DF.OLS.t)) ||
    (MDF.t < (sap.ur1 * sap.cv.ur.PY.k0 * cv.MDF.t))
    )
    ur2.id.sa <- as.numeric(
    (MDF.GLS < (sap.ur2 * sap.cv.ur.PY.k0 * cv.MDF.GLS.lib)) ||
    (MDF.OLS < (sap.ur2 * sap.cv.ur.PY.k0 * cv.MDF.OLS.lib))
    )
    ur3.id.sa <- as.numeric(MDF.OLS < sap.cv.ur.PY.k0 * cv.MDF.OLS.lib)
    UR.PY <-
        (1 - t.lambda.id) * (1 - t.alpha.id) * ur.id.sa +
        (1 - t.lambda.id) * t.alpha.id * ur1.id.sa +
        t.lambda.id * (1 - t.alpha.id) * ur2.id.sa +
        t.lambda.id * t.alpha.id * ur3.id.sa

    ## A-HLT
    t.lambda.id <- as.numeric(t.HLT > cv.HLT)
    ur2.id.sa <- as.numeric(
    (MDF.GLS < (sap.ur2 * sap.cv.A.HLT.k0 * cv.MDF.GLS.lib)) ||
    (MDF.OLS < (sap.ur2 * sap.cv.A.HLT.k0 * cv.MDF.OLS.lib))
    )
    ur2a.id.sa <- as.numeric(MDF.OLS < sap.cv.A.HLT.k0 * cv.MDF.OLS.lib)
    ur3.id.sa <- as.numeric(
    (MDF.GLS < (sap.ur3 * sap.cv.A.HLT.k0 * cv.MDF.GLS)) ||
    (MDF.t < (sap.ur3 * sap.cv.A.HLT.k0 * cv.MDF.t))
    )
    ur3a.id.sa <- as.numeric(MDF.t < sap.cv.A.HLT.k0 * cv.MDF.t)

    A.HLT <-
        (1 - t.lambda.id) * (1 - t.alpha.id) * ur3.id.sa +
        (1 - t.lambda.id) * t.alpha.id * ur3a.id.sa +
        t.lambda.id * (1 - t.alpha.id) * ur2.id.sa +
        t.lambda.id * t.alpha.id * ur2a.id.sa

    ## A-PY
    t.lambda.id <- as.numeric(t.PY > cv.PY[2])
    ur2.id.sa <- as.numeric(
    (MDF.GLS < (sap.ur2 * sap.cv.A.PY.k0 * cv.MDF.GLS.lib)) ||
    (MDF.OLS < (sap.ur2 * sap.cv.A.PY.k0 * cv.MDF.OLS.lib))
    )
    ur2a.id.sa <- as.numeric(MDF.OLS < sap.cv.A.PY.k0 * cv.MDF.OLS.lib)
    ur3.id.sa <- as.numeric(
    (MDF.GLS < (sap.ur3 * sap.cv.A.PY.k0 * cv.MDF.GLS)) ||
    (MDF.t < (sap.ur3 * sap.cv.A.PY.k0 * cv.MDF.t))
    )
    ur3.id.sa <- as.numeric(MDF.t < sap.cv.A.PY.k0 * cv.MDF.t)

    A.PY <- (1 - t.lambda.id) * (1 - t.alpha.id) * ur3.id.sa +
        (1 - t.lambda.id) * t.alpha.id * ur3.id.sa +
        t.lambda.id * (1 - t.alpha.id) * ur2.id.sa +
        t.lambda.id * t.alpha.id * ur2a.id.sa

    result$HLT <- list(
        stat = t.HLT,
        cv = cv.HLT
    )
    result$PY <- list(
        stat = t.PY,
        cv = cv.PY[2]
    )
    result$DF.GLS <- list(
        stat = DF.GLS.t$t.alpha,
        cv = cv.DF.GLS.t
    )
    result$DF.OLS <- list(
        stat = DF.OLS.t$t.alpha,
        cv = cv.DF.OLS.t
    )
    result$MDF.GLS <- list(
        stat = MDF.GLS,
        cv = cv.MDF.GLS
    )
    result$MDF.OLS <- list(
        stat = MDF.OLS,
        cv = cv.MDF.OLS
    )
    result$MDF.t <- list(
        stat = MDF.t,
        cv = cv.MDF.t
    )
    result$A.HLT <- A.HLT
    result$A.PY <- A.PY
    result$UR.HLT <- UR.HLT
    result$UR.PY <- UR.PY

    class(result) <- "robustUR"

    return(result)
}
