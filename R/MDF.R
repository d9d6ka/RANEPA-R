#' @importFrom zeallot %<-%
#'
#' @export
MDF <- function(y,
                const = FALSE,
                breaks = 1,
                breaks.star = 1,
                trim = 0.15,
                ZA = FALSE) {
    if (!is.matrix(y)) y <- as.matrix(y)

    ## Critical values ##
    cv.DF.GLS.t <- -2.84317
    cv.DF.OLS.t <- -3.39735
     cv.MDF.GLS1 <- -3.84632
    cv.MDF.OLS1 <- -34.7163
     cv.MDF.GLS2 <- -4.55807
    cv.MDF.OLS2 <- -47.2304
     cv.MDF.GLS3 <- -4.84651
    cv.MDF.OLS3 <- -52.1377

    sap.ur.3ols <- 1.1382
    sap.ur.3olsgls <- 1.1598
    sap.ur3.1ols <- 1.1133
    sap.ur3.1olsgls <- 1.1178
    sap.ur3.2ols <- 1.0548
    sap.ur3.2olsgls <- 1.066
    sap.ur3.3ols <- 1
    sap.ur3.3olsgls <- 1.0301

    sap.ur2ols <- 1.1137
    sap.ur2olsgls <- 1.1402
    sap.ur2.1ols <- 1.0764
    sap.ur2.1olsgls <- 1.091
    sap.ur2.2ols <- 1
    sap.ur2.2olsgls <- 1.0312

    sap.cv.ur.2 <- 1.0115
    sap.cv.ur.3 <- 1.0053

    if (ZA) {
        sap.cv.ur.2 <- 1.0053
        sap.cv.ur.3 <- 1.0115
    }

    if (const && !ZA) {
        cv.MDF.OLS1 <- -40.4767
        cv.MDF.OLS2 <- -57.5673
        cv.MDF.OLS3 <- -67.6389

        sap.ur3.ols <- 1.1176
        sap.ur3.olsgls <- 1.1505
        sap.ur3.1ols <- 1.091
        sap.ur3.1olsgls <- 1.1137
        sap.ur3.2ols <- 1.0478
        sap.ur3.2olsgls <- 1.078
        sap.ur3.3ols <- 1
        sap.ur3.3olsgls <- 1.0447

        sap.ur2.ols <- 1.0987
        sap.ur2.olsgls <- 1.1323
        sap.ur2.1ols <- 1.0539
        sap.ur2.1olsgls <- 1.0918
        sap.ur2.2ols <- 1
        sap.ur2.2olsgls <- 1.0439

        sap.cv.ur.2 <- 1.0075
        sap.cv.ur.3 <- 1.0053
    }

    if (const && ZA) {
        cv.MDF.OLS1 <- -4.55413
        cv.MDF.OLS2 <- -5.41121
        cv.MDF.OLS3 <- -5.84884

        sap.ur3.ols <- 1.075
        sap.ur3.olsgls <- 1.12
        sap.ur3.1ols <- 1.0416
        sap.ur3.1olsgls <- 1.0865
        sap.ur3.2ols <- 1.0253
        sap.ur3.2olsgls <- 1.0615
        sap.ur3.3ols <- 1
        sap.ur3.3olsgls <- 1.0337

        sap.ur2.ols <- 1.0625
        sap.ur2.olsgls <- 1.1105
        sap.ur2.1ols <- 1.0297
        sap.ur2.1olsgls <- 1.0732
        sap.ur2.2ols <- 1
        sap.ur2.2olsgls <- 1.0352

        sap.cv.ur.2 <- 1.0126
        sap.cv.ur.3 <- 1.0098
    }

    N <- nrow(y)

    max.lag <- trunc(12 * (N / 100)^(1 / 4))

    x.const <- rep(1, N)
    x.trend <- 1:N

    first.break <- trunc(trim * N) + 1
    last.break <- trunc((1 - trim) * N) + 1

    x <- cbind(x.const, x.trend)

    ## GLS case
    c(., resid.GLS.t, ., .) %<-% GLS(y, x, -13.5)

    c(beta.OLS.t, resid.OLS.t, ., .) %<-% OLS(y, x)
    DF.OLS.t <- ADF.test(resid.OLS.t,
                         const = FALSE, trend = FALSE,
                         max.lag = max.lag,
                         criterion = "aic",
                         modified.criterion = TRUE)
    k.m <- max(1, DF.OLS.m$lag)

    DF.GLS.t <- ADF.test(resid.GLS.t,
                         const = FALSE, trend = FALSE,
                         max.lag = k.t,
                         criterion = NULL)
    DF.GLS.t <- DF.GLS.t$t.alpha

    ## OLS case
    DF.OLS.t <- ADF.test(resid.OLS.t,
                         const = FALSE, trend = FALSE,
                         max.lag = k.t,
                         criterion = NULL)
    DF.OLS.t <- DF.OLS.t$t.alpha

    ## OLS-GLS ##
    ## One break
    MDF.OLS1 <- Inf
    MDF.GLS1 <- Inf

    for (tb1 in first.break:last.break) {
        DU1 <- as.numeric(x.trend > tb1)
		DT1 <- DU1 * (x.trend - tb1)

        x <- cbind(
            x.const,
            x.trend,
            if (const) DU1 else NULL,
            DT1
        )

        c(beta.OLS, resid.OLS, ., .) %<-% OLS(y, x)
        DF1.tb <- ADF.test(resid.OLS,
                           const = FALSE, trend = FALSE,
                           max.lag = max.lag,
                           criterion = "aic",
                           modified.criterion = TRUE)
        k.m <- max(1, DF1.tb$lag)

        DF1.tb <- ADF.test(resid.OLS,
                           const = FALSE, trend = FALSE,
                           max.lag = k.t,
                           criterion = NULL)
        if (!ZA) {
            denom <- 1 - sum(DF1.tb$beta) + DF1.tb$alpha
            stat.OLS <- N * DF1.tb$alpha / denom
        } else {
            stat.OLS <- DF1.tb$t.alpha
        }

        c(., resid.GLS, ., .) %<-% GLS(y, x, -17.6)
        DF1.tb <- ADF.test(resid.GLS,
                           const = FALSE, trend = FALSE,
                           max.lag = k.t,
                           criterion = NULL)

        if (stat.OLS < MDF.OLS1) MDF.OLS1 <- stat.OLS
        if (DF1.tb$t.alpha < MDF.GLS1) MDF.GLS1 <- DF1.tb$t.alpha
    }

    ## Two breaks
    MDF.OLS2 <- Inf
    MDF.GLS2 <- Inf

    for (tb1 in first.break:(last.break - first.break)) {
        for (tb2 in (first.break + first.break):last.break) {
            DU1 <- as.numeric(x.trend > tb1)
            DT1 <- DU1 * (x.trend - tb1)
            DU2 <- as.numeric(x.trend > tb2)
            DT2 <- DU2 * (x.trend - tb2)

            x <- cbind(
                x.const,
                x.trend,
                if (const) DU1 else NULL,
                DT1,
                if (const) DU2 else NULL,
                DT2
            )

            c(beta.OLS, resid.OLS, ., .) %<-% OLS(y, x)
            DF2.tb <- ADF.test(resid.OLS,
                const = FALSE, trend = FALSE,
                max.lag = max.lag,
                criterion = "aic",
                modified.criterion = TRUE
            )
            k.m <- max(1, DF2.tb$lag)

            DF2.tb <- ADF.test(resid.OLS,
                const = FALSE, trend = FALSE,
                max.lag = k.t,
                criterion = NULL
            )
            if (!ZA) {
                denom <- 1 - sum(DF2.tb$beta) + DF2.tb$alpha
                stat.OLS <- N * DF2.tb$alpha / denom
            } else {
                stat.OLS <- DF2.tb$t.alpha
            }

            c(., resid.GLS, ., .) %<-% GLS(y, x, -21.5)
            DF2.tb <- ADF.test(resid.GLS,
                const = FALSE, trend = FALSE,
                max.lag = k.t,
                criterion = NULL
            )

            if (stat.OLS < MDF.OLS2) MDF.OLS2 <- stat.OLS
            if (DF2.tb$t.alpha < MDF.GLS2) MDF.GLS2 <- DF2.tb$t.alpha
        }
    }

    ## Three breaks
    MDF.OLS3 <- Inf
    MDF.GLS3 <- Inf

    for (tb1 in first.break:(last.break - 2 * first.break)) {
        for (tb2 in (first.break + first.break):(last.break - first.break)) {
            for (tb3 in (first.break + 2 * first.break):last.break) {
                DU1 <- as.numeric(x.trend > tb1)
                DT1 <- DU1 * (x.trend - tb1)
                DU2 <- as.numeric(x.trend > tb2)
                DT2 <- DU2 * (x.trend - tb2)
                DU3 <- as.numeric(x.trend > tb3)
                DT3 <- DU3 * (x.trend - tb3)

                x <- cbind(
                    x.const,
                    x.trend,
                    if (const) DU1 else NULL,
                    DT1,
                    if (const) DU2 else NULL,
                    DT2,
                    if (const) DU3 else NULL,
                    DT3
                )

                c(beta.OLS, resid.OLS, ., .) %<-% OLS(y, x)
                DF3.tb <- ADF.test(resid.OLS,
                    const = FALSE, trend = FALSE,
                    max.lag = max.lag,
                    criterion = "aic",
                    modified.criterion = TRUE
                )
                k.m <- max(1, DF3.tb$lag)

                 DF3.tb <- ADF.test(resid.OLS,
                     const = FALSE, trend = FALSE,
                     max.lag = k.t,
                     criterion = NULL
                 )
                 if (!ZA) {
                     denom <- 1 - sum(DF3.tb$beta) + DF3.tb$alpha
                     stat.OLS <- N * DF3.tb$alpha / denom
                 } else {
                     stat.OLS <- DF3.tb$t.alpha
                 }

                c(., resid.GLS, ., .) %<-% GLS(y, x, -25.5)
                DF3.tb <- ADF.test(resid.GLS,
                    const = FALSE, trend = FALSE,
                    max.lag = k.t,
                    criterion = NULL
                )

                if (stat.OLS < MDF.OLS3) MDF.OLS3 <- stat.OLS
                if (DF3.tb$t.alpha < MDF.GLS3) MDF.GLS3 <- DF3.tb$t.alpha
            }
        }
    }

    ## Alternative break selection
    if (breaks == 2) {
        c(tb1, tb2) %<-% segments.GLS(y, const, TRUE, 2)

        DU1 <- as.numeric(x.trend > tb1)
        DT1 <- DU1 * (x.trend - tb1)
        DU2 <- as.numeric(x.trend > tb2)
        DT2 <- DU2 * (x.trend - tb2)

        x <- cbind(
            x.const,
            x.trend,
            if (const) DU1 else NULL,
            DT1,
            if (const) DU2 else NULL,
            DT2
        )

        c(bb, rr, ., .) %<-% OLS(y, x)
        t.alpha <- bb[1] / sqrt(drop(t(rr) %*% rr) / N)
        t.alpha.2.id <- as.numeric(t.alpha > 1)
    }
    if (breaks == 3) {
        c(tb1, tb2) %<-% segments.GLS(y, const, TRUE, 3)

        DU1 <- as.numeric(x.trend > tb1)
        DT1 <- DU1 * (x.trend - tb1)
        DU2 <- as.numeric(x.trend > tb2)
        DT2 <- DU2 * (x.trend - tb2)
        DU3 <- as.numeric(x.trend > tb3)
        DT3 <- DU3 * (x.trend - tb3)

        x <- cbind(
            x.const,
            x.trend,
            if (const) DU1 else NULL,
            DT1,
            if (const) DU2 else NULL,
            DT2,
            if (const) DU3 else NULL,
            DT3
        )

        c(bb, rr, ., .) %<-% OLS(y, x)
        t.alpha <- bb[1] / sqrt(drop(t(rr) %*% rr) / N)
        t.alpha.3.id <- as.numeric(t.alpha > 1)
    }

    ## breaks.star
    if (breaks.star == 0)
        tbb <- 0
    else
        tbb <- segments.GLS(y, const, TRUE, breaks.star)

    if (breaks == 2) {
        ur2.ols.sa <- as.numeric(
        (DF.OLS.t < (sap.cv.ur.2 * sap.ur2.ols * cv.DF.OLS.t)) ||
        (MDF.OLS1 < (sap.cv.ur.2 * sap.ur2.ols * cv.MDF.OLS1)) ||
        (MDF.OLS2 < (sap.cv.ur.2 * sap.ur2.ols * cv.MDF.OLS2))
        )

        ur2.olsgls.sa <- as.numeric(
        (DF.GLS.t < (sap.cv.ur.2 * sap.ur2.olsgls * cv.DF.GLS.t)) ||
        (MDF.GLS1 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.GLS1)) ||
        (DF.OLS.t < (sap.cv.ur.2 * sap.ur2.olsgls * cv.DF.OLS.t)) ||
        (MDF.OLS1 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.OLS1)) ||
        (MDF.GLS2 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.GLS2)) ||
        (MDF.OLS2 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.OLS2))
        )

        ur2.1ols.sa <- as.numeric(
        (MDF.OLS1 < (sap.cv.ur.2 * sap.ur2.1ols * cv.MDF.OLS1)) ||
        (MDF.OLS2 < (sap.cv.ur.2 * sap.ur2.1ols * cv.MDF.OLS2))
        )

        ur2.1olsgls.sa <- as.numeric(
        (MDF.GLS1 < (sap.cv.ur.2 * sap.ur2.1olsgls * cv.MDF.GLS1)) ||
        (MDF.OLS1 < (sap.cv.ur.2 * sap.ur2.1olsgls * cv.MDF.OLS1)) ||
        (MDF.GLS2 < (sap.cv.ur.2 * sap.ur2.1olsgls * cv.MDF.GLS2)) ||
        (MDF.OLS2 < (sap.cv.ur.2 * sap.ur2.1olsgls * cv.MDF.OLS2))
        )

        ur2.2ols.sa <- as.numeric(
            MDF.OLS2 < (sap.cv.ur.2 * sap.ur2.2ols * cv.MDF.OLS2)
        )

        ur2.2olsgls.sa <- as.numeric(
        (MDF.GLS2 < (sap.cv.ur.2 * sap.ur2.2olsgls * cv.MDF.GLS2)) ||
        (MDF.OLS2 < (sap.cv.ur.2 * sap.ur2.2olsgls * cv.MDF.OLS2))
        )

        UR <- as.numeric(m.star == 2) *
            ((1 - t.alpha.2.id) * ur2.olsgls.sa + t.alpha.2.id * ur2.ols.sa) +
            as.numeric(m.star == 1) *
            ((1 - t.alpha.2.id) * ur2.1olsgls.sa + t.alpha.2.id * ur2.1ols.sa) +
            as.numeric(m.star == 0) *
            ((1-t.alpha.2.id) * ur2.2olsgls.sa + t.alpha.2.id * ur2.2ols.sa)

        ## without pre-test	for breaks

        ur2.ols.sa <- as.numeric(
        (DF.OLS.t < (sap.cv.ur.2 * sap.ur2.ols * cv.DF.OLS.t)) ||
        (MDF.OLS1 < (sap.cv.ur.2 * sap.ur2.ols * cv.MDF.OLS1)) ||
        (MDF.OLS2 < (sap.cv.ur.2 * sap.ur2.ols * cv.MDF.OLS2))
        )

        ur2.olsgls.sa <- as.numeric(
        (DF.GLS.t < (sap.cv.ur.2 * sap.ur2.olsgls * cv.DF.GLS.t)) ||
        (MDF.GLS1 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.GLS1)) ||
        (DF.OLS.t < (sap.cv.ur.2 * sap.ur2.olsgls * cv.DF.OLS.t)) ||
        (MDF.OLS1 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.OLS1)) ||
        (MDF.GLS2 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.GLS2)) ||
        (MDF.OLS2 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.OLS2))
        )

        UR1 <- (1 - t.alpha.2.id) * ur2.olsgls.sa + t.alpha.2.id * ur2.ols.sa
    } else if (breaks == 3) {
        ur3.ols.sa <- as.numeric(
        (DF.OLS.t < (sap.cv.ur.3 * sap.ur3.ols * cv.DF.OLS.t)) ||
        (MDF.OLS1 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS1)) ||
        (MDF.OLS2 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS2)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS3))
        )

        ur3.olsgls.sa <- as.numeric(
        (DF.GLS.t < (sap.cv.ur.3 * sap.ur3.olsgls * cv.DF.GLS.t)) ||
        (MDF.GLS1 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS1)) ||
        (DF.OLS.t < (sap.cv.ur.3 * sap.ur3.olsgls * cv.DF.OLS.t)) ||
        (MDF.OLS1 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS1)) ||
        (MDF.GLS2 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS2)) ||
        (MDF.OLS2 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS2)) ||
        (MDF.GLS3 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS3)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS3))
        )

        ur3.1ols.sa <- as.numeric(
        (MDF.OLS1 < (sap.cv.ur.3 * sap.ur3.1ols * cv.MDF.OLS1)) ||
        (MDF.OLS2 < (sap.cv.ur.3 * sap.ur3.1ols * cv.MDF.OLS2)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.1ols * cv.MDF.OLS3))
        )

        ur3.1olsgls.sa <- as.numeric(
        (MDF.GLS1 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.GLS1)) ||
        (MDF.OLS1 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.OLS1)) ||
        (MDF.GLS2 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.GLS2)) ||
        (MDF.OLS2 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.OLS2)) ||
        (MDF.GLS3 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.GLS3)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.OLS3))
        )

        ur3.2ols.sa <- as.numeric(
        (MDF.OLS2 < (sap.cv.ur.3 * sap.ur3.2ols * cv.MDF.OLS2)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.2ols * cv.MDF.OLS3))
        )

        ur3.2olsgls.sa <- as.numeric(
        (MDF.GLS2 < (sap.cv.ur.3 * sap.ur3.2olsgls * cv.MDF.GLS2)) ||
        (MDF.OLS2 < (sap.cv.ur.3 * sap.ur3.2olsgls * cv.MDF.OLS2)) ||
        (MDF.GLS3 < (sap.cv.ur.3 * sap.ur3.2olsgls * cv.MDF.GLS3)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.2olsgls * cv.MDF.OLS3))
        )

        ur3.3ols.sa <- as.numeric(
            MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.3ols * cv.MDF.OLS3)
        )

        ur3.3olsgls.sa <- as.numeric(
        (MDF.GLS3 < (sap.cv.ur.3 * sap.ur3.3olsgls * cv.MDF.GLS3)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.3olsgls * cv.MDF.OLS3))
        )

        UR <- as.numeric(m.star == 3) *
            ((1 - t.alpha.3.id) * ur3.olsgls.sa + t.alpha.3.id * ur3.ols.sa) +
            as.numeric(m.star == 2) *
            ((1-t.alpha.3.id) * ur3.1olsgls.sa + t.alpha.3.id * ur3.1ols.sa) +
            as.numeric(m.star == 1) *
            ((1-t.alpha.3.id) * ur3.2olsgls.sa + t.alpha.3.id * ur3.2ols.sa) +
            as.numeric(m.star == 0) *
            ((1-t.alpha.3.id) * ur3.3olsgls.sa + t.alpha.3.id * ur3.3ols.sa)

        ## without pre-test for breaks

        ur3.ols.sa <- as.numeric(
        (DF.OLS.t < (sap.cv.ur.3 * sap.ur3.ols * cv.DF.OLS.t)) ||
        (MDF.OLS1 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS1)) ||
        (MDF.OLS2 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS2)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS3))
        )

        ur3.olsgls.sa <- as.numeric(
        (DF.GLS.t < (sap.cv.ur.3 * sap.ur3.olsgls * cv.DF.GLS.t)) ||
        (MDF.GLS1 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS1)) ||
        (DF.OLS.t < (sap.cv.ur.3 * sap.ur3.olsgls * cv.DF.OLS.t)) ||
        (MDF.OLS1 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS1)) ||
        (MDF.GLS2 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS2)) ||
        (MDF.OLS2 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS2)) ||
        (MDF.GLS3 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS3)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS3))
        )

        UR1 <- (1 - t.alpha.3.id) * ur3.olsgls.sa + t.alpha.3.id * ur3.ols.sa
    }

    result <- list()
    result$breaks.star <- breaks.star
    result$breaks.tbb <- tbb
    result$breaks <- breaks

    result$MDF.GLS.1 <- list(
        stat = MDF.GLS1,
        cv = cv.MDF.GLS1
    )
    result$MDF.GLS.2 <- list(
        stat = MDF.GLS2,
        cv = cv.MDF.GLS2
    )
    if (breaks == 3) {
        result$MDF.GLS.3 <- list(
            stat = MDF.GLS3,
            cv = cv.MDF.GLS3
        )
    }

    result$MDF.OLS.1 <- list(
        stat = MDF.OLS1,
        cv = cv.MDF.OLS1
    )
    result$MDF.OLS.2 <- list(
        stat = MDF.OLS2,
        cv = cv.MDF.OLS2
    )
    if (breaks == 3) {
        result$MDF.OLS.3 <- list(
            stat = MDF.OLS3,
            cv = cv.MDF.OLS3
        )
    }

    result$UR1 <- UR1
    result$UR <- UR

    class(result) <- "robustURN"
}
