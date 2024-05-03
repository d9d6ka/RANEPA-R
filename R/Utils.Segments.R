#' @title
#' Procedure to minimize the SSR for 1 break point
#'
#' @param beg Start of the sample.
#' @param end End of the sample.
#' @param first.break First possible break point.
#' @param last.break Last possible break point.
#' @param len Total number of observations.
#' @param SSR.data The matrix of recursive SSR values.
#'
#' @return A list of:
#' * `SSR`: Optimal SSR value,
#' * `break.point`: The point of possible break.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “Testing the Null of Cointegration with Structural Breaks.”
#' Oxford Bulletin of Economics and Statistics 68, no. 5 (October 2006): 623–46.
#' https://doi.org/10.1111/j.1468-0084.2006.00180.x.
#'
#' @keywords internal
segments.OLS.1.break <- function(beg,
                                 end,
                                 first.break,
                                 last.break,
                                 len,
                                 SSR.data) {
    tmp.result <- matrix(data = Inf, nrow = len, ncol = 1)

    for (tb in first.break:last.break) {
        tmp.result[tb] <- SSR.data[beg, tb] +
            SSR.data[tb + 1, end]
    }

    min.ssr <- min(tmp.result[first.break:last.break])
    break.point <- (first.break - 1) +
        which.min(tmp.result[first.break:last.break])

    return(
        list(
            SSR = min.ssr,
            break.point = break.point
        )
    )
}


#' @title
#' Procedure to minimize the SSR for 2 break points
#'
#' @param y A time series of interest.
#' @param model A scalar equal to
#' * 1: for the AA (without trend) model,
#' * 2: for the AA (with trend) model,
#' * 3: for the BB model,
#' * 4: for the CC model,
#' * 5: for the AC-CA model.
#'
#' @return A list of
#' * resid: (Tx1) vector of estimated OLS residuals,
#' * tb1: The first break point,
#' * tb2: The second break point.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “The KPSS Test with Two Structural Breaks.”
#' Spanish Economic Review 9, no. 2 (May 16, 2007): 105–27.
#' https://doi.org/10.1007/s10108-006-9017-8.
#'
#' @keywords internal
segments.OLS.2.breaks <- function(y,
                                  model) {
    if (!is.matrix(y)) y <- as.matrix(y)

    n.obs <- nrow(y)

    res.r <- 0
    cur.ssr <- Inf
    res.tb1 <- 0
    res.tb2 <- 0

    if (1 <= model && model <= 4) {
        for (i in 2:(n.obs - 4)) {
            for (j in (i + 2):(n.obs - 2)) {
                z <- determinants.KPSS.2.breaks(model, n.obs, c(i, j))
                resids <- OLS(y, z)$residuals
                ssr <- drop(t(resids) %*% resids)
                if (ssr < cur.ssr) {
                    res.r <- resids
                    cur.ssr <- ssr
                    res.tb1 <- i
                    res.tb2 <- j
                }
            }
        }
    } else if (5 <= model && model <= 7) {
        for (i in 2:(n.obs - 4)) {
            for (j in (i + 2):(n.obs - 2)) {
                z <- determinants.KPSS.2.breaks(model, n.obs, c(i, j))
                resids <- OLS(y, z)$residuals
                ssr <- drop(t(resids) %*% resids)
                if (ssr < cur.ssr) {
                    res.r <- resids
                    cur.ssr <- ssr
                    res.tb1 <- i
                    res.tb2 <- j
                }
            }
        }
        for (j in 2:(n.obs - 4)) {
            for (i in (j + 2):(n.obs - 2)) {
                z <- determinants.KPSS.2.breaks(model, n.obs, c(i, j))
                resids <- OLS(y, z)$residuals
                ssr <- drop(t(resids) %*% resids)
                if (ssr < cur.ssr) {
                    res.r <- resids
                    cur.ssr <- ssr
                    res.tb1 <- j
                    res.tb2 <- i
                }
            }
        }
    }

    return(
        list(
            residuals = res.r,
            tb1 = res.tb1,
            tb2 = res.tb2
        )
    )
}


#' @title
#' Find \eqn{m + 1} optimal partitions
#'
#' @param y (Tx1)-vector of the dependent variable.
#' @param x (Txk)-vector of the explanatory stochastic regressors.
#' @param m Number of breaks.
#' @param width Minimum spacing between the breaks.
#' @param SSR.data Optional matrix of recursive SSR's.
#'
#' @return A list of:
#' * optimal SSR,
#' * the vector of break points.
#'
#' @references
#' Bai, Jushan, and Pierre Perron.
#' “Computation and Analysis of Multiple Structural Change Models.”
#' Journal of Applied Econometrics 18, no. 1 (2003): 1–22.
#' https://doi.org/10.1002/jae.659.
#'
#' @keywords internal
segments.OLS.N.breaks <- function(y,
                                  x,
                                  m = 1,
                                  width = 2,
                                  SSR.data = NULL) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    n.obs <- nrow(y)

    if (is.null(SSR.data)) {
        SSR.data <- SSR.matrix(y, x, width)
    }

    if (m == 1) {
        tmp.result <- segments.OLS.1.break(
            1, n.obs,
            width, n.obs - width,
            n.obs, SSR.data
        )
        res.ssr <- tmp.result$SSR
        res.break <- tmp.result$break.point
    } else {
        variants <- n.obs - (m + 1) * width + 1
        cur.ssr <- matrix(
            data = Inf,
            nrow = variants,
            ncol = 1
        )
        cur.breaks <- matrix(
            data = 0,
            nrow = variants,
            ncol = m
        )
        for (step in 1:m) {
            tmp.ssr <- matrix(
                data = Inf,
                nrow = variants,
                ncol = 1
            )
            if (step == 1) {
                for (v in 1:variants) {
                    step.end <- 2 * width + v - 1
                    tmp.res <- segments.OLS.1.break(
                        1,
                        step.end,
                        width,
                        step.end - width,
                        step.end, SSR.data
                    )
                    cur.ssr[v, 1] <- tmp.res$SSR
                    cur.breaks[v, 1] <- tmp.res$break.point
                }
            } else if (step == m) {
                for (v in 1:variants) {
                    tmp.ssr[v, 1] <- cur.ssr[v, 1] +
                        SSR.data[step * width + v, n.obs]
                }
                res.ssr <- min(tmp.ssr)
                res.index <- which.min(tmp.ssr)
                res.break <- cur.breaks[res.index, ]
                res.break[m] <- step * width + res.index - 1
            } else {
                new.ssr <- matrix(
                    data = Inf,
                    nrow = variants,
                    ncol = 1
                )
                new.breaks <- matrix(
                    data = 0,
                    nrow = variants,
                    ncol = m
                )
                for (step.end in ((step + 1) * width):(n.obs - (m - step) * width)) { # nolint
                    new.v <- step.end - (step + 1) * width + 1
                    for (v in 1:variants) {
                        tmp.ssr[v, 1] <- cur.ssr[v, 1] +
                            SSR.data[step * width + v, step.end]
                    }
                    new.ssr[new.v, 1] <- min(tmp.ssr)
                    new.index <- which.min(tmp.ssr)
                    new.breaks[new.v, 1:m] <- cur.breaks[new.index, ]
                    new.breaks[new.v, step] <- step * width + new.index - 1
                }
                cur.ssr <- new.ssr
                cur.breaks <- new.breaks
            }
        }
    }

    return(
        list(
            SSR = res.ssr,
            break.point = res.break
        )
    )
}


#' @title
#' Procedure to minimize the GLS-SSR for 1 break point
#'
#' @param y Variable of interest.
#' @param const Whether there is a break in the constant.
#' @param trend Whether there is a break in the trend.
#' @param breaks Number of breaks.
#' @param first.break First possible break point.
#' @param last.break Last possible break point.
#' @param trim Trim value to calculate `first.break` and `last.break`
#' if not provided.
#'
#' @return The point of possible break.
#'
#' @references
#' Skrobotov, Anton.
#' “On Trend Breaks and Initial Condition in Unit Root Testing.”
#' Journal of Time Series Econometrics 10, no. 1 (2018): 1–15.
#' https://doi.org/10.1515/jtse-2016-0014.
#'
#' @keywords internal
segments.GLS <- function(y,
                         const = FALSE,
                         trend = FALSE,
                         breaks = 1,
                         first.break = NULL,
                         last.break = NULL,
                         trim = 0.15) {
    if (!is.matrix(y)) y <- as.matrix(y)

    if (breaks < 1) {
        warning("At least one break is needed!")
        breaks <- 1
    }
    if (breaks > 3) {
        warning("More than three breaks are not supported at the moment!")
        breaks <- 3
    }

    n.obs <- nrow(y)
    x.const <- rep(1, n.obs)
    x.trend <- 1:n.obs

    if (is.null(first.break) || is.null(last.break)) {
        first.break <- floor(trim * n.obs) + 1
        last.break <- floor((1 - trim) * n.obs) + 1
    }
    width <- first.break - 1

    steps <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.975, 1)

    res.SSR <- Inf
    res.tb <- rep(0, breaks)

    for (alpha in steps) {
        if (breaks == 1) {
            for (tb1 in first.break:last.break) {
                DU1 <- as.numeric(x.trend > tb1)
                DT1 <- DU1 * (x.trend - tb1)

                x <- cbind(
                    x.const,
                    x.trend,
                    if (const) DU1 else NULL,
                    if (trend) DT1 else NULL
                )

                c_bar <- n.obs * (alpha - 1)
                resids <- GLS(y, x, c_bar)$residuals

                tmp.SSR <- drop(t(resids) %*% resids)

                if (tmp.SSR < res.SSR) {
                    res.SSR <- tmp.SSR
                    res.tb <- c(tb1)
                }
            }
        } else if (breaks == 2) {
            for (tb1 in first.break:(last.break - width)) {
                for (tb2 in (tb1 + width):last.break) {
                    DU1 <- as.numeric(x.trend > tb1)
                    DT1 <- DU1 * (x.trend - tb1)
                    DU2 <- as.numeric(x.trend > tb2)
                    DT2 <- DU2 * (x.trend - tb2)

                    x <- cbind(
                        x.const,
                        x.trend,
                        if (const) DU1 else NULL,
                        if (trend) DT1 else NULL,
                        if (const) DU2 else NULL,
                        if (trend) DT2 else NULL
                    )

                    c_bar <- n.obs * (alpha - 1)
                    resids <- GLS(y, x, c_bar)$residuals

                    tmp.SSR <- drop(t(resids) %*% resids)

                    if (tmp.SSR < res.SSR) {
                        res.SSR <- tmp.SSR
                        res.tb <- c(tb1, tb2)
                    }
                }
            }
        } else if (breaks == 3) {
            for (tb1 in first.break:(last.break - 2 * width)) {
                for (tb2 in (tb1 + width):(last.break - width)) {
                    for (tb3 in (tb2 + width):last.break) {
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
                            if (trend) DT1 else NULL,
                            if (const) DU2 else NULL,
                            if (trend) DT2 else NULL,
                            if (const) DU3 else NULL,
                            if (trend) DT3 else NULL
                        )

                        c_bar <- n.obs * (alpha - 1)
                        resids <- GLS(y, x, c_bar)$residuals

                        tmp.SSR <- drop(t(resids) %*% resids)

                        if (tmp.SSR < res.SSR) {
                            res.SSR <- tmp.SSR
                            res.tb <- c(tb1, tb2, tb3)
                        }
                    }
                }
            }
        }
    }

    return(res.tb)
}
