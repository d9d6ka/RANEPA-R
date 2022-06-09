#' @title
#' Procedure to minimize the GLS-SSR for 1 break point
#'
#' @details
#' See Skrobotov (2018) for further details.
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
#' @importFrom zeallot %<-%
segments.GLS <- function(y,
                         const = FALSE, trend = FALSE,
                         breaks = 1,
                         first.break = NULL, last.break = NULL,
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

    N <- nrow(y)
    x.const <- rep(1, N)
    x.trend <- 1:N

    if (is.null(first.break) || is.null(last.break)) {
        first.break <- floor(trim * N)
        last.break <- floor((1 - trim) * N)
    }

    steps <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.975, 1)

    res.SSR <- Inf
    res.tb <- rep(0, breaks)

    for (alpha in steps) {
        if (const && !trend) {
            for (tb1 in first.break:last.break) {
                DU1 <- as.numeric(x.trend > tb1)
                DT1 <- DU1 * (x.trend - tb1)

                x <- cbind(
                    x.const,
                    x.trend,
                    if (const) DU1 else NULL,
                    if (trend) DT1 else NULL
                )

                c.bar <- N * (alpha - 1)
                c(., resid, ., .) %<-% GLS(y, x, c.bar)

                tmp.SSR <- drop(t(resid) %*% resid)

                if (tmp.SSR < res.SSR) {
                    res.SSR <- tmp.SSR
                    res.tb <- c(tb1)
                }
            }
        } else if (!const && trend) {
            for (tb1 in first.break:(last.break - first.break)) {
                for (tb2 in (tb1 + first.break):last.break) {
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

                    c.bar <- N * (alpha - 1)
                    c(., resid, ., .) %<-% GLS(y, x, c.bar)

                    tmp.SSR <- drop(t(resid) %*% resid)

                    if (tmp.SSR < res.SSR) {
                        res.SSR <- tmp.SSR
                        res.tb <- c(tb1, tb2)
                    }
                }
            }
        } else if (const && trend) {
            for (tb1 in first.break:(last.break - 2 * first.break)) {
                for (tb2 in (tb1 + first.break):(last.break - first.break)) {
                    for (tb3 in (tb2 + first.break):last.break) {
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

                        c.bar <- N * (alpha - 1)
                        c(., resid, ., .) %<-% GLS(y, x, c.bar)

                        tmp.SSR <- drop(t(resid) %*% resid)

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


#' @title
#' Procedure to minimize the SSR for 1 break point
#'
#' @details
#' See Carrion-i-Silvestre and Sansó (2006) for further details.
#'
#' @param beg Sample begin.
#' @param end Sample end.
#' @param first.break First possible break point.
#' @param last.break Last possible break point.
#' @param len Total number of observations.
#' @param SSR.data The matrix of recursive SSR values.
#'
#' @return List containing
#' \describe{
#' \item{SSR}{Optimal SSR value.}
#' \item{break_point}{The point of possible break.}
#' }
segments.OLS.single <- function(beg, end,
                                first.break, last.break,
                                len, SSR.data) {
    tmp_result <- matrix(data = Inf, nrow = len, ncol = 1)

    for (p_break in first.break:last.break) {
        tmp_result[p_break] <- SSR.data[beg, p_break] +
            SSR.data[p_break + 1, end]
    }

    min.SSR <- min(tmp_result[first.break:last.break])
    break.point <- (first.break - 1) +
        which.min(tmp_result[first.break:last.break])

    return(
        list(
            SSR = min.SSR,
            break.point = break.point
        )
    )
}


#' @title
#' Procedure to minimize the SSR for 2 break points
#'
#' @details
#' See Carrion-i-Silvestre and Sansó (2007) for further details.
#'
#' @param y (Tx1)-vector of time series
#' @param model \describe{
#' \item{1}{for the AA (without trend) model.}
#' \item{2}{for the AA (with trend) model.}
#' \item{3}{for the BB model.}
#' \item{4}{for the CC model.}
#' \item{5}{for the AC-CA model.}
#' }
#'
#' @return List containing
#' \describe{
#' \item{resid}{(Tx1) vector of estimated OLS residuals.}
#' \item{tb1}{The first break point.}
#' \item{tb2}{The second break point.}
#' }
#'
#' @importFrom zeallot %<-%
segments.OLS.double <- function(y, model) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)

    r.min <- 0
    SSR.min <- Inf
    tb1.min <- 0
    tb2.min <- 0

    if (1 <= model & model <= 4) {
        for (i in 2:(N - 4)) {
            for (j in (i + 2):(N - 2)) {
                z <- determinants.KPSS.2.breaks(model, N, c(i, j))
                c(., resid, p, .) %<-% OLS(y, z)
                SSR <- drop(t(resid) %*% resid)
                if (SSR < SSR.min) {
                    r.min <- resid
                    SSR.min <- SSR
                    tb1.min <- i
                    tb2.min <- j
                }
            }
        }
    } else if (5 <= model & model <= 7) {
        for (i in 2:(N - 4)) {
            for (j in (i + 2):(N - 2)) {
                z <- determinants.KPSS.2.breaks(model, N, c(i, j))
                c(., resid, p, .) %<-% OLS(y, z)
                SSR <- drop(t(resid) %*% resid)
                if (SSR < SSR.min) {
                    r.min <- resid
                    SSR.min <- SSR
                    tb1.min <- i
                    tb2.min <- j
                }
            }
        }
        for (j in 2:(N - 4)) {
            for (i in (j + 2):(N - 2)) {
                z <- determinants.KPSS.2.breaks(model, N, c(i, j))
                c(., resid, p, .) %<-% OLS(y, z)
                SSR <- drop(t(resid) %*% resid)
                if (SSR < SSR.min) {
                    r.min <- resid
                    SSR.min <- SSR
                    tb1.min <- j
                    tb2.min <- i
                }
            }
        }
    }

    return(
        list(
            resid = r.min,
            tb1 = tb1.min,
            tb2 = tb2.min
        )
    )
}


#' @title
#' Find \eqn{m + 1} optimal partitions
#'
#' @details
#' Based on Bai and Perron (2003).
#'
#' @param y (Tx1)-vector of the dependent variable.
#' @param x (Txk)-vector of the explanatory stochastic regressors.
#' @param m Number of breaks.
#' @param width Minimum spacing between the breaks.
#' @param SSR.data Optional matrix of recursive SSR's.
#'
#' @return List of 2 elements: optimal SSR and the vector of break points.
#'
#' @export
segments.OLS <- function(y, x, m = 1, width = 2, SSR.data = NULL) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)

    if (is.null(SSR.data)) {
        SSR.data <- SSR.matrix(y, x, width)
    }

    if (m == 1) {
        tmp.result <- segments.OLS.single(
            1, N,
            width, N - width,
            N, SSR.data
        )
        optimal.SSR <- tmp.result$SSR
        optimal.break <- tmp.result$break.point
    } else {
        variants <- N - (m + 1) * width + 1
        step.SSR <- matrix(
            data = Inf,
            nrow = variants,
            ncol = 1
        )
        step.break <- matrix(
            data = 0,
            nrow = variants,
            ncol = m
        )
        for (step in 1:m) {
            temp.SSR <- matrix(
                data = Inf,
                nrow = variants,
                ncol = 1
            )
            if (step == 1) {
                for (v in 1:variants) {
                    step.end <- 2 * width + v - 1
                    tmp_res <- segments.OLS.single(
                        1,
                        step.end,
                        width,
                        step.end - width,
                        step.end, SSR.data
                    )
                    step.SSR[v, 1] <- tmp_res$SSR
                    step.break[v, 1] <- tmp_res$break.point
                }
            } else if (step == m) {
                for (v in 1:variants) {
                    temp.SSR[v, 1] <- step.SSR[v, 1] +
                        SSR.data[step * width + v, N]
                }
                optimal.SSR <- min(temp.SSR)
                optimal.index <- which.min(temp.SSR)
                optimal.break <- step.break[optimal.index, ]
                optimal.break[m] <- step * width + optimal.index - 1
            } else {
                next.SSR <- matrix(
                    data = Inf,
                    nrow = variants,
                    ncol = 1
                )
                next.break <- matrix(
                    data = 0,
                    nrow = variants,
                    ncol = m
                )
                for (step.end in ((step + 1) * width):(N - (m - step) * width)) {
                    next.v <- step.end - (step + 1) * width + 1
                    for (v in 1:variants) {
                        temp.SSR[v, 1] <- step.SSR[v, 1] +
                            SSR.data[step * width + v, step.end]
                    }
                    next.SSR[next.v, 1] <- min(temp.SSR)
                    next.index <- which.min(temp.SSR)
                    next.break[next.v, 1:m] <- step.break[next.index, ]
                    next.break[next.v, step] <- step * width + next.index - 1
                }
                step.SSR <- next.SSR
                step.break <- next.break
            }
        }
    }

    return(
        list(
            SSR = optimal.SSR,
            break.point = optimal.break
        )
    )
}
