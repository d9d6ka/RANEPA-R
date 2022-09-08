#' @title
#' Procedure to minimize the GLS-SSR for 1 break point
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
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

                c_bar <- n.obs * (alpha - 1)
                resids <- GLS(y, x, c_bar)$residuals

                tmp.SSR <- drop(t(resids) %*% resids)

                if (tmp.SSR < res.SSR) {
                    res.SSR <- tmp.SSR
                    res.tb <- c(tb1)
                }
            }
        } else if (!const && trend) {
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
        } else if (const && trend) {
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
