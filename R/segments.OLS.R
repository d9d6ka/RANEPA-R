#' @title
#' Find \eqn{m + 1} optimal partitions
#'
#' @param y (Tx1)-vector of the dependent variable.
#' @param x (Txk)-vector of the explanatory stochastic regressors.
#' @param m Number of breaks.
#' @param width Minimum spacing between the breaks.
#' @param SSR.data Optional matrix of recursive SSR's.
#'
#' @return List of 2 elements: optimal SSR and the vector of break points.
#'
#' @references
#' Bai, Jushan, and Pierre Perron.
#' “Computation and Analysis of Multiple Structural Change Models.”
#' Journal of Applied Econometrics 18, no. 1 (2003): 1–22.
#' https://doi.org/10.1002/jae.659.
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
        res.ssr <- tmp.result$SSR
        res.break <- tmp.result$break.point
    } else {
        variants <- N - (m + 1) * width + 1
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
                    tmp.res <- segments.OLS.single(
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
                        SSR.data[step * width + v, N]
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
                for (step.end in ((step + 1) * width):(N - (m - step) * width)) {
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
