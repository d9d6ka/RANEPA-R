#' @title
#' Generalized supremum ADF test.
#'
#' @param y The input time series of interest.
#' @param trim Trimming parameter to determine the lower and upper bounds.
#' @param const Whether the constant needs to be included.
#' @param add.p.value Whether the p-value is to be returned. This argument is
#' needed to suppress the calculation of p-values during the precalculation of
#' tables needed for the p-values estimating.
#'
#' @references
#' Kurozumi, Eiji, Anton Skrobotov, and Alexey Tsarev.
#' “Time-Transformed Test for the Explosive Bubbles under
#' Non-Stationary Volatility.”
#' arXiv, November 15, 2021. http://arxiv.org/abs/2012.13937.
#'
#' @export
GSADF.test <- function(y,
                       trim = 0.01 + 1.8 / sqrt(length(y)),
                       const = TRUE,
                       add.p.value = TRUE) {
    N <- length(y)

    if (const == FALSE) {
        y <- y - y[1]
    }

    t.values <- c()
    m <- 1
    for (i in 1:(N - floor(trim * N) + 1)) {
        for (j in (i + floor(trim * N) - 1):N) {
            model <- ADF.test(y[i:j], const = const)
            t.values[m] <- model$t.alpha
            m <- m + 1
        }
    }

    ## Take the maximum of the calculated t-statistics.
    GSADF.value <- max(t.values)

    if (add.p.value) {
        if (const == TRUE) {
            cr.values <- .cval_GSADF_with_const
        } else {
            cr.values <- .cval_GSADF_without_const
        }

        p.value <- p.values.SADF(GSADF.value, N, cr.values)
    }

    result <- c(
        list(
            y = y,
            trim = trim,
            const = const,
            t.values = t.values,
            GSADF.value = GSADF.value
        ),
        if (add.p.value) {
            list(p.value = p.value)
        } else {
            NULL
        }
    )

    class(result) <- "sadf"

    return(result)
}
