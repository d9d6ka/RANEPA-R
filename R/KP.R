#' @title
#' Kejrival-Perron procedure of breaks number detection.
#'
#' @param y The input series of interest.
#' @param const Allowing the break in constant.
#' @param breaks Number of breaks.
#' @param criterion Needed information criterion: aic, bic, hq or lwz.
#' @param trim A trimming value for a possible break date bounds.
#'
#' @references
#' Kejriwal, Mohitosh, and Pierre Perron.
#' “A Sequential Procedure to Determine the Number of Breaks in Trend
#' with an Integrated or Stationary Noise Component:
#' Determination of Number of Breaks in Trend.”
#' Journal of Time Series Analysis 31, no. 5 (September 2010): 305–28.
#' https://doi.org/10.1111/j.1467-9892.2010.00666.x.
#'
#' @importFrom zeallot %<-%
#'
#' @export
KP <- function(y,
               const = FALSE,
               breaks = 1,
               criterion = "aic",
               trim = 0.15) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)
    k.max <- trunc(12 * (N / 100)^(1 / 4))

    model <- as.numeric(const) + 1
    trim.pos <- which(c(0.01, 0.05, 0.1, 0.15, 0.25) == trim)

    res <- 0

    for (l in 0:(breaks - 1)) {
        test.stat <- PY.sequential(y, const, l, criterion, trim, k.max)
        c.v <- .cval_KP[[model]][[trim.pos]]

        if (test.stat < c.v[2, l + 1]) {
            res <- l
            break
        }
    }

    return(res)
}
