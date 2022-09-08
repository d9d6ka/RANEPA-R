#' @title
#' A wrapping function around [MDF.single].
#'
#' @param y A series of interest.
#' @param const Whether the constant term should be included.
#' @param trend Whether the trend term should be included.
#' @param season Whether the seasonal adjustment is needed.
#' @param trim Trimming value for a possible break date bounds.
#'
#' @export
robust.tests.single <- function(y,
                                const = FALSE, trend = FALSE, season = FALSE,
                                trim = 0.15) {
    if (!is.matrix(y)) y <- as.matrix(y)

    n.obs <- nrow(y)

    x.const <- rep(1, n.obs)

    if (season) {
        SEAS <- cbind(
            x.const,
            seasonal.dummies(n.obs)
        )
        y <- OLS(y, SEAS)$residuals
    }

    result <- MDF.single(
        y = y,
        const = const,
        trend = trend,
        trim = trim
    )

    result <- append(result, list(season = season), 2)
    class(result) <- "robustUR"

    return(result)
}


#' @title
#' A wrapping function around [KP] and [MDF.multiple].
#'
#' @param y A series of interest.
#' @param const Whether the constant term should be included.
#' @param season Whether the seasonal adjustment is needed.
#' @param breaks Number of breaks.
#' @param trim Trimming value for a possible break date bounds.
#'
#' @export
robust.tests.multiple <- function(y,
                                  const = FALSE, season = FALSE,
                                  breaks = 2,
                                  trim = 0.15) {
    if (!is.matrix(y)) y <- as.matrix(y)

     ## Start ##
    n.obs <- nrow(y)

    x.const <- rep(1, n.obs)

    if (season) {
        SEAS <- cbind(
            x.const,
            seasonal.dummies(n.obs)
        )
        y <- OLS(y, SEAS)$residuals
    }

    m.star <- KP(
        y = y,
        const = const, breaks = breaks,
        criterion = "aic", trim = trim
    )

    result <- MDF.multiple(
        y = y,
        const = const,
        breaks = breaks,
        breaks.star = m.star,
        trim = trim,
        ZA = FALSE
    )

    result <- append(result, list(season = season), 1)
    class(result) <- "robustURN"

    return(result)
}
