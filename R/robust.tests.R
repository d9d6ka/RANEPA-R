#' @importFrom zeallot %<-%
#'
#' @export
robust.tests.single <- function(y,
                                const = FALSE, trend = FALSE, season = FALSE,
                                trim = 0.15) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)

    x.const <- rep(1, N)
    x.trend <- 1:N

    if (season) {
        SEAS <- cbind(
            x.const,
            seasonal.dummies(N)
        )
        c(., y, ., .) %<-% OLS(y, SEAS)
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


#' @importFrom zeallot %<-%
#'
#' @export
robust.tests.multiple <- function(y,
                                  const = FALSE, season = FALSE,
                                  breaks = 2,
                                  trim = 0.15) {
    if (!is.matrix(y)) y <- as.matrix(y)

     ## Start ##
    N <- nrow(y)

    x.const <- rep(1, N)
    x.trend <- 1:N

    if (season) {
        SEAS <- cbind(
            x.const,
            seasonal.dummies(N)
        )
        c(., y, ., .) %<-% OLS(y, SEAS)
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
