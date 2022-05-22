#' @export
print.sadf <- function(object) {
    test.name <- NULL

    if ("sadf.value" %in% names(object)) {
        test.statistic <- object$sadf.value
        test.name <- "SADF"
    }
    else if ("stadf.value" %in% names(object)) {
        test.statistic <- object$stadf.value
        test.name <- "STADF"
    }
    else if ("gsadf.value" %in% names(object)) {
        test.statistic <- object$gsadf.value
        test.name <- "GSADF"
    }
    else if ("gstadf.value" %in% names(object)) {
        test.statistic <- object$gstadf.value
        test.name <- "GSTADF"
    }

    if (!is.null(test.name))
        cat(test.name, "test statistic:", test.statistic, "\n")

    if ("supBZ.value" %in% names(object))
        cat("supBZ statistic:", object$supBZ.value, "\n")

    if ("p.value" %in% names(object))
        cat("P-value:", object$p.value, "\n")
}