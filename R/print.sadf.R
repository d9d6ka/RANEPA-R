#' @export
print.sadf <- function(object) {
    test.name <- NULL

    if ("SADF.value" %in% names(object)) {
        test.statistic <- object$SADF.value
        test.name <- "SADF"
    } else if ("stadf.value" %in% names(object)) {
        test.statistic <- object$stadf.value
        test.name <- "STADF"
    } else if ("GSADF.value" %in% names(object)) {
        test.statistic <- object$GSADF.value
        test.name <- "GSADF"
    } else if ("gstadf.value" %in% names(object)) {
        test.statistic <- object$gstadf.value
        test.name <- "GSTADF"
    }

    if (!is.null(test.name)) {
        cat(test.name, "test statistic:", test.statistic, "\n")
    }

    if ("supBZ.value" %in% names(object)) {
        cat("supBZ statistic:", object$supBZ.value, "\n")
    }

    if ("p.value" %in% names(object)) {
        cat("P-value:", object$p.value, "\n")
    }

    if ("is.explosive" %in% names(object)) {
        cat(
            "Current process is",
            ifelse(object$is.explosive == 1, "", "not"),
            "explosive\n"
        )
    }
}
