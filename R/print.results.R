#' @export
print.sadf <- function(object) {
    test_name <- NULL

    if ("SADF.value" %in% names(object)) {
        test_statistic <- object$SADF.value
        test_name <- "SADF"
    } else if ("stadf.value" %in% names(object)) {
        test_statistic <- object$stadf.value
        test_name <- "STADF"
    } else if ("GSADF.value" %in% names(object)) {
        test_statistic <- object$GSADF.value
        test_name <- "GSADF"
    } else if ("gstadf.value" %in% names(object)) {
        test_statistic <- object$gstadf.value
        test_name <- "GSTADF"
    }

    if (!is.null(test_name)) {
        cat(test_name, "test statistic:", test_statistic, "\n")
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

#' @importFrom stringr str_split
#' @export
print.mdfHLT <- function(object) {
    if (object$const && !object$trend) {
        cat("Model 0: Structural change in intercept\n")
        cat("Y{t}=a0+a1*DU+b0*t+e{t}\n",
            "    where DU=1(t>TB)\n")
    } else if (!object$const && object$trend) {
        cat("Model 1: Structural change in slope\n")
        cat("Y{t}=a0+b0*t+b1*DT+e{t}\n",
            "    where DT=1(t>TB)*(t-TB)\n")
    } else if (object$const && object$trend) {
        cat("Model 2: Structural change in both intercept and slope\n")
        cat("Y{t}=a0+a1*DU+b0*t+b1*DT+e{t}\n",
            "    where DU=1(t>TB) and DT=1(t>TB)*(t-TB)\n")
    }
    cat("\n")

    cat("The break date is estimated in", object$break.time, "\n")
    cat("\n")

    cat("Robust tests for the break with uncertaint over errors",
        "(integrated or stationary):\n\n",
        "\tstat\tc.v.\n")
    for (v in c("HLT", "PY")) {
        cat(sprintf("%-7s\t%.4f\t%.4f\n", v, object[[v]]$stat, object[[v]]$cv))
        if (object[[v]]$stat > object[[v]]$cv)
            cat("reject\n")
        else
            cat("fails to reject\n")
        cat("\n")
    }
    cat("\n")

    cat("Unit root tests:\n\n",
        "\tstat\tc.v.\n")
    for (v in c("DF.GLS", "DF.OLS", "MDF.GLS", "MDF.OLS", "MDF.t")) {
        cat(sprintf("%-7s\t%.4f\t%.4f\n", v, object[[v]]$stat, object[[v]]$cv))
        if (object[[v]]$stat < object[[v]]$cv)
            cat("reject\n")
        else
            cat("fails to reject\n")
        cat("\n")
    }
    cat("\n")

    cat("Testing strategies:\n\n")
    for (v in c("A.HLT", "A.PY", "UR.HLT", "UR.PY")) {
        tmp_str <- stringr::str_split(v, "\\.")
        cat(tmp_str[[1]][1], "*(t_", tmp_str[[1]][2], ", s_alpha): ", sep = "")
        if (object[[v]] == 1)
            cat("reject\n")
        else
            cat("fails to reject\n")
        cat("\n")
    }
}


#' @importFrom stringr str_split
#' @export
print.mdfHLTN <- function(object) {
    cat("\t\tstat\tc.v.\n\n")

    for (v in c(
        "MDF.GLS.1",
        "MDF.GLS.2",
        if (object$breaks == 3) "MDF.GLS.3" else NULL,
        "MDF.OLS.1",
        "MDF.OLS.2",
        if (object$breaks == 3) "MDF.OLS.3" else NULL
    )) {
        cat(sprintf("%-9s:\t%.4f\t%.4f\n", v, object[[v]]$stat, object[[v]]$cv))
        if (object[[v]]$stat < object[[v]]$cv)
            cat("reject\n")
        else
            cat("fails to reject\n")
        cat("\n")
    }
    cat("\n")

    cat(sprintf("UR^%d(s.alpha): ", object$breaks))
    if (object$UR1 == 1)
        cat("reject\n")
    else
        cat("fails to reject\n")

    cat(sprintf("UR^%d(s.alpha, %d): ", object$breaks, object$breaks.star))
    if (object$UR == 1)
        cat("reject\n")
    else
        cat("fails to reject\n")
    cat("\n")
}


#' @importFrom stringr str_split
#' @export
print.mdfHHLT <- function(object) {
    cat("\t\tstat\tc.v.\t wild c.v.\n\n")
    for (v in c("MZa", "MSB", "MZt", "ADF")) {
        cat(
            sprintf(
                "%s stat:\t%.4f\t%.4f\t%.4f\n",
                v,
                object[[v]]$stat,
                object[[v]]$cv,
                object[[v]]$cv.bootstrap
            )
        )
    }
}
