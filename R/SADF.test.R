#' @title
#' SADF test.
#'
#' @export
SADF.test <- function(y,
                      r0 = 0.01 + 1.8 / sqrt(length(y)),
                      const = TRUE,
                      add.p.value = TRUE) {
    N <- length(y)

    if (const == FALSE) {
        y <- y - y[1]
    }

    t.values <- c()
    m <- 1
    for (j in (floor(r0 * N)):N) {
        model <- ADF.test(y[1:j], const = const)
        t.values[m] <- model$t.alpha
        m <- m + 1
    }

    ## Take the maximum of the calculated t-statistics.
    SADF.value <- max(t.values)

    if (add.p.value) {
        if (const == TRUE) {
            cr.values <- .cval_SADF_with_const
        } else {
            cr.values <- .cval_SADF_without_const
        }

        p.value <- p.values.SADF(SADF.value, N, cr.values)
    }

    result <- c(
        list(
            y = y,
            r0 = r0,
            const = const,
            t.values = t.values,
            SADF.value = SADF.value
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
