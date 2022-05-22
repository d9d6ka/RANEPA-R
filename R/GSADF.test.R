#' @title
#' GSADF test.
#'
#' @export
GSADF.test <- function(y,
                       r0 = 0.01 + 1.8 / sqrt(length(y)),
                       const = TRUE,
                       add.p.value = TRUE) {
    N <- length(y)

    if (const == FALSE) {
        y <- y - y[1]
    }

    t.values <- c()
    m <- 1
    for (i in 1:(N - floor(r0 * N) + 1)) {
        for (j in (i + floor(r0 * N) - 1):N) {
            model <- ADF.test(y[i:j], const = const)
            t.values[m] <- model$t.beta
            m <- m + 1
        }
    }

    GSADF.value <- max(t.values)

    if (const == TRUE) {
        cr.value <- 1.524
    } else {
        cr.value <- 2.781
    }

    if (add.p.value) {
        if (const == TRUE) {
            cr.values <- .cval_GSADF_with_const
        } else {
            cr.values <- .cval_GSADF_without_const
        }
        p.value <- round(sum(cr.values > GSADF.value) / length(cr.values), 4)
    }

    is.explosive <- ifelse(GSADF.value > cr.value, 1, 0)

    result <- c(
        list(
            y = y,
            r0 = r0,
            const = const,
            t.values = t.values,
            GSADF.value = GSADF.value,
            cr.value = cr.value,
            is.explosive = is.explosive
        ),
        if (add.p.value) {
            list(p.value = p.value)
        } else NULL
    )

    class(result) <- "sadf"

    return(result)
}
