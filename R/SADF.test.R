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
        t.values[m] <- model$t.beta
        m <- m + 1
    }

    sadf.value <- max(t.values)

    if (const == TRUE) {
        cr.value <- 2.2 # modify
    } else {
        cr.value <- 3.36 # modify
    }
    p.value <- NULL

    if (add.p.value) {
        if (const == TRUE) {
            cr.values <- .cval_SADF_with_const
        } else {
            cr.values <- .cval_SADF_without_const
        }
        p.value <- round(sum(cr.values > sadf.value) / length(cr.values), 4)
    }

    is.explosive <- ifelse(sadf.value > cr.value, 1, 0)

    return(
        list(
            y = y,
            r0 = r0,
            const = const,
            t.values = t.values,
            sadf.value = sadf.value,
            p.value = p.value,
            is.explosive = is.explosive
        )
    )
}
