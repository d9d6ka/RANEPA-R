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

    SADF.value <- max(t.values)

    if (const == TRUE) {
        cr.value <- 2.2 # modify
    } else {
        cr.value <- 3.36 # modify
    }

    if (add.p.value) {
        if (const == TRUE) {
            cr.values <- .cval_SADF_with_const
        } else {
            cr.values <- .cval_SADF_without_const
        }

        log.N = log(N / 100, base = 2)

        i.0 <- max(floor(log.N), 0)
        p.0 <- sum(cr.values[[i.0 + 1]] > SADF.value) /
            length(cr.values[[i.0 + 1]])

        i.1 <- min(ceiling(log.N), 4)
        p.1 <- sum(cr.values[[i.1 + 1]] > SADF.value) /
            length(cr.values[[i.1 + 1]])

        p.value <- p.0 + (p.1 - p.0) * (N - 2^i.0 * 100) /
            ((2^i.1 - 2^i.0) * 100)
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
