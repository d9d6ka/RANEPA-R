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

    ## Take the maximum of the calculated t-statistics.
    SADF.value <- max(t.values)

    if (add.p.value) {
        if (const == TRUE) {
            cr.values <- .cval_SADF_with_const
        } else {
            cr.values <- .cval_SADF_without_const
        }

        N.list <- as.numeric(names(cr.values))

        if (N < min(N.list)) {
            warning("Too little number of observations, using data for T = 30")
            i.0 <- min(N.list)
        } else {
            i.0 <- max(N.list[N.list <= N])
        }
        p.0 <- sum(cr.values[[as.character(i.0)]] > SADF.value) /
            length(cr.values[[as.character(i.0)]])

        if (N > max(N.list)) {
            i.1 <- max(N.list)
        } else {
            i.1 <- min(N.list[N.list >= N])
        }
        p.1 <- sum(cr.values[[as.character(i.1)]] > SADF.value) /
            length(cr.values[[as.character(i.1)]])

        if (i.0 != i.1) {
            p.value <- p.0 + (p.1 - p.0) * (N - i.0) / (i.1 - i.0)
        } else {
            p.value <- p.0
        }
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
