#' @title
#' Calculate sup sign-based SADF statistic (HLZ, 2019)).
#'
#' @importFrom zeallot %<-%
#'
#' @export
supSBADF.statistic <- function(y,
                               r0 = 0.01 + 1.8 / sqrt(length(y)),
                               generalized = FALSE) {
    N <- length(y)

    # Calculate C.t.
    C.t <- cumsum(sign(diff(y)))

    SBADF.values <- c()
    m <- 1

    if (!generalized) {
        for (j in (floor(r0 * N)):N) {
            c(., ., ., t.beta) %<-% OLS(diff(C.t)[1:j], C.t[1:j])
            SBADF.values[m] <- drop(t.beta)
            m <- m + 1
        }
    } else {
        for (i in 1:(N - floor(r0 * N) + 1)) {
            for (j in (i + floor(r0 * N) - 1):N) {
                c(., ., ., t.beta) %<-% OLS(diff(C.t)[i:j], C.t[i:j])
                SBADF.values[m] <- drop(t.beta)
                m <- m + 1
            }
        }
    }

    supSBADF.value <- max(SBADF.values)

    return(
        list(
            y = y,
            r0 = r0,
            C.t = C.t,
            SBADF.values = SBADF.values,
            supSBADF.value = supSBADF.value
        )
    )
}
