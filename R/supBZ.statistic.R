#' @title
#' supBZ.statistic - Calculate supBZ statistic (HLZ, 2018).
#'
#' @export
supBZ.statistic <- function(y,
                            r0 = 0.01 + 1.8 / sqrt(length(y)),
                            sigma.sq = "auto",
                            generalized = FALSE) {
    N <- length(y)

    if (sigma.sq == "auto") {
        # NW estimation.
        my <- (diff(y))^2
        mx <- rep(1, N - 1)
        nw.loocv.model <- NW.loocv(my, mx, kernel = "gauss")
        h.est <- nw.loocv.model$h
        nw.model <- NW.estimation(my, mx,
            kernel = "gauss",
            h = NW.loocv.model$h
        )
        u.hat <- nw.model$u.hat

        # Calculate sigma.sq.
        sigma.sq <- my - nw.model$u.hat
    }

    y <- y - y[1]
    d.y <- diff(y)
    l.y <- y[1:(N - 1)]

    BZ.values <- c()
    m <- 1

    if (!generalized) {
        for (j in (floor(r0 * N)):N) {
            BZ.values[m] <-
                sum(d.y[1:(j - 1)] * l.y[1:(j - 1)] / sigma.sq[1:(j - 1)]) /
                    (sum(l.y[1:(j - 1)]^2 / sigma.sq[1:(j - 1)]))^0.5
            m <- m + 1
        }
    } else {
        for (i in 1:(N - floor(r0 * N) + 1)) {
            for (j in (i + floor(r0 * N) - 1):N) {
                BZ.values[m] <-
                    sum(d.y[i:(j - 1)] * l.y[i:(j - 1)] / sigma.sq[i:(j - 1)]) /
                        (sum(l.y[i:(j - 1)]^2 / sigma.sq[i:(j - 1)]))^0.5
                m <- m + 1
            }
        }
    }

    supBZ.value <- max(BZ.values)

    return(
        list(
            y = y,
            r0 = r0,
            h.est = h.est,
            u.hat = u.hat,
            sigma.sq = sigma.sq,
            BZ.values = BZ.values,
            supBZ.value = supBZ.value
        )
    )
}
