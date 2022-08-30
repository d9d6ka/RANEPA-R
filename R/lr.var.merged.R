.lr.var.kernel <- function(kernel, alpha, N) {
    if (kernel == "truncated") {
        f.limit <- function(y, l) {
            return(0.6611 * (N * alpha(y, l)$q2)^(1 / 5))
        }
        f.weight <- function(i, l) {
            if (abs(i / l) <= 1) {
                return(1)
            }
            return(0)
        }
    } else if (kernel == "bartlett") {
        f.limit <- function(y, l) {
            return(1.1447 * (N * alpha(y, l)$q1)^(1 / 3))
        }
        f.weight <- function(i, l) {
            x <- i / l
            if (abs(x) <= 1) {
                return(1 - abs(x))
            }
            return(0)
        }
    } else if (kernel == "parzen") {
        f.limit <- function(y, l) {
            return(2.6614 * (N * alpha(y, l)$q2)^(1 / 5))
        }
        f.weight <- function(i, l) {
            x <- i / l
            if (abs(x) <= 0.5) {
                return(1 - 6 * x^2 + 6 * abs(x)^3)
            } else if (abs(x) <= 1) {
                return(2 * (1 - abs(x))^3)
            }
            return(0)
        }
    } else if (kernel == "tukey-hanning") {
        f.limit <- function(y, l) {
            return(1.7462 * (N * alpha(y, l)$q2)^(1 / 5))
        }
        f.weight <- function(i, l) {
            x <- i / l
            if (abs(x) <= 1) {
                return((1 + cos(pi * x)) / 2)
            }
            return(0)
        }
    } else if (kernel == "quadratic") {
        f.limit <- function(y, l) {
            return(1.3221 * (N * alpha(y, l)$q2)^(1 / 5))
        }
        f.weight <- function(i, l) {
            delta <- 6 * pi * (i / l) / 5
            return((3 / delta^2) * (sin(delta) / delta - cos(delta)))
        }
    }

    return(
        list(
            limit = f.limit,
            weight = f.weight
        )
    )
}

.lr.var.alpha.single <- function(y, upper.rho.limit) {
    N <- nrow(y)

    if (!is.null(N) && N > 1) {
        r <- drop(
            (t(y[1:(N - 1), 1]) %*% y[2:N, 1]) /
            (t(y[1:(N - 1), 1]) %*% y[1:(N - 1), 1])
        )

        if (r > upper.rho.limit) {
            r <- upper.rho.limit
        } else if (r < -upper.rho.limit) {
            r <- -upper.rho.limit
        }
    } else {
        r <- y
    }

    return(
        list(
            q1 = 4 * r^2 / (1 + r)^2 / (1 - r)^2,
            q2 = 4 * r^2 / (1 - r)^4
        )
    )
}

.lr.var.alpha.multi <- function(y, upper.rho.limit) {
    N <- nrow(y)
    P <- ncol(y)

    nominator_1 <- 0
    nominator_2 <- 0
    denominator <- 0

    for (i in 1:P) {
        r <- (t(y[1:(N - 1), i]) %*% y[2:N, i]) /
            (t(y[1:(N - 1), i]) %*% y[1:(N - 1), i])
        r <- drop(r)

        if (r > upper.rho.limit) {
            r <- upper.rho.limit
        } else if (r < -upper.rho.limit) {
            r <- -upper.rho.limit
        }

        resids <- y[2:N, i] - y[1:(N - 1), i] * r
        s2 <- mean(resids^2)

        nominator_1 <- nominator_1 + 4 * r^2 * s2^2 / (1 - r)^6 / (1 + r)^2
        nominator_2 <- nominator_2 + 4 * r^2 * s2^2 / (1 - r)^8
        denominator <- denominator + s2^2 / (1 - r)^4
    }

    return(
        list(
            q1 = nominator_1 / denominator,
            q2 = nominator_2 / denominator
        )
    )
}

#' TEST
#' @export
lr.var <- function(y,
                   demean = TRUE,
                   kernel = "bartlett",
                   limit.lags = FALSE,
                   limit.selector = "kpss-q",
                   upper.rho.limit = 0.97,
                   upper.lag.limit = 0.8,
                   recolor = FALSE,
                   max.lag = 0,
                   criterion = "bic") {
    if (!is.matrix(y)) y <- as.matrix(y)

    P <- ncol(y)

    if (!kernel %in% c("truncated",
                       "bartlett",
                       "parzen",
                       "tukey-hanning",
                       "quadratic")) {
        stop("ERROR! Unknown kernel")
    }
    if (!limit.selector %in% c("kpss-q", "kpss-m", "Andrews", "Kurozumi")) {
        stop("ERROR! Unknown criterion")
    }
    if (limit.selector == "Kurozumi" && is.null(upper.lag.limit)) {
        stop("ERROR! Upper limit is needed for Kurozumi proposal")
    }
    if (limit.selector == "Kurozumi" && P > 1) {
        stop("ERROR! Kurozumi proposal is for a single variable case")
    }
    if (recolor && P > 1) {
        stop("ERROR! Recolorization is for a single variable case")
    }
    if (!criterion %in% c("bic", "aic", "lwz")) {
        stop("ERROR! Unknown criterion")
    }
    if (recolor) {
        limit.lags <- TRUE
        limit.selector <- "Andrews"
        upper.lag.limit <- NULL
    }

    if (demean) {
        for (i in 1:P) {
            mean.y <- mean(y[, i])
            y[, i] <- y[, i] - mean.y
        }
    }

    N <- nrow(y)

    if (P == 1) {
        funcs <- .lr.var.kernel(kernel, .lr.var.alpha.single, N)
    } else {
        funcs <- .lr.var.kernel(kernel, .lr.var.alpha.multi, N)
    }

    if (recolor) {
        info.crit.min <- log(drop(t(y) %*% y) / (N - max.lag))

        tmp.AR <- AR(y, NULL, max.lag, criterion)

        if (info.crit.min < tmp.AR$criterion) {
            rho <- 0
            k <- 0
        } else {
            rho <- tmp.AR$beta
            k <- tmp.AR$lag
            y <- tmp.AR$residuals
        }
    } else {
        k <- 1
    }

    if (!limit.lags) {
        limit <- N - 1
    } else if (limit.selector == "kpss-q") {
        limit <- 4 * ((N / 100)^(1 / 4))
    } else if (limit.selector == "kpss-m") {
        limit <- 12 * ((N / 100)^(1 / 4))
    } else {
        if (k > 0) {
            limit <- funcs$limit(y, upper.rho.limit)

            if (limit.selector == "Kurozumi") {
                upper.lag.limit <- funcs$limit(upper.lag.limit, upper.rho.limit)
                limit <- min(limit, upper.lag.limit)
            }
        } else {
            limit <- 0
        }
    }
    limit <- trunc(limit)

    lrv <- (t(y) %*% y) / N
    if (k > 0) {
        for (i in 1:limit) {
            if (i < N - 1) {
                lrv <- lrv + funcs$weight(i, limit) * (
                    (t(y[1:(N - i), ]) %*% y[(1 + i):N, ]) / N +
                    t(t(y[1:(N - i), ]) %*% y[(1 + i):N, ]) / N
                )
            } else {
                lrv <- lrv + funcs$weight(i, limit) * as.vector(
                    (t(y[1:(N - i), ]) %*% y[(1 + i):N, ]) / N +
                    t(t(y[1:(N - i), ]) %*% y[(1 + i):N, ]) / N
                )
            }
        }
    }

    if (recolor) {
        lrv.recolored <- lrv / (1 - sum(rho))^2
        lrv <- min(lrv.recolored, N * 0.15 * lrv)
    }

    return(drop(lrv))
}
