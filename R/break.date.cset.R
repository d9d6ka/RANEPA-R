#' @title
#' Confidence sets for the break date in dointegrating regressions
#'
#' @param y A LHS variable of interest.
#' @param trend Whether the trend is to be included.
#' @param zb I(1) regressors with break.
#' @param zf I(1) regressors without break.
#' @param z.lead,z.lag Number of leads and lags of `z` regressors.
#' If any is NULL then both are estimated using informational `criterion`.
#' @param conf.level Confidence level to obtain appropriate critical values.
#' @param trim The trimming parameter to find the lower and upper bounds of
#' possible break date.
#' @param criterion A criterion for lead and lag number estimation.
#'
#' @returns
#' A list of confidence sets.
#'
#' @references
#' Kurozumi, Eiji, and Anton Skrobotov.
#' “Confidence Sets for the Break Date in Cointegrating Regressions.”
#' Oxford Bulletin of Economics and Statistics 80, no. 3 (2018): 514–35.
#' https://doi.org/10.1111/obes.12223.
break.date.cset <- function(y,
                            trend = FALSE,
                            zb = NULL,
                            zf = NULL,
                            z.lead = NULL,
                            z.lag = NULL,
                            conf.level = 0.9,
                            trim = 0.05,
                            criterion = "bic") {
    if (!is.matrix(y)) y <- as.matrix(y)

    if (is.null(z.lead) || is.null(z.lag)) {
        ll.est <- select.lead.lag(y, trend, zb, zf, trim, criterion)
        z.lead <- ll.est$lead
        z.lag <- ll.est$lag
    }

    n.obs <- nrow(y)
    td <- as.matrix(1:n.obs)

    wb <- cbind(
        rep(1, n.obs),
        if (trend) (1:n.obs) / n.obs else NULL,
        zb
    )
    z <- cbind(zb, zf)

    p.zb <- if (!is.null(zb)) ncol(zb) else 0
    p.zf <- if (!is.null(zf)) ncol(zf) else 0

    td <- as.matrix(td[2:n.obs, ])
    y <- as.matrix(y[2:n.obs, ])
    wb <- as.matrix(wb[2:n.obs, ])
    d.z <- as.matrix(z[2:n.obs, ] - z[1:(n.obs - 1), ])

    wf <- cbind(
        if (!is.null(zf)) as.matrix(zf[2:n.obs, ]) else NULL,
        d.z
    )

    N1 <- nrow(y)
    td <- as.matrix(td[(z.lag + 1):(N1 - z.lead), ])
    y <- as.matrix(y[(z.lag + 1):(N1 - z.lead), ])
    wb <- as.matrix(wb[(z.lag + 1):(N1 - z.lead), ])
    wf <- as.matrix(wf[(z.lag + 1):(N1 - z.lead), ])

    if (z.lead > 0) {
        for (k in 1:z.lead) {
            wf <- cbind(
                wf,
                d.z[((z.lag + 1) - k):((N1 - z.lead) - k), ]
            )
        }
    }
    if (z.lag > 0) {
        for (k in 1:z.lag) {
            wf <- cbind(
                wf,
                d.z[((z.lag + 1) + k):((N1 - z.lead) + k), ]
            )
        }
    }

    N2 <- nrow(y)
    first.break1 <- trunc(2 * trim * N2)
    last.break1 <- trunc((1 - 2 * trim) * N2)
    first.break2 <- trunc(trim * N2)
    last.break2 <- trunc((1 - trim) * N2)

    cset.sup <- rep(0, N2)
    cset.avg <- rep(0, N2)
    cset.exp <- rep(0, N2)
    cset.bls <- rep(0, N2)

    w <- cbind(wb, wf)
    u.hat <- OLS(y, w)$residuals
    ssr.0 <- c(t(u.hat) %*% u.hat)
    est.date <- N2

    for (tb in first.break1:last.break1) {
        wb1 <- rbind(
            matrix(0, tb, ncol(wb)),
            as.matrix(wb[(tb + 1):N2, ])
        )
        w <- cbind(wb, wb1, wf)
        u.hat <- OLS(y, w)$residuals
        ssr.1 <- c(t(u.hat) %*% u.hat)
        if (ssr.1 < ssr.0) {
            ssr.0 <- ssr.1
            est.date <- tb
        }
    }

    wb1e <- rbind(
        matrix(0, est.date, ncol(wb)),
        as.matrix(wb[(est.date + 1):N2, ])
    )
    td[est.date, 1] <- -1

    w <- cbind(wb1e, wb, wf)
    b.hat <- solve(t(w) %*% w) %*% t(w) %*% y
    u.hat <- y - w %*% b.hat

    lrv.u <- lr.var(
        u.hat,
        demean = FALSE,
        kernel = "quadratic",
        limit.lags = TRUE,
        limit.selector = "Andrews"
    )

    l.hat <- (wb[est.date, ] %*% b.hat[seq_len(ncol(wb))])^2 / lrv.u
    if (conf.level == 0.9) {
        c.bls <- 7.686962
    } else if (conf.level == 0.95) {
        c.bls <- 11.03281
    }

    bdd <- trunc(c.bls / l.hat)
    bls.l <- est.date - bdd - 1
    bls.u <- est.date + bdd + 1
    if (bls.l < 1) {
        bls.l <- 1
    }
    if (bls.u > N2) {
        bls.u <- N2
    }

    cset.bls[bls.l:bls.u] <- 1

    for (tb in first.break1:last.break1) {
        lambda.1 <- tb / N2

        wb1 <- rbind(
            matrix(0, tb, ncol(wb)),
            as.matrix(wb[(tb + 1):N2, ])
        )

        w <- cbind(wb, wb1, wf)

        y.hat <- OLS(y, w)$residuals

        if (abs(tb - est.date) > ncol(wb)) {
            we <- cbind(w, wb1e)
        } else {
            we <- w
        }
        be.hat <- solve(t(we) %*% we) %*% t(we) %*% y
        u.hat <- y - we %*% be.hat
        lrv.u2 <- lr.var(
            u.hat,
            demean = FALSE,
            kernel = "quadratic",
            limit.lags = TRUE,
            limit.selector = "Andrews"
        )

        sup.stat <- 0
        avg.stat <- 0
        exp.stat <- 0

        nbreak <- 0
        dbreak <- 0

        for (tb2 in first.break2:last.break2) {
            lambda.2 <- tb2 / N2

            if (abs(lambda.2 - lambda.1) <= 0.05) {
                dbreak <- dbreak + 1
            } else {
                wb2 <- rbind(
                    matrix(0, tb2, ncol(wb)),
                    as.matrix(wb[(tb2 + 1):N2, ])
                )

                r <- wb2 - wb1
                r.hat <- OLS(r, w)$residuals

                g <- t(r.hat) %*% y.hat
                h <- t(r.hat) %*% r.hat
                ghg <- c(t(g) %*% solve(h) %*% g) / lrv.u2

                if (sup.stat < ghg) {
                    sup.stat <- ghg
                }
                
                avg.stat <- avg.stat + ghg
                exp.stat <- exp.stat + exp(ghg / 2)

                nbreak <- nbreak + 1
            }
        }

        avg.stat <- avg.stat / (nbreak - dbreak)
        exp.stat <- log(exp.stat / (nbreak - dbreak))

        cv <- critical.values.break.date.cset(
            lambda.1,
            trend,
            conf.level,
            p.zb,
            p.zf
        )
        if (sup.stat <= cv$cval_sup) {
            cset.sup[t] <- 1
        }
        if (avg.stat <= cv$cval_avg) {
            cset.avg[t] <- 1
        }
        if (exp.stat <= cv$cval_exp) {
            cset.exp[t] <- 1
        }
    }

    if (z.lead == 0) {
        td <- rbind(
            matrix(0, (z.lag + 1), 1),
            td
        )
        cset.sup <- rbind(
            matrix(0, (z.lag + 1), 1),
            matrix(cset.sup, N2, 1)
        )
        cset.avg <- rbind(
            matrix(0, (z.lag + 1), 1),
            matrix(cset.avg, N2, 1)
        )
        cset.exp <- rbind(
            matrix(0, (z.lag + 1), 1),
            matrix(cset.exp, N2, 1)
        )
        cset.bls <- rbind(
            matrix(0, (z.lag + 1), 1),
            matrix(cset.bls, N2, 1)
        )
    } else {
        td <- rbind(
            matrix(0, (z.lag + 1), 1),
            td,
            matrix(0, z.lead, 1)
        )
        cset.sup <- rbind(
            matrix(0, (z.lag + 1), 1),
            matrix(cset.sup, N2, 1),
            matrix(0, z.lead, 1)
        )
        cset.avg <- rbind(
            matrix(0, (z.lag + 1), 1),
            matrix(cset.avg, N2, 1),
            matrix(0, z.lead, 1)
        )
        cset.exp <- rbind(
            matrix(0, (z.lag + 1), 1),
            matrix(cset.exp, N2, 1),
            matrix(0, z.lead, 1)
        )
        cset.bls <- rbind(
            matrix(0, (z.lag + 1), 1),
            matrix(cset.bls, N2, 1),
            matrix(0, z.lead, 1)
        )
    }

    return(
        list(
            td = td,
            cset.sup = cset.sup,
            cset.avg = cset.avg,
            cset.exp = cset.exp,
            cset.bls = cset.bls
        )
    )
}
