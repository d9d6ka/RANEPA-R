# Carrion-i-Silvestre and Sansó 2006 OBES

require(zeallot)

## `[` <- function(...) base::`[`(..., drop = FALSE) # nolint

coi_kpss <- function(y, x, model, tb, k2, cri) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    t <- nrow(y)

    if (model < 0 & model > 6)
        stop("ERROR: Try to specify the deterministic component again")

    if (cri[1] == 0) {
        if (model == 0)
            xt <- x
        else if (1 <= model & model <= 4) {
            deter <- determi(model, t, tb)
            xt <- cbind(deter, x)
        }
        else if (model == 5) {
            deter <- determi(1, t, tb)
            xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
            xt <- cbind(deter, x, xdu)
        }
        else if (model == 6) {
            deter <- determi(4, t, tb)
            xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
            xt <- cbind(deter, x, xdu)
        }

        beta <- solve(t(xt) %*% xt) %*% t(xt) %*% y
        u <- y - xt %*% beta
        t_b <- tb
    }
    else if (cri[1] == 1) {
        k <- cri[2]
        bic_op <- 100000000
        for (i in k:1) {
            c(beta, u, bic, t_b) %<-% dols(y, x, model, k, k, tb)
            if (bic < bic_op) {
                bic_op <- bic
                beta_op <- beta
                t_b_op <- t_b
                u_op <- u
            }
        }
        u <- u_op
        beta <- beta_op
        t_b <- t_b_op
    }

    sg <- alrvr(u)
    tests <- solve((t^2) * sg) * apply(apply(u, 2, cumsum)^2, 2, sum)

    return(
        list(
            beta = beta,
            tests = tests,
            u = u,
            t_b = t_b
        )
    )
}

determi <- function(model, t, tb) {
    du <- rbind(
        matrix(data = 0, nrow = tb, ncol = 1),
        matrix(data = 1, nrow = t - tb, ncol = 1)
    )
    dt <- rbind(
        matrix(data = 0, nrow = tb, ncol = 1),
        matrix(data = 1:(t - tb), nrow = t - tb, ncol = 1)
    )

    if (model == 1) {
        xt <- cbind(
            matrix(data = 1, nrow = t, ncol = 1),
            du
        )
    }
    else if (model == 2) {
        xt <- cbind(
            matrix(data = 1, nrow = t, ncol = 1),
            du,
            matrix(data = 1:t, nrow = t, ncol = 1)
        )
    }
    else if (model == 3) {
        xt <- cbind(
            matrix(data = 1, nrow = t, ncol = 1),
            matrix(data = 1:t, nrow = t, ncol = 1),
            dt
        )
    }
    else if (model == 4) {
        xt <- cbind(
            matrix(data = 1, nrow = t, ncol = 1),
            du,
            matrix(data = 1:t, nrow = t, ncol = 1),
            dt
        )
    }
    else {
        stop("ERROR: Try to specify the deterministic component again")
    }

    return(xt)
}

alrvr <- function(e) {
    if (!is.matrix(e)) e <- as.matrix(e)

    t <- nrow(e)
    k <- 0.8
    a <- solve(t(e[1:(t - 1)]) %*% e[1:(t - 1)]) %*% t(e[1:(t - 1)]) %*% e[2:t]
    l <- min(
        1.1447 * (4 * a^2 * t / ((1 + a)^2 * (1 - a)^2))^(1 / 3), # nolint
        1.1447 * (4 * k^2 * t / ((1 + k)^2 * (1 - k)^2))^(1 / 3)  # nolint
    )
    l <- trunc(l)
    lrv <- as.numeric((t(e) %*% e) / t)
    for (i in 1:l) {
        w <- (1 - i / (l + 1))
        lrv <- lrv + 2 * as.numeric(t(e[1:(t - i)]) %*% e[(1 + i):t]) * w / t
    }
    return(lrv)
}

lagn <- function(x, i) {
    if (!is.matrix(x)) x <- as.matrix(x)
    t <- nrow(x)
    k <- ncol(x)
    if (i > 0)
        return(
            rbind(
                matrix(data = NA, nrow = i, ncol = k),
                x[1:(t - i), , drop = FALSE]
            )
        )
    else
        return(
            rbind(
                x[(1 + abs(i)):t, , drop = FALSE],
                matrix(data = NA, nrow = abs(i), ncol = k)
            )
        )
}

dols <- function(y, x, model, klags, kleads, tb) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    t <- nrow(y)

    streg <- x
    d_streg <- streg[2:t, , drop = FALSE] - streg[1:(t - 1), , drop = FALSE]
    d_streg_step <- d_streg
    d_streg_r <- d_streg[(t - 1):1, , drop = FALSE]
    d_streg_r_step <- d_streg_r

    for (i in 1:klags) {
        d_streg <- cbind(
            d_streg,
            lagn(d_streg_step, i)
        )
    }

    for (i in 1:kleads) {
        d_streg_r <- cbind(
            d_streg_r,
            lagn(d_streg_r_step, i)
        )
    }

    if (klags != 0 & kleads != 0) {
        lags <- d_streg
        leads <- d_streg_r[(t - 1):1, (ncol(streg) + 1):(ncol(d_streg_r)), drop = FALSE]
        ll <- cbind(lags, leads)[(klags + 1):(t - 1 - kleads), , drop = FALSE]
    }
    else if (klags != 0 & kleads == 0) {
        lags <- d_streg
        ll <- lags[(klags + 1):(t - 1), , drop = FALSE]
    }
    else if (klags == 0 & kleads != 0) {
        lags <- d_streg
        leads <- d_streg_r[(t - 1):1, (ncol(streg) + 1):(ncol(d_streg_r)), drop = FALSE]
        ll <- cbind(lags, leads)[1:(t - 1 - kleads), , drop = FALSE]
    }
    else if (klags == 0 & kleads == 0) {
        ll <- d_streg
    }

    if (model == 0) {
        xreg <- cbind(
            streg[(klags + 2):(t - kleads), , drop = FALSE],
            ll
        )
    }
    else if (model >= 1 & model <= 4) {
        deter <- determi(model, t, tb)
        xreg <- cbind(
            deter[(klags + 2):(t - kleads), , drop = FALSE],
            streg[(klags + 2):(t - kleads), , drop = FALSE],
            ll
        )
    }
    else if (model == 5) {
        deter <- determi(1, t, tb)
        xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
        xreg <- cbind(
            deter[(klags + 2):(t - kleads), , drop = FALSE],
            streg[(klags + 2):(t - kleads), , drop = FALSE],
            xdu[(klags + 2):(t - kleads), , drop = FALSE],
            ll
        )
    }
    else if (model == 6) {
        deter <- determi(4, t, tb)
        xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
        xreg <- cbind(
            deter[(klags + 2):(t - kleads), , drop = FALSE],
            streg[(klags + 2):(t - kleads), , drop = FALSE],
            xdu[(klags + 2):(t - kleads), , drop = FALSE],
            ll
        )
    }

    beta <- solve(t(xreg) %*% xreg) %*% t(xreg) %*%
        y[(klags + 2):(t - kleads), 1, drop = FALSE]
    e <- y[(klags + 2):(t - kleads), 1, drop = FALSE] - xreg %*% beta
    s2 <- as.numeric(t(e) %*% e) / (nrow(xreg) - ncol(xreg))
    t_beta <- sweep(beta, 1, sqrt(diag(s2 * solve(t(xreg) %*% xreg))))

    bic <- log(s2) + ncol(xreg) * log(nrow(xreg)) / nrow(xreg)

    return(
        list(
            beta = beta,
            e = e,
            bic = bic,
            t_beta = t_beta
        )
    )
}

valors <- function(model, lambda, k) {
    load("data/coi_end_vc.data")
    m_vc <- coi_end_vc[[model]][[k]]
    if (0 < lambda & lambda <= 0.15)
        v_vc <- m_vc[, 1, drop = FALSE]
    else if (0.15 < lambda & lambda <= 0.25)
        v_vc <- m_vc[, 2, drop = FALSE]
    else if (0.25 < lambda & lambda <= 0.35)
        v_vc <- m_vc[, 3, drop = FALSE]
    else if (0.35 < lambda & lambda <= 0.45)
        v_vc <- m_vc[, 4, drop = FALSE]
    else if (0.45 < lambda & lambda <= 0.55)
        v_vc <- m_vc[, 5, drop = FALSE]
    else if (0.55 < lambda & lambda <= 0.65)
        v_vc <- m_vc[, 6, drop = FALSE]
    else if (0.65 < lambda & lambda <= 0.75)
        v_vc <- m_vc[, 7, drop = FALSE]
    else if (0.75 < lambda & lambda <= 0.85)
        v_vc <- m_vc[, 8, drop = FALSE]
    else if (0.85 < lambda & lambda < 1)
        v_vc <- m_vc[, 9, drop = FALSE]
    else
        stop("ERROR! Try to specify the value of lambda again")
}
