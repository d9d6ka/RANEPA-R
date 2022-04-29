## Procedure to minimise the SSR
## Sintaxis:
##  {e, Tb_vec} = min_ssr(y,model);
## Input:
##      y       (Tx1)-vector of time series
##      model   = 1 for the AA (without trend) model
##              = 2 for the AA (with trend) model
##              = 3 for the BB model
##              = 4 for the CC model
##              = 5 for the AC-CA model
##  Output:
##      e       (Tx1) vector of estimated OLS residuals
##      Tb_vec  (2x1) vector with the estimated break points
#' @importFrom zeallot %<-%
min_ssr <- function(y, model) {
    if (!is.matrix(y)) y <- as.matrix(y)

    t <- nrow(y)

    r_min <- 0
    ssr_min <- 1000000
    tb1_min <- 0
    tb2_min <- 0

    if (1 <= model & model <= 4) {
        for (i in 2:(t - 4)) {
            for (j in (i + 2):(t - 2)) {
                z <- determi_kpss_2p(model, t, i, j)
                c(b, r, p) %<-% olsqr(y, z)
                ssr <- as.numeric(t(r) %*% r)
                if (ssr < ssr_min) {
                    r_min <- r
                    ssr_min <- ssr
                    tb1_min <- i
                    tb2_min <- j
                }
            }
        }
    }
    else if (model == 5) {
        for (i in 2:(t - 4)) {
            for (j in (i + 2):(t - 2)) {
                z <- determi_kpss_2p(model, t, i, j)
                c(b, r, p) %<-% olsqr(y, z)
                ssr <- as.numeric(t(r) %*% r)
                if (ssr < ssr_min) {
                    r_min <- r
                    ssr_min <- ssr
                    tb1_min <- i
                    tb2_min <- j
                }
            }
        }
        for (j in 2:(t - 4)) {
            for (i in (j + 2):(t - 2)) {
                z <- determi_kpss_2p(model, t, i, j)
                c(b, r, p) %<-% olsqr(y, z)
                ssr <- as.numeric(t(r) %*% r)
                if (ssr < ssr_min) {
                    r_min <- r
                    ssr_min <- ssr
                    tb1_min <- j
                    tb2_min <- i
                }
            }
        }
    }

    return(
        list(u = r_min, tb1 = tb1_min, tb2 = tb2_min)
    )
}

olsqr <- function(y, x) {
    b <- solve(t(x) %*% x) %*% t(x) %*% y
    p <- x %*% b
    r <- y - p
    return(
        list(b = b, r = r, p = p)
    )
}
