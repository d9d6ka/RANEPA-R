## Procedure to compute deterministic terms
## for KPSS with 2 structural breaks.
## Sintax is as follows:
##       detho = dekpss2(model,N,tb1,tb2); # nolint
## INPUT:
##    model     Scalar, determines the deterministic model that
##              is considered under the null hypothesis. If:
##                 model = 1   Model AAn,
##                       = 2   Model AA,
##                       = 3   Model BB,
##                       = 4   Model CC,
##                       = 5   Model AC-CA.
##    N         Number of observations.
##    tb1       Position for the first structural break
##              (respective to the origin which is 1).
##    tb2       Position for the second structural break
##              (respective to the origin which is 1).
## OUTPUT:
##    detho     Matrix of deterministic terms.
determi_kpss_2p <- function(model, N, tb1, tb2) {
    cns <- matrix(data = 1, nrow = N, ncol = 1)
    tend <- matrix(data = 1:N, nrow = N, ncol = 1)

    if (model == 1) {
        du1 <- rbind(
            matrix(data = 0, nrow = tb1, ncol = 1),
            matrix(data = 1, nrow = N - tb1, ncol = 1)
        )
        du2 <- rbind(
            matrix(data = 0, nrow = tb2, ncol = 1),
            matrix(data = 1, nrow = N - tb2, ncol = 1)
        )
        xt <- cbind(cns, du1, du2)
    }
    else if (model == 2) {
        du1 <- rbind(
            matrix(data = 0, nrow = tb1, ncol = 1),
            matrix(data = 1, nrow = N - tb1, ncol = 1)
        )
        du2 <- rbind(
            matrix(data = 0, nrow = tb2, ncol = 1),
            matrix(data = 1, nrow = N - tb2, ncol = 1)
        )
        xt <- cbind(cns, tend, du1, du2)
    }
    else if (model == 3) {
        dt1 <- rbind(
            matrix(data = 0, nrow = tb1, ncol = 1),
            matrix(data = 1:(N - tb1), nrow = N - tb1, ncol = 1)
        )
        dt2 <- rbind(
            matrix(data = 0, nrow = tb2, ncol = 1),
            matrix(data = 1:(N - tb2), nrow = N - tb2, ncol = 1)
        )
        xt <- cbind(cns, tend, dt1, dt2)
    }
    else if (model == 4) {
        du1 <- rbind(
            matrix(data = 0, nrow = tb1, ncol = 1),
            matrix(data = 1, nrow = N - tb1, ncol = 1)
        )
        du2 <- rbind(
            matrix(data = 0, nrow = tb2, ncol = 1),
            matrix(data = 1, nrow = N - tb2, ncol = 1)
        )
        dt1 <- rbind(
            matrix(data = 0, nrow = tb1, ncol = 1),
            matrix(data = 1:(N - tb1), nrow = N - tb1, ncol = 1)
        )
        dt2 <- rbind(
            matrix(data = 0, nrow = tb2, ncol = 1),
            matrix(data = 1:(N - tb2), nrow = N - tb2, ncol = 1)
        )
        xt <- cbind(cns, tend, du1, dt1, du2, dt2)
    }
    else if (model == 5) {
        du1 <- rbind(
            matrix(data = 0, nrow = tb1, ncol = 1),
            matrix(data = 1, nrow = N - tb1, ncol = 1)
        )
        du2 <- rbind(
            matrix(data = 0, nrow = tb2, ncol = 1),
            matrix(data = 1, nrow = N - tb2, ncol = 1)
        )
        dt1 <- rbind(
            matrix(data = 0, nrow = tb1, ncol = 1),
            matrix(data = 1:(N - tb1), nrow = N - tb1, ncol = 1)
        )
        xt <- cbind(cns, tend, du1, dt1, du2)
    }
    else
        stop("Try to speciy another model")

    return(xt)
}
