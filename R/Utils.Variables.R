#' @title
#' Construct determinant variables for [KPSS.1.break]
#'
#' @details
#' Procedure to compute deterministic terms
#' for KPSS with 1 structural break.
#'
#' @param model A scalar equal to
#' * 1: Model with trend, break in const,
#' * 2: Model with const and trend, break in const,
#' * 3: Model with const and trend, break in trend,
#' * 4: Model with const and trend, break in const and trend.
#' @param n.obs Number of observations.
#' @param break.point Break point.
#'
#' @return Matrix of determinant variables.
#'
#' @keywords internal
determinants.KPSS.1.break <- function(model,
                                      n.obs,
                                      break.point) {
    du <- rbind(
        matrix(data = 0, nrow = break.point, ncol = 1),
        matrix(data = 1, nrow = n.obs - break.point, ncol = 1)
    )
    dt <- rbind(
        matrix(data = 0, nrow = break.point, ncol = 1),
        matrix(
            data = 1:(n.obs - break.point),
            nrow = n.obs - break.point,
            ncol = 1
        )
    )
    const <- matrix(data = 1, nrow = n.obs, ncol = 1)
    trend <- matrix(data = 1:n.obs, nrow = n.obs, ncol = 1)

    if (model == 1) {
        xt <- cbind(
            trend,
            du
        )
    } else if (model == 2) {
        xt <- cbind(
            const,
            du,
            trend
        )
    } else if (model == 3) {
        xt <- cbind(
            const,
            trend,
            dt
        )
    } else if (model == 4) {
        xt <- cbind(
            const,
            du,
            trend,
            dt
        )
    } else {
        stop("ERROR: Try to specify the deterministic component again")
    }

    return(xt)
}


#' @title
#' Construct determinant variables for [KPSS.2.breaks]
#'
#' @details
#' Procedure to compute deterministic terms
#' for KPSS with 2 structural breaks.
#'
#' @param model A scalar equal to
#' * 1: for the AA (without trend) model,
#' * 2: for the AA (with trend) model,
#' * 3: for the BB model,
#' * 4: for the CC model,
#' * 5: for the AC-CA model,
#' * 6: for the AC-CA model,
#' * 7: for the AC-CA model.
#' @param n.obs Number of observations.
#' @param break.point Positions for the first and second structural breaks
#'            (respective to the origin which is 1).
#'
#' @return Matrix of deterministic terms.
#'
#' @keywords internal
determinants.KPSS.2.breaks <- function(model,
                                       n.obs,
                                       break.point) {
    cns <- matrix(data = 1, nrow = n.obs, ncol = 1)
    tend <- matrix(data = 1:n.obs, nrow = n.obs, ncol = 1)

    if (model == 1) {
        du1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(data = 1, nrow = n.obs - break.point[1], ncol = 1)
        )
        du2 <- rbind(
            matrix(data = 0, nrow = break.point[2], ncol = 1),
            matrix(data = 1, nrow = n.obs - break.point[2], ncol = 1)
        )
        xt <- cbind(cns, du1, du2)
    } else if (model == 2) {
        du1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(data = 1, nrow = n.obs - break.point[1], ncol = 1)
        )
        du2 <- rbind(
            matrix(data = 0, nrow = break.point[2], ncol = 1),
            matrix(data = 1, nrow = n.obs - break.point[2], ncol = 1)
        )
        xt <- cbind(cns, tend, du1, du2)
    } else if (model == 3) {
        dt1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(
                data = 1:(n.obs - break.point[1]),
                nrow = n.obs - break.point[1], ncol = 1
            )
        )
        dt2 <- rbind(
            matrix(data = 0, nrow = break.point[2], ncol = 1),
            matrix(
                data = 1:(n.obs - break.point[2]),
                nrow = n.obs - break.point[2], ncol = 1
            )
        )
        xt <- cbind(cns, tend, dt1, dt2)
    } else if (model == 4) {
        du1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(data = 1, nrow = n.obs - break.point[1], ncol = 1)
        )
        du2 <- rbind(
            matrix(data = 0, nrow = break.point[2], ncol = 1),
            matrix(data = 1, nrow = n.obs - break.point[2], ncol = 1)
        )
        dt1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(
                data = 1:(n.obs - break.point[1]),
                nrow = n.obs - break.point[1], ncol = 1
            )
        )
        dt2 <- rbind(
            matrix(data = 0, nrow = break.point[2], ncol = 1),
            matrix(
                data = 1:(n.obs - break.point[2]),
                nrow = n.obs - break.point[2], ncol = 1
            )
        )
        xt <- cbind(cns, tend, du1, dt1, du2, dt2)
    } else if (model == 5) {
        du1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(data = 1, nrow = n.obs - break.point[1], ncol = 1)
        )
        dt2 <- rbind(
            matrix(data = 0, nrow = break.point[2], ncol = 1),
            matrix(
                data = 1:(n.obs - break.point[2]),
                nrow = n.obs - break.point[2], ncol = 1
            )
        )
        xt <- cbind(cns, tend, du1, dt2)
    } else if (model == 6) {
        du1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(data = 1, nrow = n.obs - break.point[1], ncol = 1)
        )
        du2 <- rbind(
            matrix(data = 0, nrow = break.point[2], ncol = 1),
            matrix(data = 1, nrow = n.obs - break.point[2], ncol = 1)
        )
        dt1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(
                data = 1:(n.obs - break.point[1]),
                nrow = n.obs - break.point[1], ncol = 1
            )
        )
        xt <- cbind(cns, tend, du1, dt1, du2)
    } else if (model == 7) {
        du1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(data = 1, nrow = n.obs - break.point[1], ncol = 1)
        )
        dt1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(
                data = 1:(n.obs - break.point[1]),
                nrow = n.obs - break.point[1], ncol = 1
            )
        )
        dt2 <- rbind(
            matrix(data = 0, nrow = break.point[2], ncol = 1),
            matrix(
                data = 1:(n.obs - break.point[2]),
                nrow = n.obs - break.point[2], ncol = 1
            )
        )
        xt <- cbind(cns, tend, du1, dt1, dt2)
    } else {
        stop("Try to speciy another model")
    }

    return(xt)
}


#' @title
#' Deterministic terms for [KPSS.N.breaks]
#'
#' @description
#' Procedure to compute deterministic terms for KPSS with \eqn{m}
#' structural breaks.
#'
#' @details
#' **model** should be either a scalar or a vector of the same size as the
#' **break.point**. If scalar **model** will be repeated till the length of
#' **break.point** is achieved.
#'
#' @param model A scalar or vector of
#' * 1: for the break in const,
#' * 2: for the break in trend,
#' * 3: for the break in const and trend.
#' @param n.obs Number of observations.
#' @param break.point Array of structural breaks.
#' @param const,trend Include constant and trend if `TRUE`.
#'
#' @return Matrix of deterministic terms.
#'
#' @keywords internal
determinants.KPSS.N.breaks <- function(model,
                                       n.obs,
                                       break.point,
                                       const = FALSE,
                                       trend = FALSE) {
    ntb <- length(break.point)
    if (length(model) == 1) {
        model <- rep(model, ntb)
    } else if (length(model) != ntb) {
        stop("ERROR! Inconsistent sizes of model and break.point")
    }

    xt <- NULL
    if (const) {
        xt <- matrix(data = 1, nrow = n.obs, ncol = 1)
    }
    if (trend) {
        xt <- cbind(
            xt,
            matrix(data = 1:n.obs, nrow = n.obs, ncol = 1)
        )
    }

    for (i in 1:ntb) {
        if (model[i] == 1) {
            xt <- cbind(
                xt,
                rbind(
                    matrix(
                        data = 0,
                        nrow = break.point[i],
                        ncol = 1
                    ),
                    matrix(
                        data = 1,
                        nrow = n.obs - break.point[i],
                        ncol = 1
                    )
                )
            )
        } else if (model[i] == 2) {
            xt <- cbind(
                xt,
                rbind(
                    matrix(
                        data = 0,
                        nrow = break.point[i],
                        ncol = 1
                    ),
                    matrix(
                        data = 1:(n.obs - break.point[i]),
                        nrow = n.obs - break.point[i],
                        ncol = 1
                    )
                )
            )
        } else if (model[i] == 3) {
            xt <- cbind(
                xt,
                rbind(
                    matrix(
                        data = 0,
                        nrow = break.point[i],
                        ncol = 1
                    ),
                    matrix(
                        data = 1,
                        nrow = n.obs - break.point[i],
                        ncol = 1
                    )
                ),
                rbind(
                    matrix(
                        data = 0,
                        nrow = break.point[i],
                        ncol = 1
                    ),
                    matrix(
                        data = 1:(n.obs - break.point[i]),
                        nrow = n.obs - break.point[i],
                        ncol = 1
                    )
                )
            )
        }
    }

    return(xt)
}


#' @title
#' Generating monthly seasonal dummy variables
#'
#' @param n.obs number of observations.
#'
#' @return The matrix of values od seasonal dummies.
#'
#' @keywords internal
seasonal.dummies <- function(n.obs) {
    s1 <- c(1 - 1 / 12, rep(-1 / 12, 11))

    result <- NULL
    for (i in 0:10) {
        result <- cbind(
            result,
            c(
                rep(-1 / 12, i),
                rep(s1, length.out = n.obs - i)
            )
        )
    }

    return(result)
}
