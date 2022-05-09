#' @title
#' Andrews-Kim (2006) test
#'
#' @description
#' Test for structural break at the sample end
#'
#' @details
#' See Andrews and Kim (2006) for the detailed description.
#'
#' @param eq Base model formula
#' @param m Post-break period length
#' @param dataset Source of the data
#'
#' @return The list, containing
#' \itemize{
#' \item **m**
#' \item estimated values of P- and R-tests.
#' \item sequences of auxiliary statistics \eqn{P_j} and \eqn{R_j}.
#' \item the corresponding p-values.
#' }
#'
#' @export
andrews_kim_test <- function(eq, m, dataset) {
    result <- list(m = m)

    dep_var <- eq[[2]]
    reg_var <- all.vars(eq)
    reg_var <- reg_var[reg_var != dep_var]

    N <- dim(dataset)[1] - m # nolint

    A <- matrix(NA, m, m) # nolint
    for (i in 1:m) for (j in 1:m) A[i, j] <- min(i, j)

    tmp_model <- lm(formula = eq, data = dataset)
    tmp_resid <- as.matrix(tmp_model$residuals[(N + 1):(N + m)])

    result$P <- drop(t(tmp_resid) %*% tmp_resid) # nolint
    result$R <- drop(t(tmp_resid) %*% A %*% tmp_resid) # nolint

    result$Pj <- seq(0, 0, length.out = N - m + 1) # nolint
    result$Rj <- seq(0, 0, length.out = N - m + 1) # nolint

    for (j in seq_len(N - m + 1)) {
        low <- j
        high <- j + ceiling(m / 2) - 1
        tmp_data <- dataset[1:N, ][-(low:high), ] # nolint

        tmp_model <- lm(formula = eq, data = tmp_data)

        tmp_resid <- dataset[[dep_var]][j:(j + m - 1)] -
            predict(tmp_model, newdata = dataset[j:(j + m - 1), ])
        tmp_resid <- as.matrix(tmp_resid)

        result$Pj[j] <- drop(t(tmp_resid) %*% tmp_resid)
        result$Rj[j] <- drop(t(tmp_resid) %*% A %*% tmp_resid)
    }

    result$p_value <- (1 / (N - m + 1)) * sum(I(result$P <= result$Pj))
    result$r_value <- (1 / (N - m + 1)) * sum(I(result$R <= result$Rj))

    return(result)
}
