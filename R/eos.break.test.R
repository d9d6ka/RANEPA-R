#' @title
#' Andrews-Kim (2006) test
#'
#' @description
#' Test for structural break at the end of the sample.
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
eos.break.test <- function(eq, m, dataset) {
    result <- list(m = m)

    dep.var <- eq[[2]]
    reg.vars <- all.vars(eq)
    reg.vars <- reg.vars[reg.vars != dep.var]

    N <- dim(dataset)[1] - m

    A <- matrix(NA, m, m)
    for (i in 1:m) for (j in 1:m) A[i, j] <- min(i, j)

    temp.model <- lm(formula = eq, data = dataset)
    temp.resid <- as.matrix(temp.model$residuals[(N + 1):(N + m)])

    result$P <- drop(t(temp.resid) %*% temp.resid)
    result$R <- drop(t(temp.resid) %*% A %*% temp.resid)

    result$Pj <- seq(0, 0, length.out = N - m + 1)
    result$Rj <- seq(0, 0, length.out = N - m + 1)

    for (j in seq_len(N - m + 1)) {
        low <- j
        high <- j + ceiling(m / 2) - 1
        temp.data <- dataset[1:N, ][-(low:high), ]

        temp.model <- lm(formula = eq, data = temp.data)

        temp.resid <- dataset[[dep.var]][j:(j + m - 1)] -
            predict(temp.model, newdata = dataset[j:(j + m - 1), ])
        temp.resid <- as.matrix(temp.resid)

        result$Pj[j] <- drop(t(temp.resid) %*% temp.resid)
        result$Rj[j] <- drop(t(temp.resid) %*% A %*% temp.resid)
    }

    result$p_value <- (1 / (N - m + 1)) * sum(I(result$P <= result$Pj))
    result$r_value <- (1 / (N - m + 1)) * sum(I(result$R <= result$Rj))

    return(result)
}
