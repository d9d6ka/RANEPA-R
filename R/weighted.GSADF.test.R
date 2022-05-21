# weighted.GSADF.test - Weighted SADF test (HLZ, 2018).
weighted.GSADF.test <- function(y, r0 = 0.01 + 1.8/sqrt(length(y)), const = TRUE,
                                B = 200, alpha = 0.05, cpu = 4, urs = TRUE,
                                seed = round(10^4*sd(y))) {
  result <- list()
  result$y <- y
  N <- length(y)
  result$r0 <- r0
  result$const <- const
  result$B <- B
  result$alpha <- alpha
  result$cpu <- cpu
  result$urs <- urs
  result$seed <- seed

  # Find supBZ_value.
  supBZ.model <- supBZ.statistic(y, r0)
  sigma_sq <- supBZ.model$sigma_sq
  result$sigma_sq <- sigma_sq
  result$BZ_values <- supBZ.model$BZ_values
  result$supBZ_value <- supBZ.model$supBZ_value

  # Auxiliary function for parallel.
  calc_sadf_supBZ_value <- function(cpu, y, r0, const, B, seed, urs, sigma_sq) {
    set.seed(seed + cpu)
    N <- length(y)
    result <- list()
    for(i in 1:B) {
      y_star <- cumsum(c(0, rnorm(N-1)*diff(y)))
      if (urs == TRUE) {
        gsadf.model <- GSADF.test(y_star, r0, const)
        result$sadf_value[i] <- gsadf.model$sadf_value
      }
      supBZ.model <- supBZ.statistic(y_star, r0, sigma_sq)
      result$supBZ_value[i] <- supBZ.model$supBZ_value
    }
    return(result)
  }

  # Do parallel.
  cl <- makePSOCKcluster(rep("localhost", cpu))
  clusterExport(cl, c("GSADF.test", "ADF.test", "supBZ.statistic",
                      "GSADF_cr_values_with_const",
                      "GSADF_cr_values_without_const"))
  sadf_supBZ_bootstrap_values <- parLapply(cl, 1:cpu, calc_sadf_supBZ_value, y,
                                           r0, const, B, seed, urs, sigma_sq)
  stopCluster(cl)

  # Get sadf_supBZ_bootstrap_values.
  result$supBZ_bootstrap_values <- c()
  for(i in 1:cpu) {
      result$supBZ_bootstrap_values <- c(result$supBZ_bootstrap_values,
        sadf_supBZ_bootstrap_values[[i]]$supBZ_value)
  }

  # Find critical value.
  result$supBZ_cr_value <- as.numeric(quantile(result$supBZ_bootstrap_values,
                                               1-alpha))

  # A union of rejections strategy.
  if (urs == TRUE) {

    # Find sadf_value.
    gsadf.model <- GSADF.test(y, r0, const)
    result$t.values <- gsadf.model$t.values
    result$sadf_value <- gsadf.model$sadf_value

    # Get sadf_supBZ_bootstrap_values.
    result$sadf_bootstrap_values <- c()
    for(i in 1:cpu) {
      result$sadf_bootstrap_values <- c(result$sadf_bootstrap_values,
                                        sadf_supBZ_bootstrap_values[[i]]$sadf_value)
    }

    # Find critical value.
    result$sadf_cr_value <- as.numeric(quantile(result$sadf_bootstrap_values,
                                                1-alpha))

    # Calculate U value.
    result$U_value <- max(result$sadf_value,
                          result$sadf_cr_value / result$supBZ_cr_value *
                            result$supBZ_value)

    # Find U_bootstrap_values.
    result$U_bootstrap_values <- c()
    for(b in 1:(cpu*B)) {
      result$U_bootstrap_values[b] <- max(result$sadf_bootstrap_values[b],
                                         result$sadf_cr_value /
                                           result$supBZ_cr_value *
                                           result$supBZ_bootstrap_values[b])
    }

    # Find critical value.
    result$U_cr_value <- as.numeric(quantile(result$U_bootstrap_values, 1-alpha))

    result$p_value <- round(sum(result$U_bootstrap_values >
                                  result$U_value)/(cpu*B), 4)

    result$is_explosive <- ifelse(result$U_value > result$U_cr_value, 1, 0)
  } else {
    result$p_value <- round(sum(result$supBZ_bootstrap_values >
                                  result$supBZ_value)/(cpu*B), 4)

    result$is_explosive <- ifelse(result$supBZ_value > result$supBZ_cr_value, 1, 0)
  }

  return(result)
}