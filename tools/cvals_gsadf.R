hbm.sim <- function(N,
                    mu = 0,
                    tau,
                    sigma0,
                    sigma1,
                    x0 = 0,
                    vs = "shift",
                    tau.break,
                    delta.break) {
    e <- rnorm(N)

    if (vs == "shift") {
        sigma <- sigma0 + (sigma1 - sigma0) *
            ifelse((1:N) / N > tau, 1, 0)
    } else if (vs == "double_shift") {
        sigma <- sigma0 + (sigma1 - sigma0) *
            ifelse(((1:N) / N > 0.4) & ((1:N) / N <= 0.6), 1, 0)
    } else if (vs == "lst") {
        sigma <- sigma0 + (sigma1 - sigma0) *
            1 / (1 + exp(-50 * ((1:N) / N - 0.5)))
    } else if (vs == "trend") {
        sigma <- sigma0 + (sigma1 - sigma0) *
            ((1:N) / N)
    }
    eps <- sigma * e
    u <- c(x0)
    for (i in 2:(N + 1)) {
        if (i <= floor(tau.break[1] * N)) {
            u[i] <- u[i - 1] + eps[i - 1]
        } else if (i > floor(tau.break[1] * N) & i <= floor(tau.break[2] * N)) {
            u[i] <- (1 + delta.break[1]) * u[i - 1] + eps[i - 1]
        } else if (i > floor(tau.break[2] * N) & i <= floor(tau.break[3] * N)) {
            u[i] <- (1 - delta.break[2]) * u[i - 1] + eps[i - 1]
        } else {
            u[i] <- u[i - 1] + eps[i - 1]
        }
    }
    y <- mu + u
    if (y[floor(tau.break[2] * N)] < y[floor(tau.break[1] * N)]) {
        y <- (-1) * y
    }
    return(
        list(
            N = N,
            mu = mu,
            sigma0 = sigma0,
            sigma1 = sigma1,
            x0 = x0,
            tau.break = tau.break,
            delta.break = delta.break,
            tau = tau,
            sigma = sigma,
            eps = eps,
            u = u,
            y = y
        )
    )
}

# Calculate finite sample critical values.
####
set.seed(42)
# N_sim - Number of observations.
N_sim <- 1000
# tau_sim - Moment of volatility shift.
tau_sim <- 0.5
# sigma0, sigma1.
sigma0 <- 1
sigma1 <- 1
# delta_b1_sim - Bubble size.
delta_b1_sim <- 0
# NN_sim - Number of simulations.
NN_sim <- 999

progress.bar <- txtProgressBar(max = NN_sim, style = 3)
progress <- function(n) setTxtProgressBar(progress.bar, n)

# tau.break - Bubble period.
.cval_SADF_without_const <- c()
for (i in 1:NN_sim) {
    y <- hbm.sim(
        N = N_sim, tau = tau_sim, sigma0 = sigma0, sigma1 = sigma1,
        tau.break = c(0.4, 0.6, 1), delta.break = c(delta_b1_sim, 0)
    )

    model <- STADF.test(y$y, r0 = 0.1, const = FALSE, add.p.value = FALSE)

    .cval_SADF_without_const[i] <- model$stadf.value

    progress(i)
}

.cval_SADF_with_const <- c()
for (i in 1:NN_sim) {
    y <- hbm.sim(
        N = N_sim, tau = tau_sim, sigma0 = sigma0, sigma1 = sigma1,
        tau.break = c(0.4, 0.6, 1), delta.break = c(delta_b1_sim, 0)
    )

    model <- STADF.test(y$y, r0 = 0.1, const = TRUE, add.p.value = FALSE)

    .cval_SADF_with_const[i] <- model$stadf.value

    progress(i)
}

# tau.break - Bubble period.
.cval_GSADF_without_const <- c()
for (i in 1:NN_sim) {
    y <- hbm.sim(
        N = N_sim, tau = tau_sim, sigma0 = sigma0, sigma1 = sigma1,
        tau.break = c(0.4, 0.6, 1), delta.break = c(delta_b1_sim, 0)
    )

    model <- GSTADF.test(y$y, r0 = 0.1, const = FALSE, add.p.value = FALSE)

    .cval_GSADF_without_const[i] <- model$gstadf.value

    progress(i)
}

.cval_GSADF_with_const <- c()
for (i in 1:NN_sim) {
    y <- hbm.sim(
        N = N_sim, tau = tau_sim, sigma0 = sigma0, sigma1 = sigma1,
        tau.break = c(0.4, 0.6, 1), delta.break = c(delta_b1_sim, 0)
    )

    model <- GSTADF.test(y$y, r0 = 0.1, const = TRUE, add.p.value = FALSE)

    .cval_GSADF_with_const[i] <- model$gstadf.value

    progress(i)
}
