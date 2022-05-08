Y <- rnorm(100)
Y2 <- Y
Y2[51:100] <- Y2[51:100] + 2
data <- data.frame(Y, Y2)

ssr.1 <- ssr_matrix(data$Y2, 1:100, width = 4)
res.1.2 <- ssr_partition_mp(data$Y2, 1:100, 1, width = 4, ssr_data = ssr.1)
kpss.1.2b <- kpss_known_mp(data$Y2,
                           NULL,
                           model = 1,
                           break.point = res.1.2$break_point,
                           const = FALSE,
                           trend = TRUE,
                           ll.init = 4,
                           corr.max = 0,
                           kernel = NULL,
                           weakly.exog = TRUE)
kpss.1.2b <- kpss_known_mp(data$Y,
                           NULL,
                           model = 1,
                           break.point = res.1.2$break_point,
                           const = FALSE,
                           trend = TRUE,
                           ll.init = 4,
                           corr.max = 0,
                           kernel = NULL,
                           weakly.exog = TRUE)
kpss.1.2boot <- bootstrap_kpss_mp(data$Y,
                                  NULL,
                                  model = 1,
                                  break.point = res.1.2$break_point,
                                  const = FALSE,
                                  trend = TRUE,
                                  ll.init = 4,
                                  corr.max = 0,
                                  kernel = NULL,
                                  weakly.exog = TRUE)
