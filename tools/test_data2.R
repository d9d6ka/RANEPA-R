Y <- rnorm(100)
Y2 <- Y
Y2[51:100] <- Y2[51:100] + 2
data <- data.frame(Y, Y2)

SSR.1 <- SSR.matrix(data$Y2, 1:100, width = 4)
res.1.2 <- segs.SSR.N.breaks(data$Y2, 1:100, 1, width = 4, SSR.data = SSR.1)
KPSS.1.2b <- KPSS.N.breaks(data$Y2,
                           NULL,
                           model = 1,
                           break.point = res.1.2$break.point,
                           const = FALSE,
                           trend = TRUE,
                           ll.init = 4,
                           corr.max = 0,
                           kernel = NULL,
                           weakly.exog = TRUE)
KPSS.1.2b <- KPSS.N.breaks(data$Y,
                           NULL,
                           model = 1,
                           break.point = res.1.2$break.point,
                           const = FALSE,
                           trend = TRUE,
                           ll.init = 4,
                           corr.max = 0,
                           kernel = NULL,
                           weakly.exog = TRUE)
KPSS.1.2boot <- KPSS.N.breaks.bootstrap(data$Y,
                                        NULL,
                                        model = 1,
                                        break.point = res.1.2$break.point,
                                        const = FALSE,
                                        trend = TRUE,
                                        ll.init = 4,
                                        corr.max = 0,
                                        kernel = NULL,
                                        weakly.exog = TRUE)
