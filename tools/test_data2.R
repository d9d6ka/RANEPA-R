Y <- rnorm(100)
Y2 <- Y
Y2[51:100] <- Y2[51:100] + 2
data <- data.frame(Y, Y2)

ssr.1 <- ssr.matrix(data$Y2, 1:100, width = 4)
res.1.2 <- segs.ssr.N.breaks(data$Y2, 1:100, 1, width = 4, ssr.data = ssr.1)
kpss.1.2b <- kpss.N.breaks(data$Y2,
                           NULL,
                           model = 1,
                           break.point = res.1.2$break.point,
                           const = FALSE,
                           trend = TRUE,
                           ll.init = 4,
                           corr.max = 0,
                           kernel = NULL,
                           weakly.exog = TRUE)
kpss.1.2b <- kpss.N.breaks(data$Y,
                           NULL,
                           model = 1,
                           break.point = res.1.2$break.point,
                           const = FALSE,
                           trend = TRUE,
                           ll.init = 4,
                           corr.max = 0,
                           kernel = NULL,
                           weakly.exog = TRUE)
kpss.1.2boot <- kpss.N.breaks.bootstrap(data$Y,
                                        NULL,
                                        model = 1,
                                        break.point = res.1.2$break.point,
                                        const = FALSE,
                                        trend = TRUE,
                                        ll.init = 4,
                                        corr.max = 0,
                                        kernel = NULL,
                                        weakly.exog = TRUE)
