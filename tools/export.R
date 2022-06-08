source("tools/cvals_kpss_1p.R")
source("tools/cvals_kpss_2p.R")
source("tools/cvals_gsadf.R")
source("tools/cvals_kp.R")

save(.cval_kpss_1p,
     .cval_kpss_2p,
     .cval_SADF_without_const,
     .cval_SADF_with_const,
     .cval_GSADF_without_const,
     .cval_GSADF_with_const,
     .cval_KP,
     file = "sysdata.rda")
