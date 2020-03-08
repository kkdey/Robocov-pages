
types = c("Robospan-Blood-Mean", "Robospan-Mean", "pRobospan-Mean", "pRobospan-Blood-Mean")
annots = rep(c("100kb", "5kb", "ABC", "Coding", "eQTL", "5kg2", "Promoter", "Promoter2",
               "Road", "TSS", "Yoshida"), 4)
p = c(0.0048, 0.016, 0.26, 0.048, 0.057, 0.032, 0.023, 0.02, 0.0017, 0.16, 0.028, 0.061,
      0.00025, 0.0036, 4e-04, 3e-03, 0.04, 0.019, 0.017, 0.045, 0.01, 0.073, 0.0043, 0.0079,
      0.0066, 0.0037, 5.4e-05, 0.0023, 0.0016, 0.026, 0.013, 0.032, 5e-04, 2.5e-03, 1e-03, 0.036,
      1e-04, 0.17, 0.013, 0.4, 0.026, 0.28, 0.015, 0.075, 0.0063, 0.032, 0.03, 0.18)

library(qvalue)

qq = qvalue(p, pi0 = 1)
which(qq$qvalues < 0.01)

## marginally significant annotations

## Robospan-Mean+100kb
## Robospan-Mean+ABC
## pRobospan-Mean + ABC
## pRobospan-Mean + Promoter-Gazal
## pRobospan-Mean + TSS
## pRobospan-Blood-Mean + 100kb


