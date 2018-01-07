
# An Introduction to the 'dlmBE' package

<img src="BE.png" alt="Built environment" width="400" height="281">

_Simulated features of the built environment. Each cube represents the
location(s) of environment features; each orange dot represents a participant
location._


<img src="BE_circ.png" alt="Built environment" width="400" height="281">


## Fitting and understanding DL models

```R
## (x, y) positions for subjects
> head(subj.xy)
   x  y
1 18 46
2 46 47
3 27 22
4 28 26
5 45 19
6 47 48

## (x, y) positions for BE features
> head(feat.xy)
   x y
1  1 1
2  8 1
3  8 1
4 31 1
5 37 1
6 38 1
```

We begin by computing, for each participant, the radial distance to each
environmental feature. We follow the literature and count the total number
of features at each available distance on the (50 x 50) grid.

```R
count.features <- function(xy, feature.xy, radii) {
  .dist <- function(x) sqrt(sum(x^2))
  dxy <- apply(sweep(feature.xy, 2, xy), 1, .dist)
  table(cut(dxy, radii, include.lowest = TRUE))
}

## each available radius
lag <- 1:50

## count of features for each subject (row) and radius (column)
Conc <- t(apply(subj.xy, 1, count.features,
                feature.xy = feat.xy, radii = c(0, lag))
          )
```

```R
> anova(fit0, fit1)
refitting model(s) with ML (instead of REML)
Data: NULL
Models:
fit0: y ~ cr(lag, Conc)
fit1: y ~ c.age + I(c.age^2) + cr(lag, Conc)
     Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
fit0  5 672.03 688.53 -331.02   662.03
fit1  7 551.96 575.04 -268.98   537.96 124.08      2  < 2.2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


