
# An Introduction to the 'dlmBE' package

The goal of this package is to provide researchers with a convenient interface
to fit and summarize distributed lag models (DLMs) using the R programming
language. Among so many modeling strategies in the toolbox, users might reach
for DLMs when they have an outcome that is related to distance-profiled
predictors through some unknown function. A typical goal could then be to
learn about the shape of that distance-profiled response function. For
example, this type of model might be applied when a researcher wants to learn:

 - (Epidemiology) the relationship between health outcomes and proximity of
   subjects to certain environmental features.
 - (Cognitive Neuroscience) the shape of the blood-oxygen response
   in the time folling some sort of stimulus in functional MR images of the
   brain.
 - ...

And many more. For the purposes of this walkthrough, we will focus on just
one type of example from our own past research: health outcomes and
subjects' proximity to features of the built environment
(Baek, J *et al.* 2016, 2017).

We simulate data on a (50 x 50) grid representing an imaginary cityscape.
Within our city we simulate a number of built environment features with a
mild spatial correlation, and *N* = 200 subjects with homes distributed
uniformly over our plot of land.


<img src="BE.png" alt="Built environment" width="400" height="281">

_Simulated features of the built environment. Each cube represents the
location(s) of environment features; each orange dot represents a participant
location._


<img src="BE_circ.png" alt="Built environment" width="400" height="281">

_Same environment as above Radii are 4, 8, 12, and 16 units._


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

lag <- 1:50  # each available radius

## count of features for each subject (row) and radius (column)
Conc <- t(apply(subj.xy, 1, count.features,
                feature.xy = feat.xy, radii = c(0, lag)))

## basic model--only DL term
fit0 <- dlm(y ~ cr(lag, Conc))

## Model including a function of age
fit1 <- dlm(y ~ c.age + I(c.age^2) + cr(lag, Conc))

## Model including interaction between DL term and gender
fit2 <- dlm(y ~ c.age + I(c.age^2) + cr(lag, Conc) * female)
```


## Fitting and understanding DL models

```R
> anova(fit0, fit1, fit2)
refitting model(s) with ML (instead of REML)
Data: NULL
Models:
fit0: y ~ cr(lag, Conc)
fit1: y ~ c.age + I(c.age^2) + cr(lag, Conc)
fit2: y ~ c.age + I(c.age^2) + cr(lag, Conc) * female
     Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)
fit0  5 923.08 939.57 -456.54   913.08
fit1  7 896.68 919.77 -441.34   882.68  30.391      2  2.516e-07 ***
fit2 11 618.79 655.07 -298.39   596.79 285.897      4  < 2.2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## References
Baek J, Sanchez BN, Berrocal VJ, & Sanchez-Vaznaugh EV (2016) Epidemiology
27(1):116-24. ([PubMed](https://www.ncbi.nlm.nih.gov/pubmed/26414942))

Baek J, Hirsch JA, Moore K, Tabb LP, et al. (2017) Epidemiology 28(3):403-11.
([PubMed](https://www.ncbi.nlm.nih.gov/pubmed/28145983))
