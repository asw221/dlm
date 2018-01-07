
# An Introduction to the 'dlmBE' package

The goal of this package is to provide researchers with a convenient interface
to fit and summarize distributed lag models (DLMs) using the R programming
language. Among so many modeling strategies in the toolbox, users might reach
for DLMs when they have an outcome that is related to distance-profiled
predictors through some unknown function. A typical goal could then be to
learn about the shape of that distance-profiled response function. For
example, this type of model might be applied when a researcher wants to learn:

 - (Epidemiology) the relationship between health outcomes and proximity of
   subjects to certain environmental features
 - (Cognitive Neuroscience) the shape of the blood-oxygen response
   in the time following some sort of stimulus in functional MR images of the
   brain
 - ...

And many more. For the purposes of this walkthrough, we will focus on just
one type of example from our own past research: health outcomes and
subjects' proximity to features of the built environment
(Baek, J, *et al*, 2016, 2017).

We simulate data on a (50 x 50) grid representing an imaginary cityscape.
Within our city we simulate a number of built environment features with a
mild spatial correlation, and *N* = 200 subjects with homes distributed
uniformly over our plot of land. We have data in the form of (*x*, *y*)
coordinate pairs for each subject and environmental feature, and descriptive
information about the gender and age of each subject.


<img src="BE.png" alt="Built environment" width="400" height="281">

_Simulated features of the built environment. Each cube represents the
location(s) of environment features; each orange dot represents a participant
location._

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

> table(female)
female
  0   1
 99 101

> summary(age)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  18.00   32.00   47.00   46.56   61.00   75.00
```

We're interested in the case where the features we've built into our cityscape
have some average measurable impact on participant outcomes, but this impact
changes as a function of distance. Imagine each "feature" is, say, a fast-food
restaurant and we want to learn about what kind of effect living near this
type of convenience has on subjects' body-mass index (BMI). One analytical
approach could be to count how many restaurants are within a given radial
distance of each subject's home and include that count in a regression model
as a predictor of our outcome (BMI; `y`).


<img src="BE_circ.png" alt="Built environment" width="400" height="281">

_Same environment as above, but focused on the first participant in our
data set and the number of features within some units of her home. Radii shown
are 4, 8, 12, and 16 units._


In practice, we probably don't often know the appropriate radius to pick for
this type of problem (or if different sub-populations react differently over
different radii, etc.) although extant literature or domain-specific knowledge
may lead us to a few reasonable guesses. The DLM framework provides an
alternative solution to this type of problem when the analyst is comfortable
making the extra assumption that the underlying function of distance is
continuous, or approximately so. In this example, DLMs free the user from
the responsability of selecting a single radius or distance-threshold.
Instead, the analyst can input many radii and rely on the DLM to infer
the shape of a continuous function that links them all.

Although there are multiple different options to allow users to estimate
arbitrary functions of distance, we focus on the use of splines as a
flexible and interpretable semi-parametric method.
Our implementation in the `dlmBE` package relies on the excellent `lme4`
package to penalize the spline terms using mixed effects modeling and
provide numerically stable results, even for large numbers of radii.


## Fitting and understanding DL models

```R
count.features <- function(xy, feature.xy, radii) {
  .dist <- function(x) sqrt(sum(x^2))  # Euclidean distance
  dxy <- apply(sweep(feature.xy, 2, xy), 1, .dist)
  table(cut(dxy, radii, include.lowest = TRUE))
}

lag <- 1:50  # each available radius

## count of features for each subject (row) and radius (column)
Conc <- t(apply(subj.xy, 1, count.features,
                feature.xy = feat.xy, radii = c(0, lag)))

## basic model--only DL term
fit0 <- dlm(y ~ cr(lag, Conc))
```

We begin by computing, for each participant, the radial distance to each
environmental feature. We follow our prior work and count the total number
of features at each available distance on the (50 x 50) grid. In general,
we advocate an analysis strategy of starting simple and gradually allowing
for more complexity, so we fit an initial model with only one DL function
of distance and number of fast-food locations.

```R
> summary(fit0)
Linear mixed model fit by REML ['dlMod']
Formula: y ~ cr(lag, Conc)

REML criterion at convergence: 936.7

Scaled residuals:
     Min       1Q   Median       3Q      Max
-2.85065 -0.68897  0.00298  0.78527  2.26968

Random effects:
 Groups        Name   Variance  Std.Dev.
 cr(lag, Conc) (mean) 1.013e-08 0.0001007
 Residual             5.584e+00 2.3630757
Number of obs: 200, groups:  cr(lag, Conc), 48

Fixed effects:
            Estimate Std. Error t value
(Intercept)   24.069      6.731   3.576

Correlation of Fixed Effects:
<0 x 0 matrix>
```


<img src="fit0_resids.png" alt="fit0 diagnostics" width="600" height="268">

_Quick residual diagnostics for the model with only one DL function.
There appears to be a non-constant variance pattern, and age is clearly
correlated with the residuals from this fit._

```R
qplot(fitted(fit0), residuals(fit0)) +
  geom_hline(yintercept = 0, col = "gray40")

qplot(age, residuals(fit0)) +
  geom_hline(yintercept = 0, col = "gray40")
```

Standard `summary()` methods are available for `dlmBE::dlMod` objects
(the output type of the `dlm()` function), but the printout is designed
mostly for easy interpretation of fixed effects covariates. Here, it's a
little more informative to explore model summaries graphically.
`dlmBE` uses the `ggplot2` package for its default plotting methods, so
we continue in that vein for exploratory and diagnostic data visualization.


```R
## Model including a function of age
fit1 <- dlm(y ~ c.age + I(c.age^2) + cr(lag, Conc))

## Model including interaction between DL term and gender
fit2 <- dlm(y ~ c.age + I(c.age^2) + cr(lag, Conc) * female)
```


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
