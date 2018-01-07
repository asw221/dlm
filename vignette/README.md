
# An Introduction to the 'dlmBE' package

<img src="BE.png" alt="Built environment" style="width: 200px;"/>
*Simulated features of the built environment. Each cube represents the
location(s) of environment features; each orange dot represents a participant
location.*



## Fitting and understanding DL models

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


