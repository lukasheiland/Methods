##########################################################################################
# SEM  ---------------------------------------
##########################################################################################

# Notes -------------------------------------------------------------------
## Unterschied zwischen piecewiseSEM und lavaan ist:
##  piecewiseSEM (nichts anderes als mehrer Regressionen und ihre Tabellen zusammengefasst): decomposing the whole graph into d-seperation claims, with individual models for base pairs, Gesamt-p-Wert (kann der Gesamt-Graph rejected werden?)
##      but also possibilities like random effects etc.
##  lavaan: fittet eine Kovarianzmatrix, kann deshalb latente Variablen!
##    but bo random effects etc.
##    fittet by default 'freie Kovarianzen' für Knoten von denen nur Pfeile ausgehen, ist numerisch das Gleiche wie kausale Kovarianz, aber wird nicht kausal interpretiert. Die sind drin, weil man ja nicht weiß wie die losen Enden korreliert sind.
##    Alles innerhalb der Annahme multivariat normalverteilt.

# Library -----------------------------------------------------------------
library(lavaan)
library(piecewiseSEM)
# vignette("piecewiseSEM")

# Assumptions -----------------------------------------------------------------
dat <- data.frame(x1 = runif(50), y1 = runif(50), y2 = runif(50), y3 = runif(50))

## Hypothesis
# DAG: x1 -> y2 -> y3, x1  -> y2  -> y1
model <- psem(lm(y1 ~ x1, dat),
              lm(y1 ~ y2, dat),
              lm(y2 ~ x1, dat),
              lm(y3 ~ y1, dat))

model <- psem(lm(y1 ~ x1 + y2, dat),
              lm(y2 ~ x1, dat),
              lm(y3 ~ y1, dat))

summary(model)
coefs(model, standardize = 'scale')
basisSet(model) # Acquires the set of independence claims–or the 'basis set'–for use in evaluating the goodness-of-fit for piecewise structural equation models.

dSep(model) # The tests of directed separation evaluate this hypothesis: that we are justified in excluding relationships.
fisherC(model) #
