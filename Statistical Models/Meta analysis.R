##########################################################################################
# Meta analysis ---------------------------------------------------------------
##########################################################################################

## 1. 'Averaging'
outcome <- rnorm(5)
ci_studyresults <-  runif(5, 0, 5)
weight <- 1/ci_studyresults
summary(metastudy <- lm(outcome ~ 1, weights = weight))

## 2. Bias correction
# A. File drawer problem: bias towards not pulbicising non-significant results
#   -> Correction
# B. Publication bias: bias towards publication of bigger effects (~ journal)
#   -> Account for by analysing the drop in small effects. Tannenbaumplot with quantile regression.