library(mlr)
library(keras)
library(tidyverse)

data = EcoData::titanic
data = data %>% 
  select(pclass, survived, sex, age, sibsp, parch, fare, body)

data[,-2] = missRanger::missRanger(data[,-2])



learner = makeLearner("classif.ranger", importance = "impurity", predict.type = "prob")
param = makeParamSet(makeIntegerParam("mtry", 1L, 5L), 
                     makeIntegerParam("min.node.size", 2L, 50L))
task = makeClassifTask(data = data.frame(data), target = "survived", positive = "1")
task = oversample(task, rate = 2)

learnerTune = makeTuneWrapper(learner, resampling = cv10, measures = acc, par.set = param, control = makeTuneControlRandom(maxit = 40L))

resample = resample(learnerTune, task, cv5, measures = auc, models = TRUE, keep.pred = TRUE, extract = mlr::getTuneResult)


