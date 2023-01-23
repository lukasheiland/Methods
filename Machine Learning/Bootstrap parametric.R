fit = lm(log(survival$surv)~survival$dose)

res = resid(fit)

f = fitted(fit)

surv.r.mle = data.frame(f,res)
surv.r.fun = function(data) coef(lm(log(data$surv)~data$dose))

surv.r.sim = function(data, mle){
  data$surv = exp(mle$f+sample(mle$res,nrow(mle),T))
  return(data)
}

surv.r.boot = boot(survival, surv.r.fun, R = 999, sim = "parametric",
                   ran.gen = surv.r.sim, mle = surv.r.mle)
