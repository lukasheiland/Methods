##########################################################################################
# Binomial GLM ---------------------------------------------------------------
##########################################################################################


# k/N-Codierung -----------------------------------------------------------
fit_knb <- lme4::glmer(cbind(count_alive, count_dead) ~ gradient + (1 | block), family = binomial)


# Dispersion sensitive family  --------------------------------------------
# betabinomial has a parameter for variance

fit_knb <- glmmTMB::glmmTMB(cbind(count_alive, count_dead) ~ gradient + (1 | block),
                 family = betabinomial,
                 dispformula = ~ gradient)
                 # WARNING: there is treatment coding for the varianc





# Note regarding residual checks ------------------------------------------
# In 0/1 data residual problems can only occur in group structures
# (with one draw they are necessarily binomial!)
# -> group for residual checks

## 1. residuals vs. fitted
# res_mm <- DHARMa::simulateResiduals(m_m)
# res_mm2 <- DHARMa::recalculateResiduals(res_mm, group = TD$ID, aggregateBy = sum)
# plot(res_mm)
# plot(res_mm2)
# testDispersion(res_mm)


## 2. always residuals vs. predictor!
# DHARMa::plotResiduals(TD$Visit, res_mm$scaledResiduals)