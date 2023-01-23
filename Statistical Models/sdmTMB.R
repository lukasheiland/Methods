# Library -----------------------------------------------------------------

library(sdmTMB)


## Example

# SPDE Mesh -----------------------------------------------------------------
skimr::skim(pcod)
pcod_spde <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 25) # a coarse mesh for example speed
plot(pcod_spde)

## for spatial autocorrelation: approximating a Gaussian Field, i.e. a continuos coordinate space of gaussians,
## with stochastic partial differential equations (SPDE) by dividing the coordinate space into certain triangles with associated errors
## that allow for estimating a sparse correlation matrix
## See: https://becarioprecario.bitbucket.io/spde-gitbook/ch-intro.html#sec:spde




# Fit ---------------------------------------------------------------------
# Tweedie:
## the distribution can be parameterized by a a mean, dispersion, and shape parameter "p". For p=0 it is normal, for p=1 poisson and p=2 is gamma.
## The "interesting" distributions with point mass at zero and continuous positive support are when 1<p<2 (aka "compound poisson").
## rtweedie(1000, mu = 14, phi = 10, power = 1.5)
m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
            data = pcod_2011, time = "year", spde = pcod_spde, family = tweedie(link = "log"))
print(m)
tidy(m, conf.int = TRUE)
tidy(m, effects = "ran_par", conf.int = TRUE)


# Bernoulli:
pcod_binom <- pcod_2011
pcod_binom$present <- ifelse(pcod_binom$density > 0, 1L, 0L)
m_bin <- sdmTMB(present ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
                data = pcod_binom, time = "year", spde = pcod_spde,
                family = binomial(link = "logit"))
print(m_bin)

# Fit a spatial-only model (by not specifying `time`):
m <- sdmTMB(
  density ~ depth_scaled + depth_scaled2, data = pcod_2011,
  spde = pcod_spde, family = tweedie(link = "log"))
print(m)

# Gaussian:
pcod_gaus <- subset(pcod_2011, density > 0 & year >= 2013)
pcod_spde_gaus <- make_mesh(pcod_gaus, c("X", "Y"), cutoff = 30)
m_pos <- sdmTMB(log(density) ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
                data = pcod_gaus, time = "year", spde = pcod_spde_gaus)
print(m_pos)

# With splines via mgcv.
# Make sure to pre-specify an appropriate basis dimension (`k`) since
# the smoothers are not penalized in the current implementation.
# See ?mgcv::choose.k
m_gam <- sdmTMB(log(density) ~ 0 + as.factor(year) + s(depth_scaled, k = 4),
                data = pcod_gaus, time = "year", spde = pcod_spde_gaus)
print(m_gam)

# With IID random intercepts:
# Simulate some data:
set.seed(1)
x <- runif(500, -1, 1)
y <- runif(500, -1, 1)
loc <- data.frame(x = x, y = y)
spde <- make_mesh(loc, c("x", "y"), n_knots = 50, type = "kmeans")
s <- sdmTMB_sim(x = x, y = y, betas = 0, time = 1L,
                phi = 0.1, range = 1.4, sigma_O = 0.2, sigma_E = 0, mesh = spde)
s$g <- gl(50, 10)
iid_re_vals <- rnorm(50, 0, 0.3)
s$observed <- s$observed + iid_re_vals[s$g]

# Fit it:
m <- sdmTMB(observed ~ 1 + (1 | g), spde = spde, data = s)
print(m)
tidy(m, "ran_pars", conf.int = TRUE) # see tau_G
theta <- as.list(m$sd_report, "Estimate")
plot(iid_re_vals, theta$RE)



## Vignette  -----------------------------------------------------------------
# https://pbs-assess.github.io/sdmTMB/articles/spatial-trend-models.html

# SPDE Mesh -----------------------------------------------------------------
skimr::skim(pcod)
pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 12)
plot(pcod_spde)

## for spatial autocorrelation: approximating a Gaussian Field, i.e. a continuos coordinate space of gaussians,
## with stochastic partial differential equations (SPDE) by dividing the coordinate space into certain triangles with associated errors
## that allow for estimating a sparse correlation matrix
## See: https://becarioprecario.bitbucket.io/spde-gitbook/ch-intro.html#sec:spde


# Fit ---------------------------------------------------------------------
m1 <- sdmTMB(density ~ 1, data = pcod,
             spde = pcod_spde, family = tweedie(link = "log"),
             spatial_trend = TRUE, time = "year",
             spatial_only = TRUE)

m2 <- sdmTMB(density ~ 1, data = pcod,
             spde = pcod_spde, family = tweedie(link = "log"),
             spatial_trend = TRUE, time = "year",
             spatial_only = FALSE)

m3 <- sdmTMB(density ~ 1, data = pcod,
             spde = pcod_spde, family = tweedie(link = "log"),
             spatial_trend = TRUE, time = "year",
             spatial_only = FALSE, fields = "AR1")

d <- pcod
d$residuals1 <- residuals(m1)
d$residuals2 <- residuals(m2)
d$residuals3 <- residuals(m3)

qqnorm(d$residuals1);abline(a = 0, b = 1)
qqnorm(d$residuals2);abline(a = 0, b = 1)
qqnorm(d$residuals3);abline(a = 0, b = 1)

plot_map_point <- function(dat, column = "est") {
  ggplot(dat, aes_string("X", "Y", colour = column)) +
    geom_point() +
    facet_wrap(~year) +
    coord_fixed()
}


plot_map_point(d, "residuals1") + scale_color_gradient2()
plot_map_point(d, "residuals2") + scale_color_gradient2()
plot_map_point(d, "residuals3") + scale_color_gradient2()

sd1 <- as.data.frame(summary(TMB::sdreport(m1$tmb_obj)))
sd2 <- as.data.frame(summary(TMB::sdreport(m2$tmb_obj)))
sd3 <- as.data.frame(summary(TMB::sdreport(m3$tmb_obj)))

r1 <- m1$tmb_obj$report()
r2 <- m2$tmb_obj$report()
r3 <- m3$tmb_obj$report()


sd3$Estimate[row.names(sd3) == "ar1_phi"]
sd3$Estimate[row.names(sd3) == "ar1_phi"] +
  c(-2, 2) * sd3$`Std. Error`[row.names(sd3) == "ar1_phi"]

plot_map_raster <- function(dat, column = "est") {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed() +
    scale_fill_viridis_c()
}


p1 <- predict(m1, newdata = qcs_grid)
p2 <- predict(m2, newdata = qcs_grid)
p3 <- predict(m3, newdata = qcs_grid)

plot_map_raster(filter(p1, year == 2003), "zeta_s")
plot_map_raster(filter(p2, year == 2003), "zeta_s")
plot_map_raster(filter(p3, year == 2003), "zeta_s")

plot_map_raster(p1, "est")
plot_map_raster(p2, "est")
plot_map_raster(p3, "est")


# And we can look at just the spatiotemporal random effects for models 2 and 3:
plot_map_raster(p2, "est_rf") + scale_fill_gradient2()
plot_map_raster(p3, "est_rf") + scale_fill_gradient2()


plot_map_raster(filter(p1, year == 2003), "omega_s")
plot_map_raster(filter(p2, year == 2003), "omega_s")
plot_map_raster(filter(p3, year == 2003), "omega_s")



# Regensburg --------------------------------------------------------------
# library(glmmTMB)
library(sdmTMB)
# devtools::install_github(repo = "florianhartig/EcoData", subdir = "EcoData", 
#                          dependencies = T, build_vignettes = T)
library(EcoData)
library(sf)
library(mapview)

plants <- plantcounts
plants$agrarea_scaled <- scale(plants$agrarea)
plants$tk <- as.factor(plants$tk)
plants_sf <- sf::st_as_sf(plants, coords = c('lon', 'lat'), crs = st_crs("+proj=longlat +ellps=bessel +towgs84=606,23,413,0,0,0,0 +no_defs"))

mapview(plants_sf["agrarea"])
mapview(plants_sf["richness"], map.types = "OpenTopoMap")


fit <-  glmmTMB(richness ~ agrarea_scaled + offset(log(area)), family = nbinom1, data = plants_sf)
summary(fit)

SPDE <- make_mesh(plants, c("lon", "lat"), n_knots = 50, type = "kmeans")
plot(SPDE)
fit2 <-  sdmTMB(richness ~ scale(agrarea) + (1 | tk), family = nbinom2(), data = plants, spde = SPDE, fields = "IID")
summary(fit2)
P <- predict(fit2)
# est: Estimate in link space (everything is in link space)
# est_non_rf: Estimate from everything that isn't a random field
# est_rf: Estimate from all random fields combined

plot(residuals(fit2, type = "sim") ~ P$est) # ?sdmTMB::residuals.sdmTMB


library(ggplot2)
plot_map <- function(dat, column = "est") {

}

# ggplot(plants_sf) +
#   geom_sf() +
#   geom_raster(data = P, aes_string("lon", "lat", fill = "est"))
