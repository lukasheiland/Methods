##########################################################################################
# Lotka-Volterra model for Predator--Prey system  ---------------------------------------
##########################################################################################

# Notes -------------------------------------------------------------------
# The predator--prey variant has exponential growth (while the competition variant is logistic)


# Library -----------------------------------------------------------------
library(deSolve)


# Model -------------------------------------------------------------------
calcPP <- function (time, state, par) {
  with(as.list(c(state, par)), {
    dx <- x * (alpha - beta*y) # y: state
    dy <- -y * (gamma - delta*x) # x: state
    return(list(c(dx, dy)))
  })
}


# Simulation --------------------------------------------------------------
par <- c(alpha = .5, beta = .5, gamma = .4, delta = .3)
state <- c(x = 10, y = 9)
time <- seq(0, 300, by = 1)

Sim <- ode(func = calcPP, y = state, parms = par, times = time)

matplot(Sim[,-1], type = "l", xlab = "time", ylab = "population")
legend("topright", c("Cute bunnies", "Rabid foxes"), lty = c(1,2), col = c(1,2), box.lwd = 0)


# Stan --------------------------------------------------------------------
library(rethinking)
data(Lynx_Hare_model)
cat(Lynx_Hare_model) # Nice formulation with observation model
