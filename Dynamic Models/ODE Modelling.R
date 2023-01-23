library("dplyr")  # for data manipulation.
library("tidyr") # for data tidying, reshape2 functionality
# mutate(), select(), filter(), summarise(), arrange(

library("deSolve") # Solving Initial Value Differential Equations
library("simecol")

### Implementation of a model
# 1. das Modell, also die Differentialgleichungen in computerlesbarer Form,
# 2. die Parameter (Konstanten) des Modells,
# 3. Anfangswerte (oder Startwerte) für die Zustandsgrößen, z.B. Startabundanz oder Nährstoffkon- zentration am Anfang der Simulation,
# 4. Randbedingungen, die die ausdrücken, wie sich die Modell-Umwelt ändert (Werden solche Um- welteinflüsse nicht betrachtet, spricht man von einem autonomen System.),
# 5. eine Simulationszeit (Für welche Zeitpunkte soll das Modell simuliert werden?),
# 6. ein Lösungsverfahren, z.B. lsoda.



## Formulierung des zu integrierenden Modells (dx/dt = µ * X => X = X_0 * e^(µ*t)) als R-Funktion

exp_growth <-               # alternativ: new("odeModel", main = )
  function(t, x, p) {
    µ <- p["µ"]
    dx.dt <- µ * x+t
    list(c(dx.dt)) # return: list, whose first element is a vector containing the derivatives of y with respect to time, ? and whose next elements are global values that are required at each point in times ???
  }

parametervektor <- c(µ = 0.3)
xstart <- c(x1 = 1) # nur Startwert für die einzige abhängige x1
times <- seq(0, 10, c(delta.t= 0.2)) # Der Vektor times enthält die Zeitpunkte mit der externen Schrittweite delta.t, für die Simulationsergebnisse ausgegeben werden sollen (die eigentliche Integrations- schrittweite wird intern von lsoda festgelegt).

# Solve the differential equation
out <- lsoda(xstart, times, exp_growth, parametervektor) # simulation
# Solver for ODE:
# lsoda(y = start-values,
#       times = times-vector-where-first-is-t_0,
#       func = name of function(t, y, parms,...),     ## where t is the current time point, y current estimate of the ODE system variables
#                     ## The return value of func should be a list, whose first element is a vector containing the derivatives of y with respect to time, and whose next elements are global values that are required at each point in times. The derivatives must be specified in the same order as the state variables y.
#
#       parms = parameter-vector-used-in-func
#)
plot(out)

out.discrete <- as.data.frame(out)
plot(out.discrete$time, out.discrete$x1)


## Resource limited growth model
rl_growth <- function(t, x, parms) {
  S <- x[1] # resource ("x1")
  X <- x[2] # abundance ("x2")
  
  S.max <- parms["S.max"]
  µ.max <- parms["µ.max"]
  Y <-  parms["Y"]
 
  µ <- µ.max * (S/S.max)
  
  dS.dt <- -µ * (1/Y) * X
  dX.dt <- µ * X
  list(c(dS.dt, dX.dt))
  }

rl.parms <- c(S.max = 8, µ.max = 0.5, Y = 20)
rl.sim <- lsoda(y = c(resource = 5, abundance = 10), # start values for two independent (is that how you call them in ODEs?) variables
      times = seq(0, 50, 0.5),
      func = rl_growth,
      parms = rl.parms)
plot(rl.sim)

## Offenes System, Durchflusssystem
# Nährlösung -> Wachstum -> Export

# dS/dt = D*S_0 - D*S - µ*(1/Y)*X # dS/dt = Import − Export − Aufnahme durch Organismen
# dX/dt = µX - D*X # dX/dt = Wachstum - Export, wobei einfach exp Wachstum, Verdünnungsrate (welcher Anteil wird pro Zeiteinheit ersetzt, In und Out) D = Q * V^-1 [m^3/d * m^-3 = 1/d] mit Q = V/t

chemostat_model <- function(t, x, p) {
  S <- x[1] # resource ("x1")
  X <- x[2] # abundance ("x2")
  
  S.0 <- p["S.0"]
  µ <- p["µ"]
  Y <-  p["Y"]
  D <- p["D"] # Dilution rate 1/s
  
  
  # µ könnte auch Substratkonzentrationsabhängig sein, z.B. Monod-Kinetik, in der es sich einem µ_max annähert
  # µ <- v.max*S/(k.max+S)
   µ <- 5*S/(100+S) # Parameter
  
  dS.dt <- D*S.0 - D*S - µ*(1/Y)*X
  dX.dt <- µ*X - D*X
  
  list(c(dS.dt, dX.dt))
}

chem.parms <- c(S.0 = 100,
              µ = 0.2,
              Y = 20,
              D = 0.1)
chem.sim <- lsoda(y = c(resource = 0, abundance = 0.001),
                times = seq(0, 50, 0.5),
                func = chemostat_model,
                parms = chem.parms)
plot(chem.sim)
# Permutation
par.m <- sapply(chem.parms, FUN = function(x) rnorm(10, mean = x, sd = x/10))
# irgendwie so könnte man Exerimentelle Ergebniswerte erhalten. Zu müde.
# mat <- apply(par.m, 1, FUN = function(x) lsoda(y = c(resource = 0, abundance = 0.001),
#                                             times = seq(0, 50, 0.5),
#                                             func = chemostat_model,
#                                             parms = x))
