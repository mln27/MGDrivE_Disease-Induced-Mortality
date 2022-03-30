# https://rdrr.io/cran/deSolve/man/events.html
# https://tpetzoldt.github.io/deSolve-forcing/deSolve-forcing.html
# https://www.rdocumentation.org/packages/deSolve/versions/1.7/topics/events

# this is how deSolve handles events
# https://github.com/cran/deSolve/blob/dcae22b4400deb28f6c01a9c9d8a1f5fe686cc53/R/checkevents.R
## =============================================================================
## EVENTS in a data.frame
## =============================================================================
library(deSolve)

# model: rate of change set to 0
eventmod <- function(t,var,parms) {
  list(dvar = rep(0,2))
}

yini <- c(v1 = 1, v2 = 2)
times <- seq(0,10, by=0.1)

#
eventdat <- data.frame(var = c("v1", "v2", "v2", "v1"), time = c(1,1,5,9) ,
  value = c(1,2,3,4), method =c("add", "mult","rep","add"))

eventdat

out <- vode(func=eventmod, y=yini, times=times, parms=NULL,
  events=list(data=eventdat))
plot(out,type="l")

#
eventdat <- data.frame(var = c(rep("v1",10),rep( "v2", 10)),
  time = c(1:10,1:10), value = runif(20), method =rep("add", 20))

eventdat
out <- ode(func=eventmod, y=yini, times=times, parms=NULL,
  events=list(data=eventdat))

plot(out,type="l")

## =============================================================================
## EVENTS in a function
## =============================================================================

# model: rate of change v1 = 0, v2 consumed at first-order rate
eventmod <- function(t,var,parms) {
   list(c(0,-0.5*var[2]))
}


# events: add 1 to v1,multiply v2 with 0.8
eventfun <- function(t,y,parms){
  with (as.list(y),{
    v1 <- v1+1
    v2 <- 5*runif(1)
    return(c(v1,v2))
  })
}



eventfun <- function(t,y,parms){
  y["v2"] <- 2*y["v1"]
  y["v1"] <- 0
  return(y)
}
eventdat <- data.frame(var = "v1", time = 1:9,
  value = 1, method =c("add"))

# the eventdat seems to get ignored - see the deSolve link above, it's because
#  they have nested returns, and if "func" is found, a return is done before data
#  is processed

yini <- c(v1 = 1, v2 = 2)
times <- seq(0,10, by=0.1)

out <- ode(func=eventmod, y=yini, times=times, parms=NULL,
  events=list(data=eventdat, func=eventfun, time=1:9) )
plot(out,type="l")

## =============================================================================
## EVENTS triggered by a root function
## =============================================================================

# derivative: simple first-order decay
func <- function(t, y, pars) {
  return(list(-0.1*y))
}

# event triggered if state variable =0.5
rootfun <- function (t, y, pars) {
  return(y-0.5)
}

# sets state vaiable = 1
eventfun <- function(t, y, pars) {
  return(y=1)
}

yini <- 2

times <- seq(0,100,0.1)

# uses lsodar to solve; root =TRUE specifies that the event is triggered by
# a root.
out <- lsodar(times=times, y=yini, func = func, parms=NULL,
  events=list(func = eventfun, root = TRUE),
  rootfun = rootfun)

plot(out,type="l")
