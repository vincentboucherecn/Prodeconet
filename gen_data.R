########################
######
####  This code generates the baseline (no fishing) ecosystems (networks).
####  It assumes that the folder "data" is present in the current path
######
########################

library(igraph)
library(Pareto)
library(deSolve)
library(numDeriv)
library(ggplot2)
library(readr)
library(parallel)
library(rslurm)
rm(list=ls())
set.seed(3212)

####################
## Load list of functions and parameter values
####################

source("function.R")
source("parameters.R")

####################
## Simulation function (for each network)
####################

simul <- function(dummy){
  print(dummy)
  ###
  # draw network and computes main properties
  ###
  
  net <- niche(nspecies,conn,err,maxtry)
  Tlev <- trolev(net)
  
  basals <- as.numeric(rowSums(net)==0)
  vert <- as.numeric(Tlev>=3)*rbinom(nspecies,1,0.6)
  
  parameters[[11]] <- net
  parameters[[12]] <- vert
  ###
  # solve the population dynamics
  ###
  
  binit <- runif(nspecies,0.05,1)
  #out <- lsoda(y=binit,times=0:tmax,func=popdyn,parms=parameters, rtol = 1e-8, atol = 1e-8)
  out <- ode(y=binit,times=0:tmax,func=popdyn,parms=parameters, method="ode45", rtol = 1e-8, atol = 1e-8)
  Btot <- out[(nrow(out)-1),2:ncol(out)] # final biomass
  
  reason <- ""
  Bad <- F
  if (max(abs(unlist(popdyn(1,Btot,parameters))))>1e-8){
    binit <- colMeans(out[(nrow(out)-1000):nrow(out),2:ncol(out)])
    out <- ode(y=binit,times=0:tmax,func=popdyn,parms=parameters, method="ode45", rtol = 1e-8, atol = 1e-8)
    Btot <- out[(nrow(out)-1),2:ncol(out)] # final biomass
    
    if (max(abs(unlist(popdyn(1,Btot,parameters))))>1e-8){
      Bad <- T
      reason <- c(reason,"No steady-state found")
    }
  }
  Btot[Btot<exth] <- 0
  if (!any(vert==1 & Btot>0)){
    Bad <- T
    reason <- c(reason,"All vertebrates are extinct")
  }
  if (sum(Btot>0)<20){
    Bad <- T
    reason <- c(reason,"fewer than 20 species")
  }
  if (max(Btot)>=1e20){
    Bad <- T
    reason <- c(reason,"upper-bound")
  }
  nettest <- net[Btot>0,Btot>0]
  if (!check2(nettest)){
    Bad <- T
    reason <- c(reason,"equilibrium network fails connectivity tests")
  }
  ###
  # computes Jacobian
  ###
  wrap <- function(b){
    return(b+unlist(popdyn(1,b,parameters)))
  }
  jm <- numDeriv::jacobian(func=wrap,x=Btot,method="simple")
  if (any(is.na(jm)) | any (abs(jm)==Inf) | Bad){
    write.table(reason, file = paste("data/reason",dummy,".txt",sep=""))
  } else {
      saveRDS(net,file=paste("data/net",dummy,".Rda",sep=""))
      saveRDS(binit,file=paste("data/binit",dummy,".Rda",sep=""))
      saveRDS(Btot,file=paste("data/Btot",dummy,".Rda",sep=""))
      saveRDS(vert,file=paste("data/vert",dummy,".Rda",sep=""))
  }
}

####################
## Simulate baseline data
####################

nc <- 32L
outsim <- mclapply(1:nsim,simul,mc.cores=nc,mc.preschedule=F)
