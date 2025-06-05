########################
######
####  Calibrates the TFP levels for the production network based on the ones for the isolated sectors.
######
########################

library(igraph)
library(deSolve)
library(numDeriv)
library(ggplot2)
library(readr)
library(Pareto)
library(parallel)
library(stringr)

rm(list=ls())
set.seed(3212)

####################
## Load list of functions and parameter values
####################

source("function.R")
source("parameters.R")

scripttime <- Sys.time()
scripttime <- str_replace_all(scripttime," ","_")
file.copy("comparisons.R",paste(scripttime,".txt",sep=""))

####################
## Set additional parameters
####################

extraction <- T # presence of externalities?
maxcaught <- 5
tmax <- 5000 # update max time


# get available networks/steady-states
lstfiles <- list.files("./data/")
numbs <- sort(parse_number(lstfiles[grep("net",lstfiles)]))

####################
## Simulation function (for each network)
####################

simul <- function(nnet,a0){
  
  use <- T
  net <- readRDS(paste("data/net",nnet,".Rda",sep="")) # load network
  Tlev <- trolev(net) # compute trophic levels
  
  basals <- as.numeric(rowSums(net)==0) # flag for basal species
  vert <- readRDS(paste("data/vert",nnet,".Rda",sep="")) # load vertebrate flag
  btot0 <- readRDS(paste("data/Btot",nnet,".Rda",sep="")) # get initial stocks

  ###
  # Production parameters
  ###
    
  parprod[[8]] <- btot0
  
  if (sum(btot0>fishth)<2){
    use = F
    spcatch <- 1
  } else {
    spcatch <- c(1:nspecies)[btot0>fishth]
  }
  spcatch <- spcatch[1:min(length(spcatch),maxcaught)]
  
  parprod[[7]] <- c(spcatch)
  
  parameters[[11]] <- net
  parameters[[12]] <- vert
  parameters[[13]] <- parprod
  
  ###
  # Calibrate based on catches for isolated sectors
  ###
  
  c0cal <- readRDS(paste("data/catch0_",nnet,"_",a0,".Rda",sep=""))
  c0cal <- c0cal[1,]
  
  binit <- btot0 # initial value is end of baseline
  objmin <- function(al){
    parprodloc <- parprod
    parprodloc[[2]] <- rep(al,n)
    loc_catches <- catches(binit,parprodloc)
    return(sum(abs(loc_catches-c0cal)))
  }
  
  out <- optimize(objmin,c(1e-8,1))
  
  return(out$minimum)
}

####################
## Define list of TFP and simulate for each network and each TFP level
####################

alist0 <- seq(from=0.01,to=0.105, by=0.005)/4 ## list of TFP for isolated sectors
alist <- rep(NA,length(alist0))

for (i in 1:length(alist)){
  simloc <- function(numbs) simul(numbs,alist0[i])
  outsim <- mclapply(numbs,simloc,mc.cores=32L,mc.preschedule = FALSE)
  alist[i] <- mean(unlist(outsim))
  print(i)
}
saveRDS(alist,"alistcal.Rda")
