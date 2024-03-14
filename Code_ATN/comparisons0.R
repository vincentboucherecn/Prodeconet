########################
######
####  Produces the raw data needed to produce Figure 6 and Figure F.1 (Isolated sectors)
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
shuffleid <- sample(1:m) # ids of suffled externalities (only used if shuffle=T)
maxcaught <- 5

# get available networks/steady-states
lstfiles <- list.files("./data/")
numbs <- sort(parse_number(lstfiles[grep("net",lstfiles)]))

####################
## Simulation function (for each network)
####################

simul <- function(nnet){
  
  use <- T
  net <- readRDS(paste("data/net",nnet,".Rda",sep="")) # load network
  Tlev <- trolev(net) # compute trophic levels
  
  basals <- as.numeric(rowSums(net)==0) # flag for basal species
  vert <- readRDS(paste("data/vert",nnet,".Rda",sep="")) # load vertebrate flag
  btot0 <- readRDS(paste("data/Btot",nnet,".Rda",sep="")) # get initial stocks
  
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
  # solve the population dynamics
  ###
  
  binit <- btot0 # initial value is end of baseline
  out <- ode(y=binit,times=0:tmax,func=popdyn,parms=parameters, method="ode45", rtol = 1e-8, atol = 1e-8)
  Btot <- out[(nrow(out)-1),2:ncol(out)] # final biomass
  Btot[Btot<exth] <- 0
  Bad <- F
  if (max(abs(unlist(popdyn(1,Btot,parameters))))>1e-8){
    Bad <- T
  }
  ###
  # computes Jacobian
  ###
  wrap <- function(b){
    return(b+unlist(popdyn(1,b,parameters)))
  }
  jm <- numDeriv::jacobian(func=wrap,x=Btot,method="simple")
  jm[Btot==0,] <- 0
  jm[, Btot==0] <- 0
  if (any(is.na(jm)) | any(abs(jm)==Inf) | use==F | Bad ){
    return(NA)
  } else {
    ###
    # computes main values (prop 6)
    ###
    
    ev <- eigen(jm)
    eval <- ev$values
    evec <- (ev$vectors)
    
    db <- data.frame(reval=Re(eval),imval=Im(eval),abs=abs(eval),Btot,binit,nnet,catches0 = catches(binit,parprod), catches=catches(Btot,parprod),vert=vert,diag=diag(jm))
    
    return(db)
  }
}


####################
## Define list of TFP, remove production network, and simulate for each network and each TFP level
####################


alist <- seq(from=0.00003, to=0.0009, by=0.00003)
W0 <- vector("list",length(alist))
W <- matrix(0,n,n)
parprod[[4]] <- W
for (i in 1:length(alist)){
  parprod[[2]] <- rep(alist[i],n)
  outsim <- mclapply(numbs,simul,mc.cores=32L,mc.preschedule = FALSE)
  db <- do.call(rbind,outsim)
  W0[[i]] <- db
  save.image("tot0.RData")
}



