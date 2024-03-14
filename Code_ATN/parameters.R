########################
######
####  This code contains many parameters values
######
########################


nspecies <- 30
conn <- 0.15
err <- 0.025
maxtry <- 1000
exth <- 1e-6
fishth <- 1e-3
tmax <- 10000
extpr <- 0.4
delta <- 0.3
nsim <- 1000
weightcor <- 0.05
extraction <- F # presence of externalities?
adjust <- T # production is a function of stocks?
shuffle <- F # shuffle externalities ? ## NOT USED
posext <- F # revert sign of externalities? ### NOT USED
m <- nspecies
n <- 30
lshare <- 0.25 # labor share
B <- matrix(0,n,m) # matrix of ecosystem services
Bpar <- 0.1 # weights of ecosystem services

## matrix for the production network
W <- matrix(rPareto(n*n,0.01,1),n,n) 
normW <- rowSums(W)
normW[normW==0] <- 1
W <- diag(1/normW)%*%W

W <- (1-lshare)*W

net <- matrix(0,nspecies,nspecies)
vert <- rep(0,nspecies)
parprod <- list(
  rep(0,n), #mu <- parprod[[1]]
  rep(0.12,n), #a <- parprod[[2]]
  B, #B <- parprod[[3]]
  W, #W <- parprod[[4]]
  rep(0,n), #lbd <- parprod[[5]]
  rep(1/n,n), #q <- parprod[[6]]
  NA, # Ccoefs <- parprod[[7]]
  NA
)
parameters <- list(
  rep(1,nspecies), # K <- parameters[[1]] #Carrying capacities
  rep(1,nspecies), # r <- parameters[[2]] # max growth rates
  c(1,0.314,0.88), # ar, ax (invertebrate), ax (vertebrate) # a <- parameters[[3]]
  c(10,10), #c(8,4),   # y (invertebrate), y(vertebrate)
  matrix(1/nspecies, nspecies,nspecies), # w <- parameters[[5]] # relative consumption rates
  0.2, # b0 <- parameters[[6]] # half-saturation
  1.2, # h <- parameters[[7]] # Hill coefficient
  rep(0,nspecies), # c <- parameters[[8]] # Predator interference
  c(0.85,0.45), # when eats non-basal, when eats basal # e <- parameters[[9]] # assimilation efficiency
  100, # Z <- parameters[[10]] # body size ratio
  net, # trophic network
  vert, # vertebrates
  parprod
)
