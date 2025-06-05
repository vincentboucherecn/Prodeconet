########################
######
####  This code contains a list of useful functions
######
########################


################################
####### Trophic network ########
################################

## check if trophic network is valid

plotnet <- function(A,x,vert,ext,clean){
  if (!clean){
    gr <- graph_from_adjacency_matrix(A)
    V(gr)$color <- vert
    V(gr)$color[c(1:nrow(A))%in%ext] <- 2
    V(gr)$size <- exp(x)*15
    V(gr)$shape <- "circle"
    V(gr)$shape[x==0] <- "square"
    plot(gr,edge.arrow.size=0.6)
  } else {
    ext <- c(1:nrow(A))%in%ext
    A <- A[x>0,x>0]
    vert <- vert[x>0]
    ext <- ext[x>0]
    x <- x[x>0]
    gr <- graph_from_adjacency_matrix(A)
    V(gr)$color <- vert
    V(gr)$color[c(1:nrow(A))%in%ext] <- 2
    V(gr)$size <- exp(x)*15
    V(gr)$shape <- "circle"
    V(gr)$shape[x==0] <- "square"
    plot(gr,edge.arrow.size=0.6)
  }
  
}



check <- function(A,err){
  basal <- as.numeric(rowSums(A)==0)
  gr <- graph_from_adjacency_matrix(A)
  dl <- distances(gr,mode="out")
  c1 <- all(sapply(c(1:nspecies)[basal==0], function(i) any(dl[i,c(1:nspecies)[basal==1]]<Inf)))
  c2 <- is.connected(gr,mode="weak")
  c3 <- (mean(A) <= (1+err)*conn) & (mean(A) >= (1-err)*conn)
  #c4 <- any(rowSums(A)+colSums(A)>0)
  return(all(c(c1,c2,c3)))
}

check2 <- function(A){
  basal <- as.numeric(rowSums(A)==0)
  mm <- length(basal)
  gr <- graph_from_adjacency_matrix(A)
  dl <- distances(gr,mode="out")
  c1 <- all(sapply(c(1:mm)[basal==0], function(i) any(dl[i,c(1:mm)[basal==1]]<Inf)))
  c2 <- is.connected(gr,mode="weak")
  return(all(c(c1,c2)))
}

## Draw a trophic network using the Niche model of Williams and Martinez (2000).

niche <- function(nspecies,conn,err,maxtry){
  j <- 1
  while (j<=maxtry){
    pos <- runif(nspecies) # position of species
    srange <- rbeta(nspecies,1,(1-2*conn)/(2*conn)) # feeding range
    srange <- srange*pos # scale range
    center <- pmin(runif(nspecies,srange/2,pos),1-srange/2) # range center
    
    ## force lowest specie to be basal (Glaum et al. 2020)
    srange[which(rank(pos)==1)] <- 0
    lb <- center - srange/2
    ub <- center + srange/2
    mpos <- matrix(rep(pos,nspecies),nspecies,nspecies,byrow=T)
    mlb <- matrix(rep(lb,nspecies),nspecies,nspecies)
    mub <- matrix(rep(ub,nspecies),nspecies,nspecies)
    net <- matrix(as.numeric(mpos<=mub & mpos>=mlb),nspecies, nspecies) ## net[i,j] means i eats j
    if (check(net,err)){
      j <- maxtry+1
    } else {
      j <- j+1
      if (j==maxtry){
        warning("Maxtry exceeded")
        net <- matrix(NA,nspecies, nspecies)
      }
    }
  }
  return(net)
}

## Compute trophic level using Levine (1980) analytic formula

trolev <- function(net){
  diag(net) <- 0
  norma <- matrix(rep(rowSums(net),nspecies),nspecies,nspecies)
  norma[norma==0] <- 1
  invmat <- try(solve(diag(nspecies)-(net/norma)),silent=T)
  if (length(invmat)==1){
    return(NA)
  } else {
    Tvec=c(rowSums(invmat))
    return(Tvec)
  }
  
}


################################
######### Production ###########
################################

drawrandomnet <- function(){
  W <- matrix(rPareto(n*n,0.01,1),n,n) 
  normW <- rowSums(W)
  normW[normW==0] <- 1
  W <- diag(1/normW)%*%W
  
  W <- (1-lshare)*W
  return(W)
}




## Computes production levels assuming Cobb-Douglas production function for a fixed production network

production <- function(rstocks, parprod){
  rstocks[rstocks<exth] <- 0
  mu <- parprod[[1]] # distortions/taxes
  a <- parprod[[2]] # TFP
  B <- parprod[[3]] # network of ecological services
  W <- parprod[[4]] # production network
  lbd <- parprod[[5]] # share of taxes rebated to consumer
  q <- parprod[[6]] # final demand weights
  n <- nrow(W) # sectors
  m <- ncol(B) # species
  spcatch <- parprod[[7]]
  DD <- rep(0,min(n,m))
  DD[spcatch] <- Bpar # sector k fishes specie k
  DD <- diag(DD)
  B[1:ncol(DD),1:ncol(DD)] <- DD
  
  dW <- diag(1/(1+mu*lbd))%*%W # dot W
  ddW <- dW + diag(mu*lbd/(1+mu*lbd))%*%matrix(1,n,n)%*%diag(q) # ddot W
  alpha <- (1+mu*lbd)/(1+mu)
  
  # log prices
  
  check_fail <- sapply(1:n, function(i) any(rstocks[B[i,]>0]==0)) # production sectors with production failure
  
  ress <- sapply(1:n, function(i) sum( B[i,B[i,]>0]*log(rstocks[B[i,]>0])))
  ress[check_fail] <- 0 # tmp value
  
  lnp <- -solve(diag(n)-W)%*%( matrix(log(a),n,1) + ress -matrix(log(1+mu),n,1) )
  lnp[check_fail] <- 0 # tmp value
  
  # log production
  lny <- log(c(solve(diag(n)-t(ddW))%*%matrix(q,n,1))) - c(lnp) - log(alpha)
  lny[check_fail] <- 0 # tmp value
  prd <- exp(lny)
  prd[check_fail] <- 1e-8
  
  return(prd)
}
  
## computes specie-specific catches

catches <- function(rstocks,parprod){
  n <- nrow(parprod[[3]])
  m <- ncol(parprod[[3]])
  spcatch <- parprod[[7]]
  
  C <- matrix(0,m,n)
  DD <- rep(0,min(n,m))
  DD[spcatch] <- 1
  DD <- diag(DD)
  C[1:ncol(DD),1:ncol(DD)] <- DD
  
  
  catch <- c(C%*%matrix(production(rstocks,parprod),n,1))
  
  catch[rstocks<=exth] <- 0
  catch <- pmin(catch,rstocks)
  return(catch)
}

#########################################
######### Population dynamics ###########
#########################################

popdyn <- function(t,Bt,parameters){
  K <- parameters[[1]] #Carrying capacities
  r <- parameters[[2]] # max growth rates
  a <- parameters[[3]] # # ar, ax (invertebrate), ax (vertebrate)
  ylev <- parameters[[4]] # y(invertebrate), y(vertebrate)
  w <- parameters[[5]] # relative consumption rates
  b0 <- parameters[[6]] # half-saturation
  h <- parameters[[7]] # Hill coefficient
  c <- parameters[[8]] # Predator interference
  etmp <- parameters[[9]] # assimilation efficiency
  Z <- parameters[[10]] # body size ratio
  net <- parameters[[11]] # trophic network: i eats j
  vert <- parameters[[12]] # vertebrates binary
  parprod <- parameters[[13]] # parameters for the production/catch function
  adj <- parprod[[8]] # production adjusts to stocks
  
  Bt[Bt<exth] <- 0 # extinction threshold
  Bt[Bt>1e20] <- 1e20 # upper bound
  
  n <- length(K)
  basal <- (rowSums(net)==0) # id for basal species
  nbasal <- (rowSums(net)>0) # id for non-basal species
  
  ## 
  exists <- matrix(1,n,n)
  exists[,Bt==0] <- 0
  exists[Bt==0,] <- 0
  
  w <- rep(0,n)
  w[nbasal] <- 1/rowSums(exists*net)[nbasal]
  w[w==Inf] <- 0
  w <- matrix(rep(w,n),n,n)*net
  
  Fmatn <- w*matrix(rep(Bt,n),n,n,byrow=T)^h
  Fmatd1 <- b0^h + b0^h*matrix(rep(c*Bt,n),n,n)
  Fmatd2 <- matrix(rep(c(rowSums(Fmatn)),n),n,n) # denominator
  Fmat <- Fmatn/(Fmatd1 + Fmatd2)
  Fmat[is.nan(Fmat)] <- 0
  Tlev <- trolev(net) # trophic levels

  e <- matrix(0,n,n) # assimilation efficiency
  e[,nbasal] <- etmp[1] # for non-basal
  e[,basal] <- etmp[2] # for basal
  
  y <- rep(0,nspecies) # max consumption rate relative to metabolic rate
  y[nbasal] <- ylev[1]
  y[vert==1] <- ylev[2]
  ax <- rep(0,n)
  ax[nbasal] <- a[2]
  ax[vert==1] <- a[3]
  
  x <- ax*(Z^(Tlev-1))^(-0.25) # mass-specific metabolic rate
  
  gain <- matrix(rep(x*y*Bt,n),n,n)*Fmat
  loss <- gain/e
  
  ch <- rep(0,n) # change (first-difference)
  ch[basal] <- r[basal]*Bt[basal]*(1-sum(Bt[basal])/K[basal])  # logistic grow rate for basal species
  ch[nbasal] <- -x[nbasal]*Bt[nbasal] # maintenance loss (non-basal)
  ch <- ch - c(colSums(loss)) # loss to predation (all)
  ch <- ch + c(rowSums(gain)) # gains from predation
  if (extraction){
    if (adjust){
      ch <- ch - catches(Bt,parprod)
    } else {
      ch <- ch - catches(adj,parprod)
    }
  }
  ch[Bt==0] <- 0 # the biomass of extinct species does not vary
  
  return(list(ch))
}



