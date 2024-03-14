########################
######
####  This code produces Figure 5
######
########################


library(combinat)
rm(list=ls())
set.seed(987654321)
library(ggplot2)

################################
######### PARAMETERS ###########
################################

N <- 2 # number of sectors
w <- c(1,1) # vector of consumption weights
w <- w/sum(w)
alphapotential <- matrix(c(0,0,1,0),N,N)
laborweithgs <- c(3,1)

mu <- rep(0,N) # vector of taxation
M <- 1 # number of resources

################################
######### THE ECONOMY ##########
################################

# List input combinations
comb <- unlist(lapply(1:N,combn, x = 1:N, simplify = FALSE), recursive = FALSE) # all possible subsets
mperm <- matrix(0,(length(comb)+1),N) # last possibility is only labor (adds empty set as possible subset)
nc <- nrow(mperm) # number of input combination choices
## write as a binary matrix
for (i in 1:(nc-1)){
  mperm[i,] <- as.numeric(1:N %in% comb[[i]])
}

### TFP for each firm, for each input combination
tech <- matrix(0,nc,N)
tech[2,] <- c(10,0) # using only firm 2 as inputs
tech[4,] <- 30 # using only labor


### dependence list
## here sector 2 is substitute to resource 1
blist <- vector("list",nc)
bmat <- matrix(0,N,M) # economic dependence on resources
blist[[1]] <- bmat # no dependence TFP=0
blist[[3]] <- bmat # no dependence TFP=0
blist[[2]] <- bmat # no dependence when using firm 2 as input
bmat[1,1] <- 0.5 # dependence on stock for firm 1 
blist[[4]] <- bmat # dependence if using only labor
 

### externatlity list
tlist <- vector("list",nc)
tmat <- matrix(0,M,N) # externalities
tlist[[1]] <- tmat
tlist[[2]] <- tmat
tlist[[3]] <- tmat
tmat[1,1] <- -1 # baseline externality when using the resource
tlist[[4]] <- tmat


#############################
##### find equilibrium ######
#############################

lsolve <- function(ich){
  
  ## finds prices, production, tfp and unitary cost for fixed input choices
  
  tlda <- amat <- matrix(0,N,N) # initialize matrices of elasticities
  tfp <- matrix(0,N,1) # initialize vector of TFP
  bmat <- matrix(0,N,M)
  for (i in 1:N){
    bmat[i,] <- blist[[ich[i]]][i,]
    winput <- which(mperm[ich[i],]>0)
    if (length(winput)>0){
      amat[i,] <- alphapotential[i,winput]/(sum(alphapotential[i,winput])+laborweithgs[i]) # alpha
    }
    tlda[i,] <- (amat[i,]+mu[i]*w[i])/(1+mu[i]) # tilde alpha
    tfp[i,1] <- tech[ich[i],i] # tfp
  }
  lp <- -solve(diag(N)-amat)%*%(log(tfp)-matrix(log(1+mu),N,1))-solve(diag(N)-amat)%*%bmat%*%matrix(lX,M,1) # equilibrium prices for fixed input choice
  ly <- log(solve(diag(N)-t(tlda))%*%matrix(w,N,1))-lp # equilibrium production level (for fixed input choice)
  
  lucost <- (1+mu)*exp(amat%*%lp)/exp(log(tfp)+bmat%*%matrix(lX,M,1)) # equilibrium unitary cost for fixed input choice
  return(cbind(ly,lp,tfp,lucost))
}

fprices <- function(i,ichi,prices){

  ### unitary cost function for fixed prices
  bmat <- matrix(0,N,M)
  bmat[i,] <- blist[[ichi]][i,]
  winput <- which(mperm[ichi,]>0)
  if (length(winput)>0){
    avec <- matrix(0,1,N)
    avec[1,winput] <- alphapotential[i,winput]/(sum(alphapotential[i,winput])+laborweithgs[i])
  } else {
    avec <- matrix(0,1,N)
  }
  tfpi <- tech[ichi,i] # tfp for sector i (fixed input choices)
  if (tfpi==0){
    lucost=Inf
  } else {
    lucost <- (1+mu[i])*exp(avec%*%log(prices))/exp(log(tfpi)+c(bmat[i,]%*%matrix(lX,M,1))) # unitary production cost for fixed input choice and prices
  }
  return(c(lucost))
}

findichoices <- function(){
  
  ## optimal choice of input
  
  p0 <- rep(1e-8,N) # initial price vector (very small)
  st <- 0
  while (st==0){ # loop until equilibrium (converges to the least eq by Tarski)
    p1 <- p0
    for (i in 1:N){ # for each sector
      cst <- sapply(1:nc,function(ch) fprices(i,ch,p0)) # unitary cost for each input choices at current prices
      p1[i] <- min(cst) # minimal unitary cost
    }
    
    if (sum((p1-p0)^2)<1e-12){ # convergence criterion
      st <- 1
    }
    p0 <- p1
  }
  ich <- rep(1,N)
  for (i in 1:N){ # for each sector
    cst <- sapply(1:nc,function(ch) fprices(i,ch,p1)) # unitary cost for each input choices at equilibrium prices
    ich[i] <- which.min(cst)[1] # optimal input choices
  }
  return(ich)
}



################################
####### THE ENVIRONMENT ########
################################

gmat <- matrix(0,M,M) # interaction matrix
diag(gmat) <- -0.2

r <- matrix(c(5),M,1) # growth rates
lX <- solve(diag(M)-gmat)%*%r # start at steady state with the economy
#lX <- -100

tmax <- 100 # simulate tmax periods
Xmat <- matrix(NA,M,tmax) # matrix of stock of resources
Ymat <- matrix(NA,N,tmax) # matrix of productions
Imat <- matrix(NA,N,tmax) # matrix of input choices
Cmat <- matrix(NA,N,tmax) # matrix of unit costs
for (t in 1:tmax){
  Xmat[,t] <- exp(lX) # store current stock of resources
  Imat[,t] <- findichoices() # compute optimal input choices
  tmat <- matrix(0,M,N)
  for (i in 1:N){
    tmat[,i] <- tlist[[Imat[i,t]]][,i]
  }
  Ymat[,t] <- exp(lsolve(Imat[,t])[,1]) # compute production at optimal input choices
  Cmat[,t] <- (lsolve(Imat[,t])[,4]) # compute production at optimal input choices
  tlX <- lX # temporary vector
  for (i in 1:M){ # for each resource
    tlX[i] <- r[i]+sum(c(gmat[i,])*lX) + sum(c(tmat[i,])*c(log(Ymat[,t]))) # update according to GLV with externalities
  }
  lX <- tlX
}

print(Imat)

dta <- as.data.frame(cbind(c(Ymat[1,1:10]),c(Xmat[1,1:10]),"Fish production",1:10),stringsAsFactors = F)
dta <- rbind(dta,as.data.frame(cbind(c(Ymat[2,1:10]),c(Xmat[1,1:10]),"Manufacturing",1:10),stringsAsFactors = F))
colnames(dta) <- c("production","stock","Sector","t")
dta$production <- as.numeric(dta$production)
dta$stock <- as.numeric(dta$stock)
dta$t <- as.numeric(dta$t)
dta$Sector <- factor(dta$Sector)
dta$Resource <- factor("Wild fish")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colorswendy <- c("#F7931D","#213F99")


p_prod <- ggplot(data = dta) + aes(x = t, y = production, shape=Sector) + 
  geom_line(linetype = "dashed",color=colorswendy[1]) + 
  geom_point(aes(shape=Sector),size=3,color=colorswendy[1]) + 
  ylab(expression("Production")) +
  scale_shape_manual(values=c(15,17)) +
  scale_x_discrete(name ="Time",limits=1:10) + theme_minimal() + 
  theme(legend.position="top", axis.text=element_text(size=11),
        axis.title=element_text(size=12),aspect.ratio=3/4) + scale_colour_manual(values=cbPalette)

p_stock <- ggplot(data = dta[dta$Sector=="Fish production",]) + aes(x = t, y = stock) + 
  geom_line(linetype = "dashed",color=colorswendy[2]) + 
  geom_point(shape=16,size=3,color=colorswendy[2]) + 
  ylab(expression("Resource Abundance")) + xlab(expression("t")) +
  scale_x_discrete(name ="Time",limits=1:10) + theme_minimal()  +  
  theme(legend.position="top", axis.text=element_text(size=11),
        axis.title=element_text(size=12),aspect.ratio=3/4) + scale_colour_manual(values=cbPalette)

grid.arrange(p_prod, p_stock, ncol=2)
