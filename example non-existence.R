########################
######
####  This code produces Figure C.1
######
########################


rm(list=ls())
set.seed(9876)
library(ggplot2)

#########
### Parameters
#########

alpha <- 0.7
phi <- 0.5
K <- 2/3
r <- 1
e <- T
epsilon <- 0.01

#########
### Difference equation
#########

tomorrow <- function(x,e){
  if (e){
    xt <- max(x + r*x*(1-x/K)-as.numeric(x>=(1-alpha)/phi)*phi*x,epsilon)
  } else {
    xt <- max(x + r*x*(1-x/K),epsilon)
  }
  
  return(xt)
}

#########
### Compute dynamics and save graphs
#########

dom <- seq(from=0,to=1,by=0.01)
ima <- sapply(dom,function(w) tomorrow(w,e))

library(ggplot2)
theme_set(theme_minimal())
dta <- data.frame(dom=dom,ima=ima)
p <- ggplot(dta, aes(x=dom,y=ima)) + geom_point() + geom_abline(slope=1)
ggsave("failure.pdf",width = 12, height = 9)
plot(p)

t <- 1:100
x <- rep(0.01,length(t))
for (tt in 2:length(x)){
  x[tt] <- tomorrow(x[(tt-1)],e)
}
dtat <- data.frame(t=t,x=x)

theme_set(theme_minimal())

p <- ggplot(dtat, aes(x=t,y=x)) + geom_line()
ggsave("failure2.pdf",width = 12, height = 9)
plot(p)



