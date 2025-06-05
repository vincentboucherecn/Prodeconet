########################
######
####  This code computes complementary variables (e.g. volatility)
######
########################



library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(stringr)
library(readr)
rm(list=ls())

lstfiles <- list.files("./data/")
numbs <- sort(parse_number(lstfiles[grep("net",lstfiles)]))
alist <- readRDS("alistcal.Rda")
alist0 <- seq(from=0.01,to=0.105, by=0.005)/4

simu <- str_detect(lstfiles, paste0("simu_", collapse = '|'))
simu <- lstfiles[simu]

logp1 <- function(inp){
  return(log(inp+1))
}

ndata <- matrix(NA,length(simu)*30,27)
pos <- 1
for (i in numbs){
  for (j in 1:length(alist)){
    tryCatch({
      s <- readRDS(paste("data/simu_",i,"_",alist[j],".Rda",sep=""))
      s0 <- readRDS(paste("data/simu0_",i,"_",alist0[j],".Rda",sep=""))
      c <- readRDS(paste("data/catch_",i,"_",alist[j],".Rda",sep=""))
      c0 <- readRDS(paste("data/catch0_",i,"_",alist0[j],".Rda",sep=""))
      sd.1 <- (sapply(2:31,function(r) sd(logp1(s[2:500,r])-logp1(s[1:499,r]))))
      sd.1[is.na(sd.1)] <- 0
      sd.2 <- (sapply(2:31,function(r) sd(logp1(s[501:1000,r])-logp1(s[500:999,r]))))
      sd.2[is.na(sd.2)] <- 0
      sd.3 <- (sapply(2:31,function(r) sd(logp1(s[1001:1500,r])-logp1(s[1000:1499,r]))))
      sd.3[is.na(sd.3)] <- 0
      sd.4 <- (sapply(2:31,function(r) sd(logp1(s[1501:2000,r])-logp1(s[1500:1999,r]))))
      sd.4[is.na(sd.4)] <- 0
      
      sd0.1 <- (sapply(2:31,function(r) sd(logp1(s0[2:500,r])-logp1(s0[1:499,r]))))
      sd0.1[is.na(sd0.1)] <- 0
      sd0.2 <- (sapply(2:31,function(r) sd(logp1(s0[501:1000,r])-logp1(s0[500:999,r]))))
      sd0.2[is.na(sd0.2)] <- 0
      sd0.3 <- (sapply(2:31,function(r) sd(logp1(s0[1001:1500,r])-logp1(s0[1000:1499,r]))))
      sd0.3[is.na(sd0.3)] <- 0
      sd0.4 <- (sapply(2:31,function(r) sd(logp1(s0[1501:2000,r])-logp1(s0[1500:1999,r]))))
      sd0.4[is.na(sd0.4)] <- 0
      
      cd.1 <- (sapply(1:30,function(r) sd(logp1(c[2:500,r])-logp1(c[1:499,r]))))
      cd.1[is.na(cd.1)] <- 0
      cd.2 <- (sapply(1:30,function(r) sd(logp1(c[501:1000,r])-logp1(c[500:999,r]))))
      cd.2[is.na(cd.2)] <- 0
      cd.3 <- (sapply(1:30,function(r) sd(logp1(c[1001:1500,r])-logp1(c[1000:1499,r]))))
      cd.3[is.na(cd.3)] <- 0
      cd.4 <- (sapply(1:30,function(r) sd(logp1(c[1501:2000,r])-logp1(c[1500:1999,r]))))
      cd.4[is.na(cd.4)] <- 0
      
      cd0.1 <- (sapply(1:30,function(r) sd(logp1(c0[2:500,r])-logp1(c0[1:499,r]))))
      cd0.1[is.na(cd0.1)] <- 0
      cd0.2 <- (sapply(1:30,function(r) sd(logp1(c0[501:1000,r])-logp1(c0[500:999,r]))))
      cd0.2[is.na(cd0.2)] <- 0
      cd0.3 <- (sapply(1:30,function(r) sd(logp1(c0[1001:1500,r])-logp1(c0[1000:1499,r]))))
      cd0.3[is.na(cd0.3)] <- 0
      cd0.4 <- (sapply(1:30,function(r) sd(logp1(c0[1501:2000,r])-logp1(c0[1500:1999,r]))))
      cd0.4[is.na(cd0.4)] <- 0
      
      ndata[pos:(pos+29),1] <- i
      ndata[pos:(pos+29),2] <- alist[j]
      ndata[pos:(pos+29),3] <- sd.1
      ndata[pos:(pos+29),4] <- sd.2
      ndata[pos:(pos+29),5] <- sd.3
      ndata[pos:(pos+29),6] <- sd.4
      
      ndata[pos:(pos+29),7] <- alist0[j]
      ndata[pos:(pos+29),8] <- sd0.1
      ndata[pos:(pos+29),9] <- sd0.2
      ndata[pos:(pos+29),10] <- sd0.3
      ndata[pos:(pos+29),11] <- sd0.4
      
      ndata[pos:(pos+29),12] <- s0[1,2:31]
      ndata[pos:(pos+29),13] <- s[1,2:31]
      ndata[pos:(pos+29),14] <- c0[1,]
      ndata[pos:(pos+29),15] <- c[1,]
      
      ndata[pos:(pos+29),16] <- cd.1
      ndata[pos:(pos+29),17] <- cd.2
      ndata[pos:(pos+29),18] <- cd.3
      ndata[pos:(pos+29),19] <- cd.4
      ndata[pos:(pos+29),20] <- cd0.1
      ndata[pos:(pos+29),21] <- cd0.2
      ndata[pos:(pos+29),22] <- cd0.3
      ndata[pos:(pos+29),23] <- cd0.4
      
      ndata[pos:(pos+29),24] <- colSums(s0[1:500,2:31])
      ndata[pos:(pos+29),25] <- colSums(s[1:500,2:31]) 
      ndata[pos:(pos+29),26] <- colSums(c0[1:500,1:30])
      ndata[pos:(pos+29),27] <- colSums(c[1:500,1:30])
      
    }, error=function(e){})
    pos <- pos + 30
  }
  print(i)
}
ndata <- as.data.frame(ndata)
colnames(ndata) <- c("net","a","sd.1","sd.2","sd.3","sd.4","a0","sd0.1","sd0.2","sd0.3","sd0.4","s01","s1","c01","c1",
                     "cd.1","cd.2","cd.3","cd.4","cd0.1","cd0.2","cd0.3","cd0.4","avs0","avs","avc0","avc")
ndata$spid <- rep(1:30,length(simu))


pdata1 <- ndata[,c(1:6,13,15,16:19,25,27,28)]
colnames(pdata1) <- c("net","a","sd1","sd2","sd3","sd4","s","c","cd1","cd2","cd3","cd4","avs",'avc',"spid")
pdata1$type <- "1"

pdata2 <- ndata[,c(1,7:12,14,20:23,24,26,28 )]
colnames(pdata2) <- c("net","a","sd1","sd2","sd3","sd4","s","c","cd1","cd2","cd3","cd4","avs",'avc',"spid")
pdata2$type <- "0"


pdata <- rbind(pdata1,pdata2)
pdata$type <- factor(pdata$type)
pdata$qc <- ntile(pdata$c,n=15)

pdata$exists <- pdata$s>0
pdata$catch <- pdata$c>0

saveRDS(pdata,"pdata.Rda")

save.image("spagplot.RData")


