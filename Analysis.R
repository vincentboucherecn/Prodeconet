########################
######
####  This code performs the analysis in Section 6
######
########################



library(ggplot2)
library(lmtest)
library(sandwich)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(gridExtra)

rm(list=ls())
dta <- readRDS("pdata.Rda")
lstfiles <- list.files("./data/")
numbs <- sort(parse_number(lstfiles[grep("net",lstfiles)]))
alist <- readRDS("alistcal.Rda")
alist0 <- seq(from=0.01,to=0.105, by=0.005)/4

###
# Creates the dataset for caught species
###

dta$qc <- NULL
dtac <- dta[dta$c>0,]
dtac$cs <- dtac$c/dtac$s

dtac$qcs_sp <- NA
for (sp in unique(dtac$spid)){
  dtac[dtac$spid==sp,"qcs_sp"] <- ntile(dtac[dtac$spid==sp,"cs"],n=20)
}


saveRDS(dtac,"dtac.Rda")

###
# Additional functions
###

                                         
prcatch <- function(m,i,tmax){
  #m <- which.min(mcheck[,3])
  aa0 <- readRDS(paste("data/catch0_",m,"_",alist0[i],".Rda",sep=""))
  aa <- readRDS(paste("data/catch_",m,"_",alist[i],".Rda",sep=""))
  aa0 <- as.data.frame(aa0)
  aa0$time <- 1:nrow(aa0)
  aa <- as.data.frame(aa)
  aa$time <- 1:nrow(aa)
  long0 <- aa0 %>%
    pivot_longer(
      cols = `V1`:`V30`,
      names_to = "specie",
      values_to = "value"
    )
  
  long <- aa %>%
    pivot_longer(
      cols = `V1`:`V30`,
      names_to = "specie",
      values_to = "value"
    )
  l <- by(long, long$specie, function(x)
    within(x, valnorm <- 100*(value-value[1])/value[1]))
  dat <- do.call(rbind, l)
  l0 <- by(long0, long0$specie, function(x)
    within(x, valnorm <- 100*(value-value[1])/value[1]))
  dat0 <- do.call(rbind, l0)
  dat$valnorm[is.nan(dat$valnorm)] <- 0
  dat0$valnorm[is.nan(dat0$valnorm)] <- 0
  p0 <- ggplot(dat0[dat0$time<=tmax,], aes(x=time,y=valnorm,color = factor(specie))) +
    geom_line(show.legend = FALSE) +
    theme_minimal() + ylab("Catches variation (%)") + xlab("Time") + coord_cartesian(ylim=c(-100,20)) +
    theme(text = element_text(size = 26), plot.title = element_text(size=24,color="#F7931D",hjust = 0.5)) + ggtitle("Isolated Sectors")
  p <- ggplot(dat[dat$time<=tmax,], aes(x=time,y=valnorm,color = factor(specie))) +
    geom_line(show.legend = FALSE) +
    theme_minimal() + ylab("catches variation (%)") + xlab("Time") + coord_cartesian(ylim=c(-100,20)) +
    theme(text = element_text(size = 26), plot.title = element_text(size=24,color="#213F99",hjust = 0.5)) + ggtitle("Production Network")
  
  grid.arrange(p0, p, ncol=2)
  pdf(paste("example_catches",i,".pdf",sep=""), width = 12, height = 9) # Open a new pdf file
  grid.arrange(p0, p, ncol=2) # Write the grid.arrange in the file
  dev.off()
}

prsimu <- function(m,i,tmax){
  #m <- which.min(mcheck[,3])
  aa0 <- readRDS(paste("data/simu0_",m,"_",alist0[i],".Rda",sep=""))
  aa <- readRDS(paste("data/simu_",m,"_",alist[i],".Rda",sep=""))
  aa0 <- as.data.frame(aa0)
  aa <- as.data.frame(aa)
  long0 <- aa0 %>%
    pivot_longer(
      cols = `1`:`30`,
      names_to = "specie",
      values_to = "value"
    )
  
  long <- aa %>%
    pivot_longer(
      cols = `1`:`30`,
      names_to = "specie",
      values_to = "value"
    )
  l <- by(long, long$specie, function(x)
    within(x, valnorm <- 100*(value-value[1])/value[1]))
  dat <- do.call(rbind, l)
  l0 <- by(long0, long0$specie, function(x)
    within(x, valnorm <- 100*(value-value[1])/value[1]))
  dat0 <- do.call(rbind, l0)
  dat$valnorm[is.nan(dat$valnorm)] <- 0
  dat0$valnorm[is.nan(dat0$valnorm)] <- 0
  p0 <- ggplot(dat0[dat0$time<=tmax,], aes(x=time,y=valnorm,color = factor(specie))) +
    geom_line(show.legend = FALSE) +
    theme_minimal() + ylab("Biomass variation (%)") + xlab("Time") + coord_cartesian(ylim=c(-100,100)) +
    theme(text = element_text(size = 26), plot.title = element_text(size=24,color="#F7931D",hjust = 0.5)) + ggtitle("Isolated Sectors")
  p <- ggplot(dat[dat$time<=tmax,], aes(x=time,y=valnorm,color = factor(specie))) +
    geom_line(show.legend = FALSE) +
    theme_minimal() + ylab("Biomass variation (%)") + xlab("Time") + coord_cartesian(ylim=c(-100,100)) +
    theme(text = element_text(size = 26), plot.title = element_text(size=24,color="#213F99",hjust = 0.5)) + ggtitle("Production Network")
  
  grid.arrange(p0, p, ncol=2)
  pdf(paste("example_biomass",i,".pdf",sep=""), width = 12, height = 9) # Open a new pdf file
  grid.arrange(p0, p, ncol=2) # Write the grid.arrange in the file
  dev.off()
}
alev <- 15
alist <- readRDS("alistcal.Rda")
alist0 <- seq(from=0.01,to=0.105, by=0.005)/4

##
# select example
##
## ** Note ** selects the ecosystem that is such that the initial catches are the closest for the the 15th level of TFP (High TFP)

wa <- c(alist0[alev],alist[alev])
dtaselect <- dtac[,c("net","a","c","type","spid")]
dtaselect <- dtaselect[dtaselect$a %in% wa,]
dtaselect2 <- aggregate(dtaselect$c, by=list(dtaselect$net,dtaselect$a,dtaselect$type),FUN=function(i)mean(i,na.rm = T))
nets <- unique(dtaselect2$Group.1)
dtaselect3 <- rep(NA,length(nets))
for (i in 1:length(nets)){
  dtaselect3[i] <-abs(dtaselect2$x[which(dtaselect2$Group.1==nets[i] & dtaselect2$Group.3==0)]-dtaselect2$x[which(dtaselect2$Group.1==nets[i] & dtaselect2$Group.3==1)])
}
idx <- sort(dtaselect3,index.return=T)
prnet <- nets[idx$ix[1]]


###
# Spag. plots
###
prcatch(prnet,alev,500)

prsimu(prnet,alev,500)

alev <- 8
prcatch(prnet,alev,500)

prsimu(prnet,alev,500)


###
# Specie-level graphs
###

dtac <- dtac[dtac$qcs_sp<15,] # print only first initial catch levels (prevents ecosystem collapse)
dtac$qcs_sp <- factor(dtac$qcs_sp)
dtac$Type <- recode(dtac$type,"0"="Isolated Sectors","1"="Production Network")


###
# Creates the dataset for ecosystem level analysis
###

dta_ag <- aggregate(dta,by=list(dta$net,dta$a,dta$type),FUN=function(x) mean(x,na.rm=T))
dta_ag$qc <- factor(ntile(dta_ag$c,n=15))
dta_ag$qc_net <- NA
for (net in unique(dta_ag$net)){
  dta_ag[dta_ag$net==net,"qc_net"] <- ntile(dta_ag[dta_ag$net==net,"c"],n=20)
}
dta_ag$qc_net <- factor(dta_ag$qc_net)
dta_ag$al <- NULL
dta_ag$al[dta_ag$Group.3=="0"] <- sapply(dta_ag$a[dta_ag$Group.3=="0"],function(i) which(alist0==i))
dta_ag$al[dta_ag$Group.3=="1"] <- sapply(dta_ag$a[dta_ag$Group.3=="1"],function(i) which(alist==i))
dta_ag$al <- factor(dta_ag$al)


p1 <- ggplot(dtac, aes(x=qcs_sp, y=cd1, fill=Type)) +
  geom_boxplot(outliers=F) +
  #    stat_summary(geom = "point", fun = \(x) quantile(x, 0.95,na.rm=T),shape=16, size=4,aes(color=type)) +
  
  theme_minimal() + scale_fill_manual(values=c("#F7931D", "#213F99")) +
  ylab("Volatility of catches") + xlab("Initial catches") + theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())+ theme(legend.position = c(0.65, .95),
                                         legend.justification = c("left", "top"),
                                         legend.box.just = "left",legend.title=element_blank(),text = element_text(size = 26))
plot(p1)
ggsave("volatility_catches.pdf", width = 12, height = 9)

p3 <- ggplot(dtac, aes(x=qcs_sp, y=avc, fill=Type)) +
  geom_boxplot(outliers=F) +
  #    stat_summary(geom = "point", fun = \(x) quantile(x, 0.95,na.rm=T),shape=16, size=4,aes(color=type)) +
  
  theme_minimal() + scale_fill_manual(values=c("#F7931D", "#213F99")) +
  ylab("Total catches") + xlab("Initial catches") + theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())+ theme(legend.position = c(0.65, .95),
                                         legend.justification = c("left", "top"),
                                         legend.box.just = "left",legend.title=element_blank(),text = element_text(size = 26))
plot(p3)
ggsave("average_catches.pdf", width = 12, height = 9)


p2 <- ggplot(dtac, aes(x=qcs_sp, y=sd1, fill=Type)) +
  geom_boxplot(outliers=F) +
  #    stat_summary(geom = "point", fun = \(x) quantile(x, 0.95,na.rm=T),shape=16, size=4,aes(color=type)) +
  
  theme_minimal() + scale_fill_manual(values=c("#F7931D", "#213F99")) +
  ylab("Volatility of biomass (caught species)") + xlab("Initial catches") + theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())+ theme(legend.position = c(0.65, .95),
                                         legend.justification = c("left", "top"),
                                         legend.box.just = "left",legend.title=element_blank(),text = element_text(size = 26))
plot(p2)
ggsave("volatility_biomass.pdf", width = 12, height = 9)

p4 <- ggplot(dtac, aes(x=qcs_sp, y=avs, fill=Type)) +
  geom_boxplot(outliers=F) +
  #    stat_summary(geom = "point", fun = \(x) quantile(x, 0.95,na.rm=T),shape=16, size=4,aes(color=type)) +
  
  theme_minimal() + scale_fill_manual(values=c("#F7931D", "#213F99")) +
  ylab("Total biomass (caught species)") + xlab("Initial catches") + theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())+ theme(legend.position = c(0.65, .95),
                                         legend.justification = c("left", "top"),
                                         legend.box.just = "left",legend.title=element_blank(),text = element_text(size = 26))
plot(p4)
ggsave("average_biomass.pdf", width = 12, height = 9)



###
# Ecosystem level analysis
###


out1 <- (lm(I(sd1/sd(dta_ag$sd1)) ~ I(c/sd(dta_ag$c))*Group.3 + factor(net),  data=dta_ag))
cout1 <- coeftest(out1, vcov = vcovHC(out1, type = 'HC0'))
print(cout1)

out2 <- (lm(I(cd1/sd(dta_ag$cd1)) ~ I(c/sd(dta_ag$c))*Group.3 + factor(net),  data=dta_ag))
cout2 <- coeftest(out2, vcov = vcovHC(out2, type = 'HC0'))
print(cout2)

out3 <- (lm(I(sd1/sd(dta_ag$sd1)) ~ I(c/sd(dta_ag$c)) + Group.3 + factor(net),  data=dta_ag))
cout3 <- coeftest(out3, vcov = vcovHC(out3, type = 'HC0'))
print(cout3)

out4 <- (lm(I(cd1/sd(dta_ag$cd1)) ~ I(c/sd(dta_ag$c)) + Group.3 + factor(net),  data=dta_ag))
cout4 <- coeftest(out4, vcov = vcovHC(out4, type = 'HC0'))
print(cout4)

