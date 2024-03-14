########################
######
####  This code produces Figure 6 and Figure F.1 
######
########################

library(ggplot2)
library(dplyr)

########################
## Load, format, and save data for the simulations with production network
########################

rm(list=ls())
load("tot.RData")
dbtot <- data.frame(A=NA,extinct=NA,abs1=NA,abs=NA,stable1=NA,stable=NA,catches=NA)
for (j in 1:i){
  db <- Wbase[[j]]
  stab <- aggregate(db$abs, by=list(db$nnet), max)
  ro <- c(alist[j],mean(db$Btot<exth,na.rm = T), mean(db$abs<1, na.rm=T),mean(db$abs,na.rm=T),mean(stab$x<1),mean(stab$x),mean(db$catches[db$catches0>0],na.rm=T))
  dbtot <- rbind(dbtot,ro)
}

saveRDS(dbtot,"netdta.Rda")

########################
## Load, format, and save data for the simulations without production network
########################

rm(list=ls())
load("tot0.RData")

dbtot0 <- data.frame(A=NA,extinct=NA,abs1=NA,abs=NA,stable1=NA,stable=NA,catches=NA)
for (j in 1:i){
  db <- W0[[j]]
  db[db=="Error in if (max(abs(unlist(popdyn(1, Btot, parameters)))) > 1e-08) { : \n  missing value where TRUE/FALSE needed\n"] <- NA
  db <- db %>% mutate_at(colnames(db), as.numeric)
  stab <- aggregate(db$abs, by=list(db$nnet), max)
  ro <- c(alist[j],mean(db$Btot<exth,na.rm = T), mean(db$abs<1, na.rm=T),mean(db$abs,na.rm=T),mean(stab$x<1),mean(stab$x),mean(db$catches[db$catches0>0],na.rm=T))
  dbtot0 <- rbind(dbtot0,ro)
}

saveRDS(dbtot0,"netdta0.Rda")

########################
## Load, format, print and save graphs
########################

rm(list=ls())
dtaW <- readRDS("netdta.Rda")
dtaW0 <- readRDS("netdta0.Rda")
dtaW$type <- "Production Network"
dtaW0$type <- "Isolated sectors"

dtaW <- dtaW[complete.cases(dtaW),]
dtaW0 <- dtaW0[complete.cases(dtaW0),]
dtaW$Ao <- order(dtaW$A)
dtaW0$Ao <- order(dtaW0$A)
dtaW <- rbind(dtaW,dtaW0)
dtaW$type <- factor(dtaW$type)
dtaW$Ao <- as.numeric(dtaW$Ao)

#p <- ggplot(dtaW, aes(catches, stable, color = type)) +
#  geom_point(shape = 16, size = 3, show.legend = FALSE) +
#  theme_minimal() + geom_smooth(method = "lm",fill="lightgrey",
#  formula = y ~ poly(x, 5)) + scale_color_manual(values=c("#F7931D", "#213F99")) +
#  ylab("Spectral radius (average)") + xlab("Seafood catches (average)") +
#  theme(
#    legend.position = c(0.05, .95),
#    legend.justification = c("left", "top"),
#    legend.box.just = "left",
#    text = element_text(size = 26)
#  ) + guides(color=guide_legend(title=NULL))
#plot(p)
#ggsave("avSigma.pdf", width = 12, height = 9)

p <- ggplot(dtaW, aes(catches, stable1, color = type)) +
  geom_point(shape = 16, size = 3, show.legend = FALSE) +
  theme_minimal() + geom_smooth(method = "lm",fill="lightgrey",
                                formula = y ~ poly(x, 5)) + scale_color_manual(values=c("#F7931D", "#213F99")) +
  ylab("Fraction of stable networks") + xlab("Seafood catches (average)") +
  theme(
    legend.position = c(0.65, .95),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    text = element_text(size = 26)
  ) + guides(color=guide_legend(title=NULL))
plot(p)
ggsave("avStable.pdf", width = 12, height = 9)

#p <- ggplot(dtaW, aes(catches, abs, color = type)) +
#  geom_point(shape = 16, size = 3, show.legend = FALSE) +
#  theme_minimal() + geom_smooth(method = "lm",fill="lightgrey",
#                                formula = y ~ poly(x, 5)) + scale_color_manual(values=c("#F7931D", "#213F99")) +
#  ylab("Eigenvalues (average, absolute value)") + xlab("Seafood catches (average)") +
#  theme(
#    legend.position = c(0.05, .25),
#    legend.justification = c("left", "top"),
#    legend.box.just = "left",
#    text = element_text(size = 26)
#  ) + guides(color=guide_legend(title=NULL))
#plot(p)
#ggsave("aveigen.pdf", width = 12, height = 9)

#p <- ggplot(dtaW, aes(catches, abs1, color = type)) +
#  geom_point(shape = 16, size = 3, show.legend = FALSE) +
#  theme_minimal() + geom_smooth(method = "lm",fill="lightgrey",
#                                formula = y ~ poly(x, 5)) + scale_color_manual(values=c("#F7931D", "#213F99")) +
#  ylab("Fraction of eigenvalues in unit disk") + xlab("Seafood catches (average)") +
#  theme(
#    legend.position = c(0.05, .25),
#    legend.justification = c("left", "top"),
#    legend.box.just = "left",
#    text = element_text(size = 26)
#  ) + guides(color=guide_legend(title=NULL))
#plot(p)
#ggsave("aveigen1.pdf", width = 12, height = 9)

p <- ggplot(dtaW, aes(catches, extinct, color = type)) +
  geom_point(shape = 16, size = 3, show.legend = FALSE) +
  theme_minimal() + geom_smooth(method = "lm",fill="lightgrey",
                                formula = y ~ poly(x, 5)) + scale_color_manual(values=c("#F7931D", "#213F99")) +
  ylab("Fraction of extinct species (average)") + xlab("Seafood catches (average)") +
  theme(
    legend.position = c(0.05, .95),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    text = element_text(size = 26)
  ) + guides(color=guide_legend(title=NULL))
plot(p)
ggsave("avext.pdf", width = 12, height = 9)

p <- ggplot(dtaW, aes(Ao, extinct, color = type)) +
  geom_point(shape = 16, size = 3, show.legend = FALSE) +
  theme_minimal() + geom_smooth(method = "lm",fill="lightgrey",
                                formula = y ~ poly(x, 5)) + scale_color_manual(values=c("#F7931D", "#213F99")) +
  ylab("Fraction of extinct species (average)")  + xlab("TFP") +
  theme(
    legend.position = c(0.05, .95),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    text = element_text(size = 26),
    axis.text.x = element_blank()
  ) + guides(color=guide_legend(title=NULL))
plot(p)
ggsave("avextAo.pdf", width = 12, height = 9)


p <- ggplot(dtaW, aes(Ao, stable1, color = type)) +
  geom_point(shape = 16, size = 3, show.legend = FALSE) +
  theme_minimal() + geom_smooth(method = "lm",fill="lightgrey",
                                formula = y ~ poly(x, 5)) + scale_color_manual(values=c("#F7931D", "#213F99")) +
  ylab("Fraction of stable networks") + xlab("TFP") +
  theme(
    legend.position = c(0.65, .95),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    text = element_text(size = 26),
    axis.text.x = element_blank()
  ) + guides(color=guide_legend(title=NULL))
plot(p)
ggsave("avStableAo.pdf", width = 12, height = 9)
