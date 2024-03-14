########################
######
####  (Information only) This code tabulates network rejection reasons
######
########################

rm(list=ls())

#####################
## gather available networks/steady-states (accepted and rejected)
#####################

lstfiles <- list.files("./data/")
numbs <- sort(parse_number(lstfiles[grep("reason",lstfiles)]))
nok <- sort(parse_number(lstfiles[grep("binit",lstfiles)]))

#####################
## Format reason data
#####################

reasons <- data.frame(reason=NULL,net=NULL)
for (i in numbs){
  a <- read.table(paste("data/reason",i,".txt",sep=""))
  a$net <- i
  colnames(a)[1] <- "reason"
  reasons <- rbind(reasons,a)
}
reasons0 <- reasons[reasons$reason!="",]

#####################
## Print some summary statistics
#####################

print(table(reasons0$reason))
print("Total:")
print((length(numbs)+length(nok)))
print("Wrong:")
print(length(numbs)/(length(numbs)+length(nok)))
reasons00 <- reasons0[reasons0$reason!="No steady-state found",]
print("Wrong0:")
print(length(unique(reasons00$net))/(length(unique(reasons00$net))+length(nok)))

