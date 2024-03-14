########################
######
####  This code produces Figure 7 and Figures F.2, F.3 and F.4
####  FAO data can be downloaded here: https://www.fao.org/fishery/statistics-query/en/global_production/global_production_quantity
######
########################


library(readr)
library(tidyverse)
library(rfishbase)
library(ggplot2)

rm(list=ls())
theme_set(theme_minimal())

######
## For World, US, Chile and China, produce figures
######


idlist <- c("World","US","Chile","China")
for (id in idlist){
  
  ### FAO data from two sources.
  
  production <- read_csv("Global_production_quantity.csv")
  production[production$PRODUCTION_SOURCE_DET.CODE!="CAPTURE","PRODUCTION_SOURCE_DET.CODE"] <- "Aquaculture"
  production[production$PRODUCTION_SOURCE_DET.CODE=="CAPTURE","PRODUCTION_SOURCE_DET.CODE"] <- "Capture"
  
  if (id=="US"){
    production <- production[production$COUNTRY.UN_CODE=="840",]  ## Keep US
  }
  if (id=="Chile"){
    production <- production[production$COUNTRY.UN_CODE=="152",]  ## Keep Chile only
  }
  if (id=="China"){
    production <- production[production$COUNTRY.UN_CODE=="156",]  ## Keep China only
  }
  
  names <- read_csv("CL_FI_SPECIES_GROUPS.csv")
  names <- names[,c("3A_Code","Scientific_Name")]
  colnames(production)[2] <- colnames(names)[1]
  datamerge <- inner_join(production,names,by='3A_Code')
  
  # check patterns
  evol <- datamerge[datamerge$MEASURE=="Q_tlw",]
  evol <- aggregate(evol$VALUE, by=list(evol$PERIOD,evol$PRODUCTION_SOURCE_DET.CODE),sum)
  
  p <- ggplot(evol, aes(x=Group.1,y=(x),color=factor(Group.2))) + geom_line(size=1.2) + theme(text = element_text(size = 26)) +
    xlab("Time") + ylab("Production (Tonnes - live weight)") + labs(color = "Source") + scale_color_manual(values=c("#213F99","#F7931D"))
  plot(p)
  ggsave(paste(id,"evol.pdf",sep=""), width = 12, height = 9)
  
  #### trophic levels (from Fishbase and sealifebase)
  ## load data
  ecbase <- ecology(server="sealifebase")
  ecbase <- ecbase[,c("SpecCode","FoodTroph","FoodSeTroph")]
  ecbase <- ecbase[complete.cases(ecbase),]
  tabase <- load_taxa(server="sealifebase")
  ecbase2 <- ecology(server="fishbase")
  ecbase2 <- ecbase2[,c("SpecCode","FoodTroph","FoodSeTroph")]
  ecbase2 <- ecbase2[complete.cases(ecbase2),]
  tabase2 <- load_taxa(server="fishbase")
  
  taxmerge1 <- inner_join(tabase,ecbase,by="SpecCode")
  taxmerge2 <- inner_join(tabase2,ecbase2,by="SpecCode")
  taxmerge1 <- taxmerge1[,c("FoodTroph","FoodSeTroph","Species")]
  taxmerge2 <- taxmerge2[,c("FoodTroph","FoodSeTroph","Species")]
  
  taxmerge <- rbind(taxmerge1,taxmerge2)
  taxmerge <- aggregate(taxmerge[,1:2],by=list(taxmerge$Species),FUN=mean)
  colnames(taxmerge)[1] <- "Species"
  
  ### Merge FAO and Fish/Sealife Base
  
  
  dataquantity <- datamerge[datamerge$MEASURE=="Q_tlw",]
  dataquantity <- aggregate(dataquantity$VALUE,
                            by=list(dataquantity$PERIOD,dataquantity$PRODUCTION_SOURCE_DET.CODE,dataquantity$Scientific_Name),sum)
  colnames(dataquantity) <- c("Year","Source","Species","Quantity")
  dataquantity <- inner_join(dataquantity,taxmerge,by="Species")
  
  dataquantity$rtroph <- "[2,3)"
  dataquantity$rtroph[dataquantity$FoodTroph>=3 & dataquantity$FoodTroph<4] <- "[3-4)"
  dataquantity$rtroph[dataquantity$FoodTroph>=4] <- "[4,5]"
  
  
  evolt <- aggregate(dataquantity$Quantity, by=list(dataquantity$Year,dataquantity$rtroph,dataquantity$Source),sum)
  colnames(evolt) <- c("Year","rtroph","Source","Quantity")
  p <- ggplot(evolt[evolt$Source!="Capture",], aes(x=Year,y=(Quantity),color=factor(rtroph))) + geom_line(size=1.2) + theme(text = element_text(size = 26)) +
    xlab("Time") + ylab("Production (Tonnes - live weight)") + labs(color = "Trophic level") + scale_color_manual(values=c("#213F99","#F7931D","#A5C22D"))
  plot(p)
  ggsave(paste(id,"Trophqt_Aqua.pdf",sep=""), width = 12, height = 9)
  
  p <- ggplot(evolt[evolt$Source=="Capture",], aes(x=Year,y=(Quantity),color=factor(rtroph))) + geom_line(size=1.2) + theme(text = element_text(size = 26)) +
    xlab("Time") + ylab("Production (Tonnes - live weight)") + labs(color = "Trophic level") + scale_color_manual(values=c("#213F99","#F7931D","#A5C22D"))
  plot(p)
  ggsave(paste(id,"Trophqt_Capt.pdf",sep=""), width = 12, height = 9)
  
}

