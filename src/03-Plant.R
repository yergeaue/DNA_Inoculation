#Plant data
setwd("C:/Users/etienne yergeau/OneDrive/Documents/IAF-INRS/Projects/20161130CharlotteGiardLaliberte/Results-DNARecruitment/plant/")
plant=read.table(file=paste(path,"/plant/plant.txt", sep=""), header=T, sep="\t",  row.names=1)
#check for normality
shapiro.test(plant$CAT_activity) #OK P=0.1306
shapiro.test(plant$SOD_activity) #OK P=0.1114
shapiro.test(plant$RWC) #OK P=0.5372
shapiro.test(plant$Leaf_WC) #OK P=0.5768
shapiro.test(plant$Fresh_biomass) #OK P=0.2275
shapiro.test(plant$Dry_biomass) #OK P=0.2072
shapiro.test(plant$PWC) #OK P=0.1477
#Check for homogenetiy of variance 
bartlett.test(plant$CAT_activity, plant$SoilWaterContent) #OK P=0.203
bartlett.test(plant$CAT_activity, plant$Inoculum) #OK P=0.4888
bartlett.test(plant$SOD_activity, plant$SoilWaterContent) #OK P=0.1788
bartlett.test(plant$SOD_activity, plant$Inoculum) #OK P=0.3866
bartlett.test(plant$RWC, plant$SoilWaterContent) #not OK, P=0.006191
bartlett.test(plant$RWC, plant$Inoculum) #OK, P=0.7172
bartlett.test(plant$Leaf_WC, plant$SoilWaterContent) #not OK P=0.005543
bartlett.test(plant$Leaf_WC, plant$Inoculum) #not OK P=0.5046
bartlett.test(plant$Fresh_biomass, plant$SoilWaterContent) #OK P=0.663
bartlett.test(plant$Fresh_biomass, plant$Inoculum) #OK P=0.9607
bartlett.test(plant$Dry_biomass, plant$SoilWaterContent) #OK P=0.3753
bartlett.test(plant$Dry_biomass, plant$Inoculum) #OK P=0.7077
bartlett.test(plant$PWC, plant$SoilWaterContent) #not OK P=0.01869
bartlett.test(plant$PWC, plant$Inoculum) #OK P=0.4854
#Anovas
summary(aov(CAT_activity~SoilWaterContent*Inoculum+Bloc, data=plant))
summary(aov(SOD_activity~SoilWaterContent*Inoculum+Bloc, data=plant))
TukeyHSD(aov(SOD_activity~SoilWaterContent*Inoculum+Bloc, data=plant))
summary(aov(RWC~SoilWaterContent*Inoculum+Bloc, data=plant))
summary(aov(Leaf_WC~SoilWaterContent*Inoculum+Bloc, data=plant))
summary(aov(Fresh_biomass~SoilWaterContent*Inoculum+Bloc, data=plant))
summary(aov(Dry_biomass~SoilWaterContent*Inoculum+Bloc, data=plant))
summary(aov(PWC~SoilWaterContent*Inoculum+Bloc, data=plant))

#ggplot for boxplot
library(ggplot2)
box.SOD=ggplot(plant, aes(x=SoilWaterContent, y=SOD_activity, fill=Inoculum))+
  geom_boxplot()+
  scale_fill_manual(values=c("black", "lightgrey", "white"))+
  theme_bw()+
  ylab("SOD activity (")

box.SOD
ggsave(file="SOD.eps", box.SOD, width = 7, units = "in")
ggsave(file="SOD.tiff", box.SOD, width = 7, units = "in", dpi=600)

box.PWC=ggplot(plant, aes(x=SoilWaterContent, y=PWC, fill=Inoculum))+
  geom_boxplot()+
  scale_fill_manual(values=c("black", "lightgrey", "white"))+
  theme_bw()+
  ylab("Plant water content (%)")

box.PWC
ggsave(file="PWC.eps", box.PWC, width = 7, units = "in")
ggsave(file="PWC.tiff", box.PWC, width = 7, units = "in", dpi=600)
