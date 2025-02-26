#Plant data
plant <- readRDS(file = here("data","intermediate","plant.RDS"))

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
summary(aov(CAT_activity~SoilWaterContent*Inoculum+Bloc, data=plant)) #Only block significant
summary(aov(SOD_activity~SoilWaterContent*Inoculum+Bloc, data=plant)) #Interaction term significant
TukeyHSD(aov(SOD_activity~SoilWaterContent*Inoculum+Bloc, data=plant)) #Nothing in P-adjust
summary(aov(RWC~SoilWaterContent*Inoculum+Bloc, data=plant)) #NS
summary(aov(Leaf_WC~SoilWaterContent*Inoculum+Bloc, data=plant)) #Soil water content significant
summary(aov(Fresh_biomass~SoilWaterContent*Inoculum+Bloc, data=plant)) #SWC and block significant
summary(aov(Dry_biomass~SoilWaterContent*Inoculum+Bloc, data=plant)) #Bloc significant
summary(aov(PWC~SoilWaterContent*Inoculum+Bloc, data=plant)) #SWC significant

#Reorder treatments
plant$Inoculum <- factor(plant$Inoculum, levels = c("ctrl", "irr", "amb"))

#Long format
plant.long <- gather(plant, variable, measurement, c(4:10))
plant.long$variable <- factor(plant.long$variable, levels = c("CAT_activity",
                                                              "SOD_activity", 
                                                              "Dry_biomass", 
                                                              "Fresh_biomass",
                                                              "Leaf_WC", "PWC",
                                                              "RWC"))

#boxplot
box.plant <- ggplot(plant.long[plant.long$variable != "RWC",], aes(x=SoilWaterContent, y=measurement, fill=Inoculum, colour = Inoculum, alpha = 0.75))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25))+
  facet_wrap(~variable, scales = "free_y", strip.position = "left", nrow = 3, 
             labeller = as_labeller(c("CAT_activity"="CAT activity (nmol/min/ml)", 
                                      "SOD_activity"="SOD activity (nmol/min/ml)",
                                      "Dry_biomass"="Dry biomass (g)", "Fresh_biomass"="Fresh biomass (g)",
                                      "Leaf_WC"="Leaf water content (%)", "PWC"= "Plant water content (%)"))) +
  scale_alpha_continuous(guide='none')+
  scale_colour_manual(guide='none', values=color6[c(1,4,5)])+
  scale_fill_manual(labels = c("Uninoculated", "Intermittent", "Continuous"), values=color6[c(1,4,5)])+
  theme_bw()+
  ylab(NULL)+
  xlab("Soil water content")+
  theme(strip.background = element_blank(), strip.placement = "outside")
box.plant
saveRDS(box.plant, file=here("data", "intermediate", "box.plant.RDS"))
