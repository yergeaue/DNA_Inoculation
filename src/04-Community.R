##Taxonomy

head(phyla)
phyla.long=gather(phyla,Taxa,RelAbund,4:18) #transform in long format for ggplot
head(phyla.long)
phyla.long$SoilWaterContent=factor(phyla.long$SoilWaterContent, c("LW", "HW", "Soil", "Inoculum"))
palette(c(brewer.pal(n=9, name="Set1"),"lightgrey", "black", "darkred", "darkblue", "darkgreen", "purple4"))
stack=ggplot(phyla.long, aes(fill=Taxa, y=RelAbund, x=Inoculum)) + 
  geom_bar( stat="identity", position="fill")+
  ylab("Relative abundance")+  
  scale_fill_manual(values=palette(), guide=guide_legend(label.theme = element_text(face="italic", size=8)))+
  facet_grid(.~SoilWaterContent, scales="free_x", space="free_x")+
  theme_bw()+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

stack
ggsave(file="stackbarchart.eps", stack, width = 7, units = "in")
ggsave(file="stackbarchart.tiff", stack, width = 7, units = "in", dpi=600)

#ANOVAs
phyla.aov=phyla[5:36,]
#check for normality
shapiro.test(phyla.aov$Nitrospirae) #P=0.000111
shapiro.test(phyla.aov$Bacteroidetes) #P=0.00285
shapiro.test(phyla.aov$Firmicutes) #P=0.631
shapiro.test(phyla.aov$Planctomycetes) #P=6.6e-05
shapiro.test(phyla.aov$Chloroflexi) #P=0.000177
shapiro.test(phyla.aov$Acidobacteria) #P=1.51e-05
shapiro.test(phyla.aov$Gemmatimonadetes) #P=3.28e-05
shapiro.test(phyla.aov$Thaumarchaeota) #P=0.193
shapiro.test(phyla.aov$Alphaproteobacteria) #P=4.02e-05
shapiro.test(phyla.aov$Betaproteobacteria) #P=3.43e-05
shapiro.test(phyla.aov$Deltaproteobacteria) #P=3.8e-05
shapiro.test(phyla.aov$Gammaproteobacteria) #P=3.91e-05
shapiro.test(phyla.aov$Actinobacteria) #P=3.21e-05

#Log transform
shapiro.test(log(phyla.aov$Nitrospirae)) #P=0.000358
shapiro.test(log(phyla.aov$Bacteroidetes)) #P=0.00366
shapiro.test(log(phyla.aov$Planctomycetes)) #P=0.00014
shapiro.test(log(phyla.aov$Chloroflexi)) #P=0.000137
shapiro.test(log(phyla.aov$Acidobacteria)) #P=5.18e-05
shapiro.test(log(phyla.aov$Gemmatimonadetes)) #P=5.9e-05
shapiro.test(log(phyla.aov$Alphaproteobacteria)) #P=9.39e-05
shapiro.test(log(phyla.aov$Betaproteobacteria)) #P=5.45e-05
shapiro.test(log(phyla.aov$Deltaproteobacteria)) #P=5.39e-05
shapiro.test(log(phyla.aov$Gammaproteobacteria)) #P=6.02e-05
shapiro.test(log(phyla.aov$Actinobacteria)) #P=2.35e-05

#Sqrt transform
shapiro.test(sqrt(phyla.aov$Nitrospirae)) #P=0.000204
shapiro.test(sqrt(phyla.aov$Bacteroidetes)) #P=0.00328
shapiro.test(sqrt(phyla.aov$Planctomycetes)) #P=9.6e-05
shapiro.test(sqrt(phyla.aov$Chloroflexi)) #P=0.000155
shapiro.test(sqrt(phyla.aov$Acidobacteria)) #P=2.72e-05
shapiro.test(sqrt(phyla.aov$Gemmatimonadetes)) #P=4.36e-05
shapiro.test(sqrt(phyla.aov$Alphaproteobacteria)) #P=6.11e-05
shapiro.test(sqrt(phyla.aov$Betaproteobacteria)) #P=4.3e-05
shapiro.test(sqrt(phyla.aov$Deltaproteobacteria)) #P=4.51e-05
shapiro.test(sqrt(phyla.aov$Gammaproteobacteria)) #P=4.84e-05
shapiro.test(sqrt(phyla.aov$Actinobacteria)) #P=2.74e-05

#Anovas
summary(aov(Firmicutes~SoilWaterContent*Inoculum+Bloc, data=phyla.aov)) 
summary(aov(Thaumarchaeota~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
summary(aov(Nitrospirae~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
summary(aov(Bacteroidetes~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
summary(aov(Planctomycetes~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
summary(aov(Chloroflexi~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
summary(aov(Acidobacteria~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
summary(aov(Gemmatimonadetes~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
summary(aov(Alphaproteobacteria~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
summary(aov(Betaproteobacteria~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
summary(aov(Deltaproteobacteria~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
summary(aov(Gammaproteobacteria~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
summary(aov(Actinobacteria~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))

#TukeyHSD
TukeyHSD(aov(Firmicutes~SoilWaterContent*Inoculum+Bloc, data=phyla.aov)) 
TukeyHSD(aov(Thaumarchaeota~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
TukeyHSD(aov(Nitrospirae~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
TukeyHSD(aov(Bacteroidetes~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
TukeyHSD(aov(Planctomycetes~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
TukeyHSD(aov(Chloroflexi~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
TukeyHSD(aov(Acidobacteria~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
TukeyHSD(aov(Gemmatimonadetes~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
TukeyHSD(aov(Alphaproteobacteria~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
TukeyHSD(aov(Betaproteobacteria~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
TukeyHSD(aov(Deltaproteobacteria~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
TukeyHSD(aov(Gammaproteobacteria~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))
TukeyHSD(aov(Actinobacteria~SoilWaterContent*Inoculum+Bloc, data=phyla.aov))


#Kruskal tests
#Soil water content
kruskal.test(Nitrospirae~SoilWaterContent, data=phyla.aov)
kruskal.test(Bacteroidetes~SoilWaterContent, data=phyla.aov)
kruskal.test(Planctomycetes~SoilWaterContent, data=phyla.aov)
kruskal.test(Chloroflexi~SoilWaterContent, data=phyla.aov)
kruskal.test(Acidobacteria~SoilWaterContent, data=phyla.aov)
kruskal.test(Gemmatimonadetes~SoilWaterContent, data=phyla.aov)
kruskal.test(Alphaproteobacteria~SoilWaterContent, data=phyla.aov)
kruskal.test(Betaproteobacteria~SoilWaterContent, data=phyla.aov)
kruskal.test(Deltaproteobacteria~SoilWaterContent, data=phyla.aov)
kruskal.test(Gammaproteobacteria~SoilWaterContent, data=phyla.aov)
kruskal.test(Actinobacteria~SoilWaterContent, data=phyla.aov)

#Inoculum
kruskal.test(Nitrospirae~Inoculum, data=phyla.aov)
kruskal.test(Bacteroidetes~Inoculum, data=phyla.aov)
kruskal.test(Planctomycetes~Inoculum, data=phyla.aov)
kruskal.test(Chloroflexi~Inoculum, data=phyla.aov)
kruskal.test(Acidobacteria~Inoculum, data=phyla.aov)
kruskal.test(Gemmatimonadetes~Inoculum, data=phyla.aov)
kruskal.test(Alphaproteobacteria~Inoculum, data=phyla.aov)
kruskal.test(Betaproteobacteria~Inoculum, data=phyla.aov)
kruskal.test(Deltaproteobacteria~Inoculum, data=phyla.aov)
kruskal.test(Gammaproteobacteria~Inoculum, data=phyla.aov)
kruskal.test(Actinobacteria~Inoculum, data=phyla.aov)

#Soil water content*Inoculum
kruskal.test(Nitrospirae~interaction(phyla.aov$SoilWaterContent,phyla.aov$Inoculum), data=phyla.aov)
kruskal.test(Bacteroidetes~interaction(phyla.aov$SoilWaterContent,phyla.aov$Inoculum), data=phyla.aov)
kruskal.test(Planctomycetes~interaction(phyla.aov$SoilWaterContent,phyla.aov$Inoculum), data=phyla.aov)
kruskal.test(Chloroflexi~interaction(phyla.aov$SoilWaterContent,phyla.aov$Inoculum), data=phyla.aov)
kruskal.test(Acidobacteria~interaction(phyla.aov$SoilWaterContent,phyla.aov$Inoculum), data=phyla.aov)
kruskal.test(Gemmatimonadetes~interaction(phyla.aov$SoilWaterContent,phyla.aov$Inoculum), data=phyla.aov)
kruskal.test(Alphaproteobacteria~interaction(phyla.aov$SoilWaterContent,phyla.aov$Inoculum), data=phyla.aov)
kruskal.test(Betaproteobacteria~interaction(phyla.aov$SoilWaterContent,phyla.aov$Inoculum), data=phyla.aov)
kruskal.test(Deltaproteobacteria~interaction(phyla.aov$SoilWaterContent,phyla.aov$Inoculum), data=phyla.aov)
kruskal.test(Gammaproteobacteria~interaction(phyla.aov$SoilWaterContent,phyla.aov$Inoculum), data=phyla.aov)
kruskal.test(Actinobacteria~interaction(phyla.aov$SoilWaterContent,phyla.aov$Inoculum), data=phyla.aov)


##COG
#Load files
COG.all <- readRDS(here("data", "intermediate", "COG.all.RDS"))
genes.all.rel <- readRDS(here("data", "intermediate", "genes.all.rel.RDS"))
#Merge annotations and relabund
sum(rownames(COG.all) == colnames(genes.all.rel), na.rm=TRUE) #
COG.relabund.all <- cbind(genes.all.rel, COG.all)

#Create summary data frame
COG.relabund.all.sum <- COG.relabund.all %>%
  group_by(cog_category) %>%
  summarise(
    n = n()#,
    #over = sum(log2 > 0, na.rm = TRUE),
    #under = sum(log2 < 0, na.rm = TRUE),
    
  ) %>%
  filter(
    n > 100
  )
COG.rhizo <- COG.rhizo[COG.rhizo$cog_category %in% COG.root$COG_category,] #to have matching legends with ROOT; NULL: 13505 6061 7444
row.names(COG.rhizo) <- COG.rhizo$cog_category
others.COG.rhizo <- data.frame("Others", (21765 - 13505 - sum(COG.rhizo$n)), (10565 - 6061 - sum(COG.rhizo$over)), (11200 - 7444 - sum(COG.rhizo$under))) #Omit NULL from Others count
colnames(others.COG.rhizo) <- colnames(COG.rhizo)
row.names(others.COG.rhizo) <- "Others"
COG.rhizo <- rbind(COG.rhizo, others.COG.rhizo) 
colSums(COG.rhizo[,2:4]) #8260 4504 3756

#Get a column for all transcripts
annot.rhizo <- annot[annot$gene_id %in% row.names(MT.rhizo),]

COG.rhizo.all = annot.rhizo %>%
  group_by(cog_category) %>%
  summarise(
    count = n(),
    
    
  )
COG.rhizo.all <- COG.rhizo.all[COG.rhizo.all$cog_category %in% COG.rhizo$cog_category,]
row.names(COG.rhizo.all) <- COG.rhizo.all$cog_category
others.COG.rhizo.all <- data.frame("Others", (1193501 - 799568 - sum(COG.rhizo.all$count))) #Omit NULL
colnames(others.COG.rhizo.all) <- colnames(COG.rhizo.all)
row.names(others.COG.rhizo.all) <- "Others"
COG.rhizo.all <- rbind(COG.rhizo.all, others.COG.rhizo.all) 
sum(COG.rhizo.all$count)

COG.rhizo.all$cog_category == COG.rhizo$cog_category
COG.rhizo <- cbind(COG.rhizo,COG.rhizo.all$count)
colnames(COG.rhizo) <- c("COG_category", "n", "DA+", "DA-", "All")

#Prepare for ggplot
COG.rhizo.long <- gather(COG.rhizo,subset,counts,3:5) #transform in long format for ggplot
comp.rhizo.cog <- data.frame(Compartment = rep("Rhizosphere",45))#Add a column for compartment
COG.rhizo.long <- cbind(COG.rhizo.long,comp.rhizo.cog)
#Put everything together
COG.long.both <- rbind(COG.rhizo.long, COG.root.long)


palette(c(brewer.pal(n = 9, name = "Set1"),"lightgrey", "black", "darkred", "darkblue", "darkgreen", "purple4", "brown3", "cyan"))
stack.COG <- ggplot(COG.long.both, aes(fill = COG_category, y = counts, x = subset)) + 
  geom_bar( stat = "identity", position = "fill") +
  ylab("Fraction of transcripts") + 
  scale_fill_manual(values = palette(), guide = guide_legend(title = "COG category")) +
  theme_bw() +
  scale_y_continuous( expand = c(0,0)) +
  scale_x_discrete(name = "Subset") +
  facet_grid(. ~ Compartment)
stack.COG 
