###Look for effects of the inoculum on the general community structure and composition
##Load files
map <- readRDS(here("data","intermediate","map.RDS"))
map.nosoilinoc <- map |> filter(SoilWaterContent != "Soil" & SoilWaterContent != "Inoculum")
genes.noninoc.rel <- readRDS(here("data","intermediate", "genes.noninoc.rel.RDS"))
annot.noninoc <- readRDS(file = here("data","intermediate", "annot.noninoc.RDS"))
genes.inoc.rel <- readRDS(here("data","intermediate", "genes.inoc.rel.RDS"))
annot.inoc <- readRDS(file = here("data","intermediate", "annot.inoc.RDS"))

##Permanova
colnames(genes.noninoc.rel)==row.names(map.nosoilinoc) #Check: all true
bray.noninoc <- vegdist(t(genes.noninoc.rel), method = "bray")
adonis2(as.matrix(bray.noninoc)~map.nosoilinoc$SoilWaterContent*map.nosoilinoc$Inoculum+map.nosoilinoc$Bloc,
        strata=map.nosoilinoc$Bloc, sqrt.dist = TRUE, by="terms")
saveRDS(bray.noninoc, here("data","intermediate", "bray.noninoc.RDS"))

##PCoA
pcoa.noninoc <- cmdscale(sqrt(as.dist(as.matrix(bray.noninoc))), k=2,eig=TRUE)
pcoa.noninoc$eig[1]/sum(pcoa.noninoc$eig)*100#Axis1 = 15.52845
pcoa.noninoc$eig[2]/sum(pcoa.noninoc$eig)*100#Axis2 = 2.825777

#Merge for ggplot 
row.names(map.nosoilinoc)==row.names(pcoa.noninoc$points) #All true
p<-cbind(map.nosoilinoc, pcoa.noninoc$points)
colnames(p) <- c("SoilWaterContent","Inoculum","Bloc","Axis1","Axis2")

#Plot using ggplot
ordi.noninoc <- ggplot(p, aes(x=Axis1, y=Axis2, colour=SoilWaterContent, shape = Inoculum)) + 
  geom_point()  + 
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Axis 2")+
  xlab("Axis 1")+
  scale_color_manual(values = palette()[c(3,1)], guide = guide_legend(title = "Soil water content"),
                     labels = c("LW", "HW"))+
  scale_shape_discrete(labels = c("Uninoculated", "Intermittent", "Continuous"))
ordi.noninoc #Major effect of inoculum - 4 outliers? 

##Taxonomy
#Make summary at the family level (top 20)
#Inoc
family.inoc <- genes.inoc.rel|>
  rownames_to_column(var = "gene_id") |>
  left_join(annot.inoc) |>
  filter(tax_family != "NULL") |>
  group_by(tax_family)|>
  summarise(across(2:37, ~ sum(.)))|>
  mutate(across(2:37, ~ ./sum(.)))
saveRDS(family.inoc, here("data","intermediate", "family.inoc.RDS"))

#Noninoc
family.noninoc <- genes.noninoc.rel|>
  rownames_to_column(var = "gene_id") |>
  left_join(annot.noninoc) |>
  filter(tax_family != "NULL") |>
  group_by(tax_family)|>
  summarise(across(2:37, ~ sum(.)))|>
  mutate(across(2:37, ~ ./sum(.)))|>
  left_join(select(family.inoc, tax_family, i.i, i.ni))|>
  rowwise() |>
  mutate(RowSums = sum(c_across(c(2,20,38,39))))|> #To have a better balance between inoculum and rhizosphere
  ungroup() |>
  slice_max(RowSums,n=19) |>
  pivot_longer(cols = -tax_family) |>
  pivot_wider(names_from = tax_family)|>
  rename("Sample" = 1) |>
  rowwise() |>
  mutate(Others = 1-sum(c_across(2:20)))|>
  ungroup() |>
  slice(-39) |>
  left_join(rownames_to_column(map, var = "Sample")) |>
  pivot_longer(cols = c(2:21), names_to = "Family", values_to = "RelAbund")
saveRDS(family.noninoc, here("data","intermediate", "family.noninoc.RDS"))

#Plot
family.noninoc$Inoculum=factor(family.noninoc$Inoculum, c("i", "ni", "ctrl"))
stack.all <- 
  ggplot(family.noninoc, aes(fill=Family, y=RelAbund, x=Inoculum)) + 
  geom_bar( stat="identity", position="fill")+
  ylab("Relative abundance")+  
  scale_fill_manual(values=palette(), guide=guide_legend(label.theme = element_text(face="italic", size=8)))+
  scale_x_discrete(labels=c("Intermittent", "Continuous", "Uninoculated"),  )+
  facet_grid(.~SoilWaterContent, scales="free_x", space="free_x", labeller = as_labeller(c("15%"="LW", "50%"="HW", "Inoculum"="Inoculum")))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

stack.all

#Test for significance of inoculation
#Normality assumption
family.noninoc %>% 
  filter(SoilWaterContent %in% c("15%","50%")) %>%
  group_by(Family) %>%
  shapiro_test(RelAbund) #Most not ok.

#Equality of variances assumption
trans.phylum.all %>% 
  filter(SoilWaterContent %in% c("15%","50%")) %>%
  group_by(Phylum) %>%
  levene_test(RelAbund~Inoculum*SoilWaterContent) #A few not OK.

#Computing the statistical test - paired t-test on effect of inoculation  
#15%-continuous
family.noninoc |> 
  filter(SoilWaterContent %in% c("15%")) |>
  filter(Inoculum %in% c("ni","ctrl")) |>
  group_by(Family) |>
  arrange(Inoculum,Sample) |>
  wilcox_test(RelAbund ~ Inoculum, paired = TRUE)
 #Nitrososphaeraceae 0.0625

#50%-continuous
family.noninoc |> 
  filter(SoilWaterContent %in% c("50%")) |>
  filter(Inoculum %in% c("ni","ctrl")) |>
  group_by(Family) |>
  arrange(Inoculum,Sample) |>
  wilcox_test(RelAbund ~ Inoculum, paired = TRUE)
#Hyphomicrobiaceae  0.0313; Mycobacteriaceae 0.0938; Nocardioidaceae 0.0625

#15%-intermittent
family.noninoc |> 
  filter(SoilWaterContent %in% c("15%")) |>
  filter(Inoculum %in% c("i","ctrl")) |>
  group_by(Family) |>
  arrange(Inoculum,Sample) |>
  wilcox_test(RelAbund ~ Inoculum, paired = TRUE)
#Nitrososphaeraceae 0.0938

#50%-continuous
family.noninoc |> 
  filter(SoilWaterContent %in% c("50%")) |>
  filter(Inoculum %in% c("i","ctrl")) |>
  group_by(Family) |>
  arrange(Inoculum,Sample) |>
  wilcox_test(RelAbund ~ Inoculum, paired = TRUE)
#None

#Means for text
print(family.noninoc |> 
  filter(SoilWaterContent %in% c("15%","50%")) |>
  filter(Family %in% c("Hyphomicrobiaceae", "Mycobacteriaceae", "Nocardioidaceae", "Nitrososphaeraceae")) |>
  group_by(Family, Inoculum, SoilWaterContent) |>
  summarise(100*mean(RelAbund)), n=24 )
