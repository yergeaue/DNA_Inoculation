##Make graphical representation and stat tests of taxonomy and functions of transferred genes vs. all assembly
genes.all <- readRDS(here("data","intermediate", "genes.all.RDS"))
annot.all <- readRDS(here("data", "intermediate", "annot.all.RDS"))
map <- readRDS(here("data", "intermediate", "map.RDS"))
map$Sample <- row.names(map)
genes.inoc.rel.absent <- readRDS(here("data", "intermediate", "gene.inoc.rel.absent.RDS"))
annot.inoc.def <- readRDS(here("data","intermediate","annot.inoc.def"))
dil <- 850/(50*100/0.5*1500/1000) #Dilution factor - added 850 ug in 15,000ug
margin <- 10

##Taxonomy
#Get number of genes for each phylum in the "all" assembly. To compare to transferred, include only one with KEGG entry available
genes.all.tax <- genes.all %>%
  mutate(gene_id = row.names(.)) %>%
  left_join(annot.all) %>%
  filter(kegg_entry != "NULL") %>%
  group_by(tax_phylum) %>%
  summarise(across(1:40, ~ sum(.x !=0))) %>%
  mutate(across(2:41, ~ ./sum(.)))

#Get the number of genes for each phylum among transferred genes - I
genes.inoc.i.tax <- genes.inoc.rel.absent %>%
  mutate(gene_id = row.names(.)) %>%
  left_join(annot.inoc.def) %>%
  group_by(tax_phylum) %>%
  summarise(across(1:38, ~ sum(.x >margin*dil*i.i))) %>%
  select(contains("tax") | contains("I1")) %>%
  mutate(across(2:13, ~ ./sum(.)))

#Get the number of genes for each phylum among transferred genes - NI
genes.inoc.ni.tax <- genes.inoc.rel.absent %>%
  mutate(gene_id = row.names(.)) %>%
  left_join(annot.inoc.def) %>%
  group_by(tax_phylum) %>%
  summarise(across(1:38, ~ sum(.x >margin*dil*i.ni))) %>%
  select(contains("tax") | contains("I2")) %>%
  mutate(across(2:13, ~ ./sum(.)))

#Create long file with all samples for plotting
genes.tax.long <- left_join(genes.inoc.i.tax,genes.inoc.ni.tax) %>%
  left_join(genes.all.tax, by = join_by(tax_phylum), suffix = c(".tr",".all")) %>%
  filter(T2.I2.3.all>0.005) %>%
  bind_rows(data.frame(tax_phylum = "Others", t(1-colSums(.[,2:65])))) %>%
  gather(sample,RelAbund,2:65) %>% #transform in long format for ggplot
  mutate(dataset=c(rep("transferred",264), rep("all",440))) %>%
  mutate(Sample=gsub("(.all)|(.tr)","", .$sample)) %>%
  left_join(map)

palette(c(brewer.pal(n=9, name="Set1"),"lightgrey", "black", "darkred", "darkblue", "darkgreen", "purple4"))
stack.trans.tax <- genes.tax.long %>% 
  filter(Inoculum != "ctrl" & SoilWaterContent %in% c("15%","50%", "Inoculum")) %>%
  ggplot(., aes(fill=tax_phylum, y=RelAbund, x=dataset)) + 
  geom_bar( stat="identity", position="fill")+
  ylab("Relative abundance")+  
  scale_fill_manual(values=palette(), guide=guide_legend(label.theme = element_text(face="italic", size=8)))+
  facet_grid(.~SoilWaterContent*Inoculum, scales="free_x", space="free_x")+
  theme_bw()+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

stack.trans.tax

#Test for significance
#Normality assumption
genes.tax.long %>% 
  filter(Inoculum != "ctrl" & SoilWaterContent %in% c("15%","50%")) %>%
  group_by(tax_phylum) %>%
  shapiro_test(RelAbund) 

#Equality of variances assumption
genes.tax.long %>% 
  filter(Inoculum != "ctrl" & SoilWaterContent %in% c("15%","50%")) %>%
  group_by(tax_phylum) %>%
  levene_test(RelAbund~Inoculum*SoilWaterContent) #ok

#Computing the statistical test - paired t-test on effect of inoculation separetely for 15% and 50% 
stat.test <- genes.tax.long %>% 
  filter(Inoculum != "ctrl" & SoilWaterContent %in% c("15%","50%")) %>%
  group_by(interaction(tax_phylum, SoilWaterContent,Inoculum)) %>%
  t_test(RelAbund ~ dataset, paired = TRUE)
stat.test #Many significant - see output

##Functions
#Get number of genes for each phylum in the "all" assembly. To compare to transferred, include only one with KEGG entry available
genes.all.COG <- genes.all %>%
  mutate(gene_id = row.names(.)) %>%
  left_join(annot.all) %>%
  filter(kegg_entry != "NULL") %>%
  group_by(cog_category) %>%
  summarise(across(1:40, ~ sum(.x !=0))) %>%
  mutate(across(2:41, ~ ./sum(.)))

#Get the number of genes for each phylum among transferred genes - I
genes.inoc.i.COG <- genes.inoc.rel.absent %>%
  mutate(gene_id = row.names(.)) %>%
  left_join(annot.inoc.def) %>%
  group_by(kegg_pathway_desc) %>%
  summarise(across(1:38, ~ sum(.x >margin*dil*i.i))) %>%
  select(contains("") | contains("I1")) %>%
  mutate(across(2:13, ~ ./sum(.)))

#Get the number of genes for each phylum among transferred genes - NI
genes.inoc.ni.COG <- genes.inoc.rel.absent %>%
  mutate(gene_id = row.names(.)) %>%
  left_join(annot.inoc.def) %>%
  group_by(cog_category) %>%
  summarise(across(1:38, ~ sum(.x >margin*dil*i.ni))) %>%
  select(contains("COG") | contains("I2")) %>%
  mutate(across(2:13, ~ ./sum(.)))

#Create long file with all samples for plotting
genes.COG.long <- left_join(genes.inoc.i.COG,genes.inoc.ni.COG) %>%
  left_join(genes.all.COG, by = join_by(cog_category), suffix = c(".tr",".all")) %>%
  filter(T2.I2.3.all>0.005) %>%
  bind_rows(data.frame(cog_category = "Others", t(1-colSums(.[,2:65])))) %>%
  gather(sample,RelAbund,2:65) %>% #transform in long format for ggplot
  mutate(dataset=c(rep("transferred",264), rep("all",440))) %>%
  mutate(Sample=gsub("(.all)|(.tr)","", .$sample)) %>%
  left_join(map)

palette(c(brewer.pal(n=9, name="Set1"),"lightgrey", "black", "darkred", "darkblue", "darkgreen", "purple4"))
stack.trans.COG <- genes.COG.long %>% 
  filter(Inoculum != "ctrl" & SoilWaterContent %in% c("15%","50%", "Inoculum")) %>%
  ggplot(., aes(fill=cog_category, y=RelAbund, x=dataset)) + 
  geom_bar( stat="identity", position="fill")+
  ylab("Relative abundance")+  
  scale_fill_manual(values=palette(), guide=guide_legend(label.theme = element_text(face="italic", size=8)))+
  facet_grid(.~SoilWaterContent*Inoculum, scales="free_x", space="free_x")+
  theme_bw()+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

stack.trans.COG

#Test for significance
#Normality assumption
genes.COG.long %>% 
  filter(Inoculum != "ctrl" & SoilWaterContent %in% c("15%","50%")) %>%
  group_by(cog_category) %>%
  shapiro_test(RelAbund) 

#Equality of variances assumption
genes.COG.long %>% 
  filter(Inoculum != "ctrl" & SoilWaterContent %in% c("15%","50%")) %>%
  group_by(cog_category) %>%
  levene_test(RelAbund~Inoculum*SoilWaterContent) #ok

#Computing the statistical test - paired t-test on effect of inoculation separetely for 15% and 50% 
stat.test <- genes.COG.long %>% 
  filter(Inoculum != "ctrl" & SoilWaterContent %in% c("15%","50%")) %>%
  group_by(interaction(cog_category, SoilWaterContent,Inoculum)) %>%
  t_test(RelAbund ~ dataset, paired = TRUE)
stat.test.COG #Many significant - see output

