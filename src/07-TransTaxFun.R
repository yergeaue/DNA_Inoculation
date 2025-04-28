##Make graphical representation and stat tests of taxonomy and functions of transferred genes vs. all genes (in inoculum assembly)
genes.inoc.rel <- readRDS(here("data","intermediate", "genes.inoc.rel.RDS"))
annot.inoc <- readRDS(here("data","intermediate", "annot.inoc.RDS"))
genes.inoc.rel.absent <- readRDS(here("data", "intermediate", "gene.inoc.rel.absent.RDS"))
annot.inoc.def <- readRDS(here("data","intermediate","annot.inoc.def"))
map <- readRDS(here("data", "intermediate", "map.RDS"))
map$Sample <- row.names(map)
dil <- 850/(50*100/0.5*1500/1000) #Dilution factor - added 850 ug in 15,000ug
margin <- 10


##Taxonomy
#Get number of genes for each phylum in the "all" assembly. To compare to transferred, include only one with KEGG entry available
genes.inoc.tax <- genes.inoc.rel %>%
  mutate(gene_id = row.names(.)) %>%
  left_join(annot.inoc) %>%
  filter(kegg_entry != "NULL") %>%
  group_by(tax_phylum) %>%
  summarise(across(1:38, ~ sum(.x !=0))) %>%
  mutate(across(2:39, ~ ./sum(.)))

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
  left_join(genes.inoc.tax, by = join_by(tax_phylum), suffix = c(".tr",".all")) %>%
  filter(T2.I2.3.all>0.005) %>%
  bind_rows(data.frame(tax_phylum = "Others", t(1-colSums(.[,2:63])))) %>%
  gather(sample,RelAbund,2:63) %>% #transform in long format for ggplot
  mutate(Dataset=c(rep("transferred",240), rep("all",380))) %>%
  mutate(Sample=gsub("(.all)|(.tr)","", .$sample)) %>%
  left_join(map)

palette(c(brewer.pal(n=9, name="Set1"),"lightgrey", "black", "darkred", "darkblue", "darkgreen", "purple4"))
stack.trans.tax <- genes.tax.long %>% 
  filter(Inoculum != "ctrl" & SoilWaterContent %in% c("15%","50%", "Inoculum")) %>%
  ggplot(., aes(fill=tax_phylum, y=RelAbund, x=Dataset)) + 
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

#Computing the statistical test - paired t-test on effect of inoculation  
stat.test.tax <- genes.tax.long %>% 
  filter(Inoculum != "ctrl" & SoilWaterContent %in% c("15%","50%")) %>%
  group_by(tax_phylum) %>%
  wilcox_test(RelAbund ~ Dataset, paired = TRUE)
stat.test.tax #Many significant - see output

##Functions
#Get number of genes for each functional category in the "inoc" assembly. To compare to transferred, include only one with KEGG entry available
genes.inoc.COG <- genes.inoc.rel %>%
  mutate(gene_id = row.names(.)) %>%
  left_join(annot.inoc) %>%
  filter(kegg_entry != "NULL") %>%
  group_by(kegg_definition) %>%
  summarise(across(1:38, ~ sum(.x !=0))) %>%
  mutate(across(2:39, ~ ./sum(.)))

#Get the number of genes for each phylum among transferred genes - I
genes.inoc.i.COG <- genes.inoc.rel.absent %>%
  mutate(gene_id = row.names(.)) %>%
  left_join(annot.inoc.def) %>%
  group_by(kegg_definition) %>%
  summarise(across(1:38, ~ sum(.x >margin*dil*i.i))) %>%
  select(contains("kegg") | contains("I1")) %>%
  mutate(across(2:13, ~ ./sum(.)))

#Get the number of genes for each phylum among transferred genes - NI
genes.inoc.ni.COG <- genes.inoc.rel.absent %>%
  mutate(gene_id = row.names(.)) %>%
  left_join(annot.inoc.def) %>%
  group_by(kegg_definition) %>%
  summarise(across(1:38, ~ sum(.x >margin*dil*i.ni))) %>%
  select(contains("kegg") | contains("I2")) %>%
  mutate(across(2:13, ~ ./sum(.)))

#Create long file with all samples for plotting
genes.COG.long <- left_join(genes.inoc.i.COG,genes.inoc.ni.COG) %>%
  left_join(genes.inoc.COG, by = join_by(kegg_definition), suffix = c(".tr",".all")) %>%
  filter(rowMeans(.[,2:63])>0.0025) %>%
  #bind_rows(data.frame(kegg_definition = "Others", t(1-colSums(.[,2:63])))) %>%
  gather(sample,RelAbund,2:63) %>% #transform in long format for ggplot
  mutate(Dataset=c(rep("transferred",336), rep("all",532))) %>%
  mutate(Sample=gsub("(.all)|(.tr)","", .$sample)) %>%
  left_join(map)

palette(c(brewer.pal(n=9, name="Set1"),"lightgrey", "black", "darkred", "darkblue", "darkgreen", "purple4", "yellow", "pink", "orange", "blue", "green", "red", "black", "grey"))
stack.trans.COG <- genes.COG.long %>% 
  filter(Inoculum != "ctrl" & SoilWaterContent %in% c("15%","50%", "Inoculum")) %>%
  ggplot(., aes(fill=kegg_definition, y=RelAbund, x=Dataset)) + 
  geom_bar( stat="identity", position="stack")+
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
  group_by(kegg_definition) %>%
  shapiro_test(RelAbund) 

#Equality of variances assumption
genes.COG.long %>% 
  filter(Inoculum != "ctrl" & SoilWaterContent %in% c("15%","50%")) %>%
  group_by(kegg_definition) %>%
  levene_test(RelAbund~Inoculum*SoilWaterContent) #ok

#Computing the statistical test - paired t-test on effect of inoculation separetely for 15% and 50% 
stat.test.COG <- genes.COG.long %>% 
  filter(Inoculum != "ctrl" & SoilWaterContent %in% c("15%","50%")) %>%
  group_by(kegg_definition) %>%
  wilcox_test(RelAbund ~ Dataset, paired = TRUE)
stat.test.COG #Many significant - see output

