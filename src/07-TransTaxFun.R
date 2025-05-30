##Make graphical representation and stat tests of taxonomy and functions of transferred genes vs. all genes
genes.inoc.rel <- readRDS(here("data","intermediate", "genes.inoc.rel.RDS"))
annot.inoc <- readRDS(here("data","intermediate", "annot.inoc.RDS"))
genes.noninoc.rel <- readRDS(here("data","intermediate", "genes.noninoc.rel.RDS"))
annot.noninoc <- readRDS(file = here("data","intermediate", "annot.noninoc.RDS"))
genes.inoc.rel.absent <- readRDS(here("data", "intermediate", "genes.inoc.rel.absent.RDS"))
#annot.inoc.def <- readRDS(here("data","intermediate","annot.inoc.def.RDS"))
tax.trans.0 <- readRDS(here("data", "intermediate", "tax.trans.0.RDS"))
map <- readRDS(here("data", "intermediate", "map.RDS"))
dil <- 850/(50*100/0.5*1500/1000) #Dilution factor - added 850 ug in 15,000 ug
margin <- 128

###Relative abundance in inoculum and rhizosphere vs. frequency in transferred genes (tax.trans.0.sum)
##Taxonomy
phylum.inoc <- genes.inoc.rel|>
  rownames_to_column(var = "gene_id") |>
  left_join(annot.inoc) |>
  #filter(tax_phylum != "NULL") |>
  group_by(tax_phylum)|>
  summarise(across(2:39, ~ sum(.)))|>
  mutate(across(2:39, ~ ./sum(.)))

phylum.noninoc <- genes.noninoc.rel|>
  rownames_to_column(var = "gene_id") |>
  left_join(annot.noninoc) |>
  #filter(tax_phylum != "NULL") |>
  group_by(tax_phylum)|>
  summarise(across(2:37, ~ sum(.)))|>
  mutate(across(2:37, ~ ./sum(.))) |>
  rowwise() |>
  mutate(HW = mean(c_across(2:19)), LW = mean(c_across(20:37)), .before = "T1.I1.1") |>
  ungroup()

tax.trans.0.sum #Compare to the two other data tables...values in the text

phylum.trans.all <- tax.trans.0.sum |>
  left_join(select(phylum.inoc, tax_phylum, i.i, i.ni), join_by(Phylum==tax_phylum)) |>
  left_join(select(phylum.noninoc, tax_phylum, LW, HW), join_by(Phylum==tax_phylum)) |>
  select(-sum) |>
  rename(Transferred=rel, Intermittent=i.i, Continuous=i.ni)|>
  filter(Transferred>0.01) |>
  pivot_longer(cols = -Phylum) |>
  pivot_wider(names_from = Phylum)|>
  rename("Dataset" = 1) |>
  rowwise() |>
  mutate(Others = 1-sum(c_across(2:11)))|>
  ungroup() |>
  pivot_longer(cols = c(2:12), names_to = "Phylum", values_to = "RelAbund")

phylum.trans.all$Dataset=factor(phylum.trans.all$Dataset, c("Transferred", "Intermittent", "Continuous", "LW", "HW"))
stack.trans.tax  <- 
  ggplot(phylum.trans.all, aes(fill=Phylum, y=RelAbund, x=Dataset)) + 
  geom_bar( stat="identity", position="stack")+
  ylab("Relative abundance")+  
  xlab("")+
  scale_fill_manual(values=palette(), guide=guide_legend(label.theme = element_text(face="italic", size=8)))+
  #scale_x_discrete(labels=c("Intermittent", "Continuous", "Uninoculated"),  )+
  #facet_grid(.~SoilWaterContent, scales="free_x", space="free_x", labeller = as_labeller(c("15%"="LW", "50%"="HW", "Inoculum"="Inoculum")))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

stack.trans.tax

##Functions
fun.inoc <- genes.inoc.rel|>
  rownames_to_column(var = "gene_id") |>
  left_join(annot.inoc) |>
  filter(kegg_pathway_desc !="NULL")|>
  filter(kegg_pathway_desc !="") |> 
  mutate(kegg_pathway_desc = gsub("==.*$","", kegg_pathway_desc)) |> #Get rid of multiple pathways - keep first (imperfect)
  group_by(kegg_pathway_desc) |>
  summarise(across(2:39, ~ sum(.)))|>
  mutate(across(2:39, ~ ./sum(.)))

fun.noninoc <- genes.noninoc.rel|>
  rownames_to_column(var = "gene_id") |>
  left_join(annot.noninoc) |>
  filter(kegg_pathway_desc !="NULL")|> 
  filter(kegg_pathway_desc !="") |> 
  mutate(kegg_pathway_desc = gsub("==.*$","", kegg_pathway_desc)) |> #Get rid of multiple pathways - keep first (imperfect)
  group_by(kegg_pathway_desc) |>
  summarise(across(2:37, ~ sum(.)))|>
  mutate(across(2:37, ~ ./sum(.))) |>
  rowwise() |>
  mutate(HW = mean(c_across(2:19)), LW = mean(c_across(20:37)), .before = "T1.I1.1") |>
  ungroup()

fun.trans.all <- fun.trans.0.sum |>
  left_join(select(fun.inoc, kegg_pathway_desc, i.i, i.ni), join_by(`KEGG pathway`==kegg_pathway_desc)) |>
  left_join(select(fun.noninoc, kegg_pathway_desc, LW, HW), join_by(`KEGG pathway`==kegg_pathway_desc)) |>
  select(-sum) |>
  rename(Transferred=rel, Intermittent=i.i, Continuous=i.ni)|>
  filter(Transferred>0.01) |>
  pivot_longer(cols = -`KEGG pathway`) |>
  pivot_wider(names_from = `KEGG pathway`)|>
  rename("Dataset" = 1) |>
  rowwise() |>
  mutate(Others = 1-sum(c_across(2:22)))|>
  ungroup() |>
  pivot_longer(cols = c(2:23), names_to = "KEGG Pathway", values_to = "RelAbund")

fun.trans.all$Dataset=factor(fun.trans.all$Dataset, c("Transferred", "Intermittent", "Continuous", "LW", "HW"))
stack.trans.fun  <- 
  ggplot(fun.trans.all, aes(fill=`KEGG Pathway`, y=RelAbund, x=Dataset)) + 
  geom_bar( stat="identity", position="stack")+
  ylab("Relative abundance")+  
  xlab("")+
  scale_fill_manual(values=c(palette(),"#9E3D22", "#2B5C8A"), guide=guide_legend(label.theme = element_text(size=8), ncol = 1))+
  #scale_x_discrete(labels=c("Intermittent", "Continuous", "Uninoculated"),  )+
  #facet_grid(.~SoilWaterContent, scales="free_x", space="free_x", labeller = as_labeller(c("15%"="LW", "50%"="HW", "Inoculum"="Inoculum")))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

stack.trans.fun




#Get relative abundance of phylum in the "inoculum" assembly.
genes.inoc.tax <- genes.inoc.rel %>%
  mutate(gene_id = row.names(.)) %>%
  left_join(annot.inoc) %>%
  #filter(kegg_entry != "NULL") %>%
  group_by(tax_phylum) %>%
  summarise(across(1:38, ~ sum(.x !=0))) %>%
  mutate(across(2:39, ~ ./sum(.)))

#Get relative abundance of phylum in the "non-inoculum" assembly.
genes.noninoc.tax <- genes.noninoc.rel %>%
  rownames_to_column(var="gene_id") %>%
  left_join(annot.noninoc) %>%
  #filter(kegg_entry != "NULL") %>%
  group_by(tax_phylum) %>%
  summarise(across(2:37, ~ sum(.x !=0))) %>%
  mutate(across(2:37, ~ ./sum(.)))

#Get fraction of genes assigned to each phylum among transferred genes - I
genes.inoc.i.tax <- genes.inoc.rel.absent %>%
  mutate(gene_id = row.names(.)) %>%
  left_join(annot.inoc) %>%
  group_by(tax_phylum) %>%
  summarise(across(3:38, ~ sum(.x >margin*dil*i.i))) %>%
  select(contains("tax") | contains("I1")) %>%
  mutate(across(2:13, ~ ./sum(.)))

#Get fraction of genes assigned to each phylum among transferred genes - NI
genes.inoc.ni.tax <- genes.inoc.rel.absent %>%
  mutate(gene_id = row.names(.)) %>%
  left_join(annot.inoc) %>%
  group_by(tax_phylum) %>%
  summarise(across(2:39, ~ sum(.x >margin*dil*i.ni))) %>%
  select(contains("tax") | contains("I2")) %>%
  mutate(across(2:13, ~ ./sum(.)))

#Create long file with all samples for plotting
genes.tax.long <- left_join(genes.inoc.i.tax,genes.inoc.ni.tax) %>%
  left_join(select(genes.inoc.tax, i.i, i.ni, tax_phylum), by = join_by(tax_phylum), suffix = c(".tr",".inoc")) %>%
  left_join(select(genes.noninoc.tax, !contains("I3")), by = join_by(tax_phylum), suffix = c("",".noninoc")) %>%
  rowwise() |>
  mutate(TransMean = mean(c_across(2:24)), .before=T1.I1.1)|>
  ungroup() |>
  filter(TransMean>0.01) %>%
  bind_rows(data.frame(tax_phylum = "Others", t(1-colSums(.[,3:52])))) %>%
  gather(sample,RelAbund,3:52) %>% #transform in long format for ggplot
  select(-TransMean) %>%
  mutate(Dataset=c(rep("transferred",288), rep ("inoculum",24), rep("rhizo",288))) %>%
  mutate(Sample=gsub("(.noninoc)|(.tr)","", .$sample)) %>%
  select(-sample) %>%
  left_join(rownames_to_column(map,var="Sample"))

genes.tax.long$Dataset=factor(genes.tax.long$Dataset, c("rhizo", "inoculum", "transferred"))
stack.trans.tax <-  ggplot(genes.tax.long, aes(fill=tax_phylum, y=RelAbund, x=Dataset)) + 
  geom_bar( stat="identity", position="fill")+
  ylab("Relative abundance")+ 
  xlab("")+
  scale_fill_manual(values=palette(), guide=guide_legend(label.theme = element_text(face="italic", size=8)))+
  facet_grid(.~Inoculum, scales="free_x", space="free_x", labeller = as_labeller(c("15%"="LW", "50%"="HW", "Inoculum"="Inoculum", "i"="Intermittent", "ni"="Continuous", "ctrl"="Uninoculated")))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank())+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

stack.trans.tax

#Test for significance vs. rhizosphere
#Normality assumption
genes.tax.long %>% 
  filter(SoilWaterContent %in% c("15%","50%")) %>%
  group_by(tax_phylum) %>%
  shapiro_test(RelAbund) #All not ok.

#Equality of variances assumption
genes.tax.long %>% 
  filter(SoilWaterContent %in% c("15%","50%")) %>%
  group_by(tax_phylum) %>%
  levene_test(RelAbund~Inoculum*SoilWaterContent) #A few not OK.

#Computing the statistical test - paired t-test on effect of inoculation  
stat.test.tax <- genes.tax.long %>% 
  filter(SoilWaterContent %in% c("15%","50%")) %>%
  group_by(tax_phylum) %>%
  wilcox_test(RelAbund ~ Dataset, paired = TRUE)
stat.test.tax #Many significant - see output

##Functions
#Get relative abundance of phylum in the "inoculum" assembly.
genes.inoc.tax <- genes.inoc.rel %>%
  mutate(gene_id = row.names(.)) %>%
  left_join(annot.inoc) %>%
  filter(kegg_pathway_desc !="NULL")|> 
  filter(kegg_pathway_desc !="") |> 
  mutate(kegg_pathway_desc = gsub("==.*$","", kegg_pathway_desc)) |> #Get rid of multiple pathways - keep first (imperfect)
  group_by(kegg_pathway_desc) |>
  summarise(across(2:39, ~ sum(.x !=0))) %>%
  mutate(across(2:39, ~ ./sum(.)))

#Get relative abundance of phylum in the "non-inoculum" assembly.
genes.noninoc.fun <- genes.noninoc.rel %>%
  rownames_to_column(var="gene_id") %>%
  left_join(annot.noninoc) %>%
  filter(kegg_pathway_desc !="NULL")|> 
  filter(kegg_pathway_desc !="") |> 
  mutate(kegg_pathway_desc = gsub("==.*$","", kegg_pathway_desc)) |> #Get rid of multiple pathways - keep first (imperfect)
  group_by(kegg_pathway_desc) |>
  summarise(across(2:37, ~ sum(.x !=0))) %>%
  mutate(across(2:37, ~ ./sum(.)))

#Get fraction of genes assigned to each phylum among transferred genes - I
genes.inoc.i.fun <- genes.inoc.rel.absent %>%
  mutate(gene_id = row.names(.)) %>%
  left_join(annot.inoc) %>%
  filter(kegg_pathway_desc !="NULL")|> 
  filter(kegg_pathway_desc !="") |> 
  mutate(kegg_pathway_desc = gsub("==.*$","", kegg_pathway_desc)) |> #Get rid of multiple pathways - keep first (imperfect)
  group_by(kegg_pathway_desc) |>
  summarise(across(2:39, ~ sum(.x >margin*dil*i.i))) %>%
  select(contains("kegg") | contains("I1")) %>%
  mutate(across(2:13, ~ ./sum(.)))

#Get fraction of genes assigned to each phylum among transferred genes - NI
genes.inoc.ni.fun <- genes.inoc.rel.absent %>%
  mutate(gene_id = row.names(.)) %>%
  left_join(annot.inoc) %>%
  filter(kegg_pathway_desc !="NULL")|> 
  filter(kegg_pathway_desc !="") |> 
  mutate(kegg_pathway_desc = gsub("==.*$","", kegg_pathway_desc)) |> #Get rid of multiple pathways - keep first (imperfect)
  group_by(kegg_pathway_desc) |>
  summarise(across(2:39, ~ sum(.x >margin*dil*i.ni))) %>%
  select(contains("kegg") | contains("I2")) %>%
  mutate(across(2:13, ~ ./sum(.)))

#Create long file with all samples for plotting
genes.fun.long <- left_join(genes.inoc.i.fun,genes.inoc.ni.fun) %>%
  left_join(select(genes.inoc.fun, i.i, i.ni, fun_phylum), by = join_by(fun_phylum), suffix = c(".tr",".inoc")) %>%
  left_join(select(genes.noninoc.fun, !contains("I3")), by = join_by(fun_phylum), suffix = c("",".noninoc")) %>%
  rowwise() |>
  mutate(TransMean = mean(c_across(2:24)), .before=T1.I1.1)|>
  ungroup() |>
  filter(TransMean>0.01) %>%
  bind_rows(data.frame(fun_phylum = "Others", t(1-colSums(.[,3:52])))) %>%
  gather(sample,RelAbund,3:52) %>% #transform in long format for ggplot
  select(-TransMean) %>%
  mutate(Dataset=c(rep("transferred",288), rep ("inoculum",24), rep("rhizo",288))) %>%
  mutate(Sample=gsub("(.noninoc)|(.tr)","", .$sample)) %>%
  select(-sample) %>%
  left_join(rownames_to_column(map,var="Sample"))

genes.fun.long$Dataset=factor(genes.fun.long$Dataset, c("rhizo", "inoculum", "transferred"))
stack.trans.fun <-  ggplot(genes.fun.long, aes(fill=fun_phylum, y=RelAbund, x=Dataset)) + 
  geom_bar( stat="identity", position="fill")+
  ylab("Relative abundance")+ 
  xlab("")+
  scale_fill_manual(values=palette(), guide=guide_legend(label.theme = element_text(face="italic", size=8)))+
  facet_grid(.~Inoculum, scales="free_x", space="free_x", labeller = as_labeller(c("15%"="LW", "50%"="HW", "Inoculum"="Inoculum", "i"="Intermittent", "ni"="Continuous", "ctrl"="Uninoculated")))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank())+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

stack.trans.fun

#Test for significance vs. rhizosphere
#Normality assumption
genes.fun.long %>% 
  filter(SoilWaterContent %in% c("15%","50%")) %>%
  group_by(fun_phylum) %>%
  shapiro_test(RelAbund) #All not ok.

#Equality of variances assumption
genes.fun.long %>% 
  filter(SoilWaterContent %in% c("15%","50%")) %>%
  group_by(fun_phylum) %>%
  levene_test(RelAbund~Inoculum*SoilWaterContent) #A few not OK.

#Computing the statistical test - paired t-test on effect of inoculation  
stat.test.fun <- genes.fun.long %>% 
  filter(SoilWaterContent %in% c("15%","50%")) %>%
  group_by(fun_phylum) %>%
  wilcox_test(RelAbund ~ Dataset, paired = TRUE)
stat.test.tax #Many significant - see output