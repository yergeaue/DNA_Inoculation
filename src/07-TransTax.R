###Make graphical representation and stat tests of taxonomy of transferred genes vs. all genes
genes.inoc.rel <- readRDS(here("data","intermediate", "genes.inoc.rel.RDS"))
annot.inoc <- readRDS(here("data","intermediate", "annot.inoc.RDS"))
genes.noninoc.rel <- readRDS(here("data","intermediate", "genes.noninoc.rel.RDS"))
annot.noninoc <- readRDS(file = here("data","intermediate", "annot.noninoc.RDS"))
genes.inoc.rel.absent <- readRDS(here("data", "intermediate", "genes.inoc.rel.absent.RDS"))
trans.bin.tax <- readRDS(file=here("data","intermediate", "trans.bin.tax.RDS"))
trans.bin.fun <- readRDS(file=here("data","intermediate", "trans.bin.fun.RDS"))
map <- readRDS(here("data", "intermediate", "map.RDS"))
palette(paletteer_d("ggthemes::Tableau_20"))

###Relative abundance in inoculum and rhizosphere vs. frequency in transferred genes (trans.bin.tax)
##Taxonomy
phylum.inoc <- genes.inoc.rel|>
  rownames_to_column(var = "gene_id") |>
  left_join(annot.inoc) |>
  filter(tax_phylum != "NULL") |>
  group_by(tax_phylum)|>
  summarise(across(2:39, ~ sum(.)))|>
  mutate(across(2:39, ~ ./sum(.)))

phylum.noninoc <- genes.noninoc.rel|>
  rownames_to_column(var = "gene_id") |>
  left_join(annot.noninoc) |>
  filter(tax_phylum != "NULL") |>
  select(!contains("I3"))|>
  group_by(tax_phylum)|>
  summarise(across(2:25, ~ sum(.)))|>
  mutate(across(2:25, ~ ./sum(.))) |>
  rowwise() |>
  mutate(HW = mean(c_across(2:13)), LW = mean(c_across(14:25)), .before = "T1.I1.1") |>
  ungroup()

#Combine all phylum tables together
trans.phylum.all <- trans.bin.tax |>
  left_join(phylum.noninoc, suffix = c("_trans","_rhizo"), by = join_by(tax_phylum)) |>
  left_join(select(phylum.inoc, tax_phylum, i.i, i.ni)) |>
  filter((i+ni+HW+LW/4)>0.05) |>
  select(!c(i,ni,HW,LW))|>
  pivot_longer(cols = -tax_phylum) |>
  pivot_wider(names_from = tax_phylum)|>
  rename("sample" = 1) |>
  rowwise() |>
  mutate(Others = 1-sum(c_across(2:10)))|>
  ungroup() |>
  pivot_longer(cols = c(2:11), names_to = "Phylum", values_to = "RelAbund")|>
  mutate(Dataset=c(rep("Transferred",240), rep ("Rhizosphere",240), rep("Inoculum",20))) |>
  mutate(Sample=gsub("(_trans)|(_rhizo)","", .data$sample)) |>
  select(-sample) |>
  left_join(rownames_to_column(map, var="Sample"))

saveRDS(trans.phylum.all, here("data", "intermediate", "trans.phylum.all"))

#phylum.trans.all$Dataset=factor(trans.phylum.all$Dataset, c("Transferred", "Intermittent", "Continuous", "LW", "HW"))
stack.trans.tax  <- 
  ggplot(trans.phylum.all, aes(fill=Phylum, y=RelAbund, x=Dataset)) + 
  geom_bar( stat="identity", position="fill")+
  ylab("Relative abundance")+  
  xlab("")+
  scale_fill_manual(values=palette(), guide=guide_legend(label.theme = element_text(face="italic", size=8)))+
  #scale_x_discrete(labels=c("Intermittent", "Continuous", "Uninoculated"),  )+
  facet_grid(.~Inoculum, scales="free_x", space="free_x", labeller = as_labeller(c("ni"="Continuous", "i"="Intermittent")))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

stack.trans.tax

#Test for significance vs. rhizosphere
#Normality assumption
trans.phylum.all %>% 
  filter(SoilWaterContent %in% c("15%","50%")) %>%
  group_by(Phylum) %>%
  shapiro_test(RelAbund) #All not ok.

#Equality of variances assumption
trans.phylum.all %>% 
  filter(SoilWaterContent %in% c("15%","50%")) %>%
  group_by(Phylum) %>%
  levene_test(RelAbund~Inoculum*SoilWaterContent) #A few not OK.

#Computing the statistical test - paired t-test on effect of inoculation  
stat.test.tax <- trans.phylum.all |> 
  filter(SoilWaterContent %in% c("15%","50%")) |>
  group_by(Phylum) |>
  arrange(Dataset,Sample) |>
  wilcox_test(RelAbund ~ Dataset, paired = TRUE)
stat.test.tax #Many significant - see output

#Means for text
trans.phylum.all |> 
  filter(SoilWaterContent %in% c("15%","50%")) |>
  group_by(Phylum, Dataset) |>
  summarise(100*mean(RelAbund))

##Functions
#Very little are defined in KEGG
#Just make a summary figure with all samples merged per inoculum
fun.inoc <- genes.inoc.rel|>
  rownames_to_column(var = "gene_id") |>
  left_join(annot.inoc) |>
  filter(kegg_pathway_desc !="NULL")|> 
  filter(kegg_pathway_desc !="") |> 
  mutate(kegg_pathway_desc = gsub("==.*$","", kegg_pathway_desc)) |> #Get rid of multiple pathways - keep first (imperfect)
  group_by(kegg_pathway_desc)|>
  summarise(across(2:39, ~ sum(.)))|>
  mutate(across(2:39, ~ ./sum(.)))

fun.noninoc <- genes.noninoc.rel|>
  rownames_to_column(var = "gene_id") |>
  left_join(annot.noninoc) |>
  filter(kegg_pathway_desc !="NULL")|>
  filter(kegg_pathway_desc !="") |> 
  mutate(kegg_pathway_desc = gsub("==.*$","", kegg_pathway_desc)) |> #Get rid of multiple pathways - keep first (imperfect)
  select(!contains("I3"))|>
  group_by(kegg_pathway_desc)|>
  summarise(across(2:25, ~ sum(.)))|>
  mutate(across(2:25, ~ ./sum(.))) |>
  rowwise() |>
  mutate(HW = mean(c_across(2:13)), LW = mean(c_across(14:25)), .before = "T1.I1.1") |>
  ungroup()

#Combine all function tables together
trans.fun.all <- trans.bin.fun |>
  left_join(fun.noninoc, suffix = c("_trans","_rhizo"), by = join_by(kegg_pathway_desc)) |>
  left_join(select(fun.inoc, kegg_pathway_desc, i.i, i.ni), by = join_by(kegg_pathway_desc), suffix = c("_TRANS","")) |>
  filter((i.ni_TRANS+i.i_TRANS+HW+LW/4)>0.1) |>
  select(!contains(c("HW","LW","_trans"), ignore.case = FALSE))|>
  pivot_longer(cols = -kegg_pathway_desc) |>
  pivot_wider(names_from = kegg_pathway_desc)|>
  rename("sample" = 1) |>
  rowwise() |>
  mutate(Others = 1-sum(c_across(2:12)))|>
  ungroup() |>
  pivot_longer(cols = c(2:13), names_to = "kegg_pathway_desc", values_to = "RelAbund")|>
  mutate(Dataset=c(rep("Transferred",24), rep ("Rhizosphere",288), rep("Inoculum",24))) |>
  mutate(Sample=gsub("(_TRANS)|(_rhizo)","", .data$sample)) |>
  select(-sample) |>
  left_join(rownames_to_column(map, var="Sample"))

saveRDS(trans.fun.all, here("data", "intermediate", "trans.fun.all"))

stack.trans.fun  <- 
  ggplot(trans.fun.all, aes(fill=kegg_pathway_desc, y=RelAbund, x=Dataset)) + 
  geom_bar( stat="identity", position="fill")+
  ylab("Relative abundance")+  
  xlab("")+
  scale_fill_manual(values=palette(), guide=guide_legend(label.theme = element_text( size=8)))+
  #scale_x_discrete(labels=c("Intermittent", "Continuous", "Uninoculated"),  )+
  facet_grid(.~Inoculum, scales="free_x", space="free_x", labeller = as_labeller(c("ni"="Continuous", "i"="Intermittent")))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

stack.trans.fun
