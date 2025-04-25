###Find genes potentially transferred
##Find genes
#Load data: all samples mapped on assembly of inoculum
genes.inoc.rel <- readRDS(here("data","intermediate", "genes.inoc.rel.RDS"))
map <- readRDS(here("data", "intermediate", "map.RDS"))
map.nosoil <- map[map$SoilWaterContent != "Soil",]
genes.inoc.rel <- genes.inoc.rel[,order(colnames(genes.inoc.rel))]
sum(rownames(map.nosoil)==colnames(genes.inoc.rel)) #38
annot.inoc <- readRDS(here("data", "intermediate", "annot.inoc.RDS"))

#Select genes that are 
#1) present in inoculum (given, since assembly if from inoculum only)
#2)Absent in controls
genes.inoc.rel.absent <- genes.inoc.rel[rowSums(genes.inoc.rel[,map.nosoil$Inoculum=="ctrl"])==0,  ] #subset to absent in non inoculated
#3)Present in rhizosphere
genes.inoc.rel.absent <- genes.inoc.rel.absent[rowSums(genes.inoc.rel.absent[,c(3:14,21:32)])!=0,] #subset to present in the rhizosphere
#Keep only genes that are defined at the kegg_entry
annot.inoc.def <- annot.inoc[annot.inoc$kegg_entry != "NULL", ]
genes.inoc.rel.absent <- genes.inoc.rel.absent[row.names(genes.inoc.rel.absent) %in% row.names(annot.inoc.def),]
genes.inoc.rel.absent$gene_id <- row.names(genes.inoc.rel.absent)
saveRDS(genes.inoc.rel.absent, here("data", "intermediate","genes.inoc.rel.absent.RDS"))

#4)Abundant at a safety margin (x) above relics
dil <- 850/(50*100/0.5*1500/1000) #Dilution factor - added 850 ug in 15,000ug
margin <- 10

#Irrigated inoculum, 15% SWHC
i.15.table <- genes.inoc.rel.absent %>%
  select(contains("T1.I1")) %>%
  mutate(count=rowSums(.>margin*dil*genes.inoc.rel.absent$i.i)) %>%
  mutate(gene_id=rownames(genes.inoc.rel.absent)) %>%
  filter(count>1) %>%
  left_join(select(annot.inoc.def, gene_id, kegg_entry, kegg_definition, tax_phylum, tax_order, tax_genus)) %>% #Add annotations (selected var)
  arrange(desc(count)) %>%
  select("Gene ID"=gene_id, "Count"=count, "KEGG entry"=kegg_entry, 
         "KEGG definition"=kegg_definition, "Phylum"=tax_phylum, "Order"=tax_order, "Genus"=tax_genus )
  
saveRDS(i.15.table, here("data","intermediate","i.15.table.RDS"))

#Look for taxonomical makeup
tax.i.15.sum <- i.15.table %>%
  group_by(Phylum) %>%
  summarise(sum = n()) %>%
  mutate(per =  round(100 *sum/sum(sum),1))

#Ambient inoculum, 15% SWHC
ni.15.table <- genes.inoc.rel.absent %>%
  select(contains("T1.I2")) %>%
  mutate(count=rowSums(.>margin*dil*genes.inoc.rel.absent$i.ni)) %>%
  mutate(gene_id=rownames(genes.inoc.rel.absent)) %>%
  filter(count>1) %>%
  left_join(select(annot.inoc.def, gene_id, kegg_entry, kegg_definition, tax_phylum, tax_order, tax_genus)) %>% #Add annotations (selected var)
  arrange(desc(count)) %>%
  select("Gene ID"=gene_id, "Count"=count, "KEGG entry"=kegg_entry, 
         "KEGG definition"=kegg_definition, "Phylum"=tax_phylum, "Order"=tax_order, "Genus"=tax_genus )

saveRDS(ni.15.table, here("data","intermediate","ni.15.table.RDS"))

#Look for taxonomical makeup
tax.ni.15.sum <- ni.15.table %>%
  group_by(Phylum) %>%
  summarise(sum = n()) %>%
  mutate(per =  round(100 *sum/sum(sum),1))

#Irrigated inoculum, 50% SWHC
i.50.table <- genes.inoc.rel.absent %>%
  select(contains("T2.I1")) %>%
  mutate(count=rowSums(.>margin*dil*genes.inoc.rel.absent$i.i)) %>%
  mutate(gene_id=rownames(genes.inoc.rel.absent)) %>%
  filter(count>1) %>%
  left_join(select(annot.inoc.def, gene_id, kegg_entry, kegg_definition, tax_phylum, tax_order, tax_genus)) %>% #Add annotations (selected var)
  arrange(desc(count)) %>%
  select("Gene ID"=gene_id, "Count"=count, "KEGG entry"=kegg_entry, 
         "KEGG definition"=kegg_definition, "Phylum"=tax_phylum, "Order"=tax_order, "Genus"=tax_genus )

saveRDS(i.50.table, here("data","intermediate","i.50.table.RDS"))

#Look for taxonomical makeup
tax.i.50.sum <- i.50.table %>%
  group_by(Phylum) %>%
  summarise(sum = n()) %>%
  mutate(per =  round(100 *sum/sum(sum),1))

#Ambient inoculum, 50% SWHC
ni.50.table <- genes.inoc.rel.absent %>%
  select(contains("T2.I2")) %>%
  mutate(count=rowSums(.>margin*dil*genes.inoc.rel.absent$i.ni)) %>%
  mutate(gene_id=rownames(genes.inoc.rel.absent)) %>%
  filter(count>1) %>%
  left_join(select(annot.inoc.def, gene_id, kegg_entry, kegg_definition, tax_phylum, tax_order, tax_genus)) %>% #Add annotations (selected var)
  arrange(desc(count)) %>%
  select("Gene ID"=gene_id, "Count"=count, "KEGG entry"=kegg_entry, 
         "KEGG definition"=kegg_definition, "Phylum"=tax_phylum, "Order"=tax_order, "Genus"=tax_genus )

saveRDS(ni.50.table, here("data","intermediate","ni.50.table.RDS"))

#Look for taxonomical makeup
tax.ni.50.sum <- ni.50.table %>%
  group_by(Phylum) %>%
  summarise(sum = n()) %>%
  mutate(per =  round(100 *sum/sum(sum),1))

##Intersect
intersect(intersect(intersect(i.15.table$`Gene ID`, i.50.table$`Gene ID`),ni.15.table$`Gene ID`),ni.50.table$`Gene ID`)
#"gene_id_69099" "gene_id_95425"  "gene_id_179887" "gene_id_196132" "gene_id_215568" "gene_id_226170" "gene_id_233774" "gene_id_268173"
length(intersect(i.15.table$`Gene ID`, i.50.table$`Gene ID`)) #47
length(intersect(i.15.table$`Gene ID`, ni.15.table$`Gene ID`)) #33
length(intersect(i.15.table$`Gene ID`, ni.50.table$`Gene ID`)) #36
length(intersect(i.50.table$`Gene ID`, ni.15.table$`Gene ID`)) #34
length(intersect(i.50.table$`Gene ID`, ni.50.table$`Gene ID`)) #38
length(intersect(ni.15.table$`Gene ID`, ni.50.table$`Gene ID`)) #47

##Use all genes that are above threshold, not only the ones appearing twice
#Irrigated inoculum, 15% SWHC
i.15.table.0 <- genes.inoc.rel.absent %>%
  select(contains("T1.I1")) %>%
  mutate(count=rowSums(.>margin*dil*genes.inoc.rel.absent$i.i)) %>%
  mutate(gene_id=rownames(genes.inoc.rel.absent)) %>%
  filter(count>0) %>%
  left_join(select(annot.inoc.def, gene_id, kegg_entry, kegg_definition, tax_phylum, tax_order, tax_genus)) %>% #Add annotations (selected var)
  arrange(desc(count)) %>%
  select("Gene ID"=gene_id, "Count"=count, "KEGG entry"=kegg_entry, 
         "KEGG definition"=kegg_definition, "Phylum"=tax_phylum, "Order"=tax_order, "Genus"=tax_genus )

saveRDS(i.15.table.0, here("data","intermediate","i.15.table.0.RDS"))

#Ambient inoculum, 15% SWHC
ni.15.table.0 <- genes.inoc.rel.absent %>%
  select(contains("T1.I2")) %>%
  mutate(count=rowSums(.>margin*dil*genes.inoc.rel.absent$i.ni)) %>%
  mutate(gene_id=rownames(genes.inoc.rel.absent)) %>%
  filter(count>0) %>%
  left_join(select(annot.inoc.def, gene_id, kegg_entry, kegg_definition, tax_phylum, tax_order, tax_genus)) %>% #Add annotations (selected var)
  arrange(desc(count)) %>%
  select("Gene ID"=gene_id, "Count"=count, "KEGG entry"=kegg_entry, 
         "KEGG definition"=kegg_definition, "Phylum"=tax_phylum, "Order"=tax_order, "Genus"=tax_genus )

saveRDS(ni.15.table.0, here("data","intermediate","ni.15.table.0.RDS"))

#Irrigated inoculum, 50% SWHC
i.50.table.0 <- genes.inoc.rel.absent %>%
  select(contains("T2.I1")) %>%
  mutate(count=rowSums(.>margin*dil*genes.inoc.rel.absent$i.i)) %>%
  mutate(gene_id=rownames(genes.inoc.rel.absent)) %>%
  filter(count>0) %>%
  left_join(select(annot.inoc.def, gene_id, kegg_entry, kegg_definition, tax_phylum, tax_order, tax_genus)) %>% #Add annotations (selected var)
  arrange(desc(count)) %>%
  select("Gene ID"=gene_id, "Count"=count, "KEGG entry"=kegg_entry, 
         "KEGG definition"=kegg_definition, "Phylum"=tax_phylum, "Order"=tax_order, "Genus"=tax_genus )

saveRDS(i.50.table.0, here("data","intermediate","i.50.table.0.RDS"))

#Ambient inoculum, 50% SWHC
ni.50.table.0 <- genes.inoc.rel.absent %>%
  select(contains("T2.I2")) %>%
  mutate(count=rowSums(.>margin*dil*genes.inoc.rel.absent$i.ni)) %>%
  mutate(gene_id=rownames(genes.inoc.rel.absent)) %>%
  filter(count>0) %>%
  left_join(select(annot.inoc.def, gene_id, kegg_entry, kegg_definition, tax_phylum, tax_order, tax_genus)) %>% #Add annotations (selected var)
  arrange(desc(count)) %>%
  select("Gene ID"=gene_id, "Count"=count, "KEGG entry"=kegg_entry, 
         "KEGG definition"=kegg_definition, "Phylum"=tax_phylum, "Order"=tax_order, "Genus"=tax_genus )

saveRDS(ni.50.table.0, here("data","intermediate","ni.50.table.0.RDS"))

#Combine all samples
tax.trans <- bind_rows(list(i.15=i.15.table.0,ni.15=ni.15.table.0,i.50=i.50.table.0,ni.50=ni.50.table.0), .id="treatment")
length(unique(tax.trans$`Gene ID`))#first see amount of unique:1169
#Create phylum level taxonomy table
tax.trans.sum <- tax.trans %>%
  distinct(`Gene ID`, .keep_all = TRUE) %>%
  group_by(Phylum) %>%
  summarise(sum = n()) %>%
  mutate(per =  round(100 *sum/sum(sum),1))



