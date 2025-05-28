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
#1)Absent in controls
genes.inoc.rel.absent <- genes.inoc.rel[rowSums(genes.inoc.rel[,map.nosoil$Inoculum=="ctrl"])==0,  ] #subset to absent in non inoculated
#2)Present in rhizosphere
genes.inoc.rel.absent <- genes.inoc.rel.absent[rowSums(genes.inoc.rel.absent[,c(3:14,21:32)])!=0,] #subset to present in the rhizosphere
#3) present in inoculum 
genes.inoc.rel.absent.i <- genes.inoc.rel.absent[genes.inoc.rel.absent$i.i !=0,]
genes.inoc.rel.absent.ni <- genes.inoc.rel.absent[genes.inoc.rel.absent$i.ni !=0,]

#NOT USED ANYMORE: Keep only genes that are defined at the kegg_entry
#annot.inoc.def <- annot.inoc[annot.inoc$kegg_entry != "NULL", ]
#genes.inoc.rel.absent <- genes.inoc.rel.absent[row.names(genes.inoc.rel.absent) %in% row.names(annot.inoc.def),]

genes.inoc.rel.absent.i$gene_id <- row.names(genes.inoc.rel.absent.i)
genes.inoc.rel.absent.ni$gene_id <- row.names(genes.inoc.rel.absent.ni)
genes.inoc.rel.absent$gene_id <- row.names(genes.inoc.rel.absent)

saveRDS(genes.inoc.rel.absent, here("data", "intermediate","genes.inoc.rel.absent.RDS"))

#4)Abundant at a safety margin (x) above relics
dil <- 850/(50*100/0.5*1500/1000) #Dilution factor - added 850 ug in 15,000ug
margin <- 128

##Use all genes that are above threshold, not only the ones appearing twice
#Irrigated inoculum, 15% SWHC
i.15.table.0 <- genes.inoc.rel.absent %>%
  select(contains("T1.I1")) %>%
  mutate(count=rowSums(.>margin*dil*genes.inoc.rel.absent.i$i.i)) %>%
  rownames_to_column(var = "gene_id") %>%
  filter(count>0) %>%
  left_join(select(annot.inoc, gene_id, kegg_entry, kegg_definition, tax_phylum, tax_family, tax_genus)) %>% #Add annotations (selected var)
  arrange(desc(count)) %>%
  select("Gene ID"=gene_id, "Count"=count, "KEGG entry"=kegg_entry, 
         "KEGG definition"=kegg_definition, "Phylum"=tax_phylum, "Family"=tax_family, "Genus"=tax_genus )

#i.15.table.0.def <- i.15.table.0[i.15.table.0$`KEGG entry`!="NULL" | i.15.table.0$Phylum!="NULL",]

#saveRDS(i.15.table.0.def, here("data","intermediate","i.15.table.0.def.RDS"))

#Ambient inoculum, 15% SWHC
ni.15.table.0 <- genes.inoc.rel.absent %>%
  select(contains("T1.I2")) %>%
  mutate(count=rowSums(.>margin*dil*genes.inoc.rel.absent.ni$i.ni)) %>%
  rownames_to_column(var = "gene_id") %>%
  filter(count>0) %>%
  left_join(select(annot.inoc, gene_id, kegg_entry, kegg_definition, tax_phylum, tax_family, tax_genus)) %>% #Add annotations (selected var)
  arrange(desc(count)) %>%
  select("Gene ID"=gene_id, "Count"=count, "KEGG entry"=kegg_entry, 
         "KEGG definition"=kegg_definition, "Phylum"=tax_phylum, "Family"=tax_family, "Genus"=tax_genus )

#ni.15.table.0.def <- ni.15.table.0[ni.15.table.0$`KEGG entry`!="NULL" | ni.15.table.0$Phylum!="NULL",]

#saveRDS(ni.15.table.0.def, here("data","intermediate","ni.15.table.0.def.RDS"))

#Irrigated inoculum, 50% SWHC
i.50.table.0 <- genes.inoc.rel.absent %>%
  select(contains("T2.I1")) %>%
  mutate(count=rowSums(.>margin*dil*genes.inoc.rel.absent$i.i)) %>%
  rownames_to_column(var = "gene_id") %>%
  filter(count>0) %>%
  left_join(select(annot.inoc, gene_id, kegg_entry, kegg_definition, tax_phylum, tax_family, tax_genus)) %>% #Add annotations (selected var)
  arrange(desc(count)) %>%
  select("Gene ID"=gene_id, "Count"=count, "KEGG entry"=kegg_entry, 
         "KEGG definition"=kegg_definition, "Phylum"=tax_phylum, "Family"=tax_family, "Genus"=tax_genus )

#i.50.table.0.def <- i.50.table.0[i.50.table.0$`KEGG entry`!="NULL" | i.50.table.0$Phylum!="NULL",]

#saveRDS(i.50.table.0.def, here("data","intermediate","i.50.table.0.def.RDS"))

#Ambient inoculum, 50% SWHC
ni.50.table.0 <- genes.inoc.rel.absent %>%
  select(contains("T2.I2")) %>%
  mutate(count=rowSums(.>margin*dil*genes.inoc.rel.absent$i.ni)) %>%
  rownames_to_column(var = "gene_id") %>%
  filter(count>0) %>%
  left_join(select(annot.inoc, gene_id, kegg_entry, kegg_definition, tax_phylum, tax_family, tax_genus)) %>% #Add annotations (selected var)
  arrange(desc(count)) %>%
  select("Gene ID"=gene_id, "Count"=count, "KEGG entry"=kegg_entry, 
         "KEGG definition"=kegg_definition, "Phylum"=tax_phylum, "Family"=tax_family, "Genus"=tax_genus )

#ni.50.table.0.def <- ni.50.table.0[ni.50.table.0$`KEGG entry`!="NULL" | ni.50.table.0$Phylum!="NULL",]

#saveRDS(ni.50.table.0.def, here("data","intermediate","ni.50.table.0.def.RDS"))

##Intersect
intersect(intersect(intersect(i.15.table.0$`Gene ID`, i.50.table.0$`Gene ID`),ni.15.table.0$`Gene ID`),ni.50.table.0$`Gene ID`)
#0 genes
length(intersect(i.15.table.0$`Gene ID`, i.50.table.0$`Gene ID`)) #6
length(intersect(i.15.table.0$`Gene ID`, ni.15.table.0$`Gene ID`)) #3
length(intersect(i.15.table.0$`Gene ID`, ni.50.table.0$`Gene ID`)) #12
length(intersect(i.50.table.0$`Gene ID`, ni.15.table.0$`Gene ID`)) #4
length(intersect(i.50.table.0$`Gene ID`, ni.50.table.0$`Gene ID`)) #23
length(intersect(ni.15.table.0$`Gene ID`, ni.50.table.0$`Gene ID`)) #7

#Combine all samples
tax.trans.0 <- bind_rows(list(i.15=i.15.table.0,ni.15=ni.15.table.0,i.50=i.50.table.0,ni.50=ni.50.table.0), .id="treatment")
length(unique(tax.trans.0$`Gene ID`))#first see amount of unique:1215
sum(tax.trans.0$Count==6)#How many count=6 : 0
sum(tax.trans.0$Count==5)#How many count=5 : 0
sum(tax.trans.0$Count==4)#How many count=4 : 3
sum(tax.trans.0$Count==3)#How many count=3 : 41
sum(tax.trans.0$Count==2)#How many count=2 : 174
sum(tax.trans.0$Count==1)#How many count=1 : 1051

saveRDS(tax.trans.0, here("data", "intermediate", "tax.trans.0.RDS"))

#Create publication table - only defined
tax.trans.0.def <- tax.trans.0 |>
  filter(`KEGG entry`!="NULL")


#Create phylum level taxonomy table
tax.trans.0.sum <- tax.trans.0 |>
  #filter(Phylum!="NULL")|> #1067 genes at this stage
  group_by(Phylum) |>
  summarise(sum = n()) |>
  mutate(per =  round(100 *sum/sum(sum),1))


#Overlap with lgt Waafle at the KEGG entry level - 3 entries
trans.lgt.overlap <- lgt.genes.pub |>
  filter(kegg_entry!="NULL")|>
  filter(kegg_entry %in% tax.trans.0$`KEGG entry`)|>
  left_join(tax.trans.0, join_by(kegg_entry==`KEGG entry`))
