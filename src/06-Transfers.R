###Find genes potentially transferred
##Find genes
#Load data: all samples mapped on assembly of inoculum
genes.inoc.rel <- readRDS(here("data","intermediate", "genes.inoc.rel.RDS"))
map <- readRDS(here("data", "intermediate", "map.RDS"))
map.nosoil <- map[map$SoilWaterContent != "Soil",]
genes.inoc.rel <- genes.inoc.rel[,order(colnames(genes.inoc.rel))]
sum(rownames(map.nosoil)==colnames(genes.inoc.rel)) #38
annot.inoc <- readRDS(here("data", "intermediate", "annot.inoc.RDS"))
lgt.genes.pub <- readRDS(here("data","intermediate","lgt.genes.pub.RDS"))

#Select genes that are 
#1)Absent in controls
genes.inoc.rel.absent <- genes.inoc.rel[rowSums(genes.inoc.rel[,map.nosoil$Inoculum=="ctrl"])==0,  ] #subset to absent in non inoculated
#2)Present in rhizosphere
genes.inoc.rel.absent <- genes.inoc.rel.absent[rowSums(genes.inoc.rel.absent[,c(3:14,21:32)])!=0,] #subset to present in the rhizosphere
genes.inoc.rel.absent <- rownames_to_column(genes.inoc.rel.absent, var = "gene_id")
#3) present in inoculum 
#Given since mapped on inoculum assembly. Careful, some are only present in one of the two inoculum
#See below for solution 

#4)Abundant at a safety margin (x) above relics
dil <- 850/(50*100/0.5*1500/1000) #Dilution factor - added 850 ug in 15,000ug
margin <- 128

#Create binary tables if gene is above threshold. Separate for the two inoculum, as threshold is related to abundance in inoculum
ni.reps <- genes.inoc.rel.absent |>
  select(contains(c("gene_id", "I2","i.ni"))) |>
  mutate(across(2:13, ~ case_when(.>margin*dil*.data$i.ni ~ 1, TRUE ~ 0))) |>
  mutate(across(2:13, ~ case_when(i.ni==0 ~ 0, .default = . ))) #Remove cases when the gene is absent in inoc
sum(colSums(ni.reps[,2:13])) #332

i.reps <- genes.inoc.rel.absent |>
  select(contains(c("gene_id","I1","i.i"))) |>
  mutate(across(2:13, ~ case_when(.>margin*dil*.data$i.i ~ 1, TRUE ~ 0))) |>
  mutate(across(2:13, ~ case_when(i.i==0 ~ 0, .default = . ))) #Remove cases when the gene is absent in inoc
sum(colSums(i.reps[,2:13])) #362

#Join the two tables
trans.bin <- left_join(ni.reps, i.reps)

#Add annotations
trans.bin.annot <- trans.bin |>
  left_join(select(annot.inoc, gene_id, kegg_entry, kegg_definition, kegg_module_desc, 
                    kegg_pathway_desc, tax_phylum, tax_family, tax_genus)) #Add annotations (selected var)

#Create phylum level taxonomy table
trans.bin.tax <- trans.bin.annot |>
  filter(tax_phylum!="NULL")|> 
  group_by(tax_phylum) |>
  summarise(across(c(2:13,15:26), ~sum(.))) |>
  mutate(across(2:25, ~ ./sum(.)))|>
  rowwise() |>
  mutate(ni = mean(c_across(2:13)), i = mean(c_across(14:25)), .before = "T1.I2.1") |>
  ungroup()

saveRDS(trans.bin.tax, here("data","intermediate", "trans.bin.tax.RDS"))

#Create kegg_pathway table - makes little sense, too few left (39 positives)
trans.bin.fun <- trans.bin.annot|>
  filter(kegg_pathway_desc !="NULL")|> #1550 genes at this stage
  filter(kegg_pathway_desc !="") |> #1270 genes at this stage
  mutate(kegg_pathway_desc = gsub("==.*$","", kegg_pathway_desc)) |> #Get rid of multiple pathways - keep first (imperfect)
  group_by(kegg_pathway_desc) |>
  summarise(across(c(2:13,15:26), ~sum(.))) |>
  mutate(across(2:25, ~ ./sum(.))) 


#Create summary table for publication
trans.pub <- trans.bin.annot |>
  filter(kegg_entry != "NULL")|>
  rowwise() |>
  mutate(count.i = sum(c_across(15:26))) |>
  mutate(count.ni = sum(c_across(2:13))) |>
  ungroup()|>
  filter(count.i+count.ni>0)|>
  mutate(kegg_pathway_desc = gsub("==.*$","", kegg_pathway_desc)) |> #Get rid of multiple pathways - keep first (imperfect)
  select("Gene ID"=gene_id, "Count intermittent"=count.i, "Count continuous"=count.ni, "KEGG entry"=kegg_entry, 
         "KEGG definition"=kegg_definition, "KEGG module"=kegg_module_desc, 
         "KEGG pathway"=kegg_pathway_desc, "Phylum"=tax_phylum, "Family"=tax_family, "Genus"=tax_genus )


#Overlap with lgt Waafle at the KEGG entry level ? - 3 entries
trans.lgt.overlap <- lgt.genes.pub |>
  filter(kegg_entry!="NULL")|>
  filter(kegg_entry %in% trans.pub$`KEGG entry`)|>
  left_join(trans.pub, join_by(kegg_entry==`KEGG entry`))
