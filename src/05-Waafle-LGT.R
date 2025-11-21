##Waafle method to find lgt
#Ran the program according to instructions on the cluster. Resulted in 111 potential LGT contigs
lgt <- readRDS(here("data","intermediate", "lgt.RDS"))
blast_output <- readRDS(here("data", "intermediate", "blast_output.RDS")) #To link with inoculum
genes.noninoc.rel <- readRDS(here("data","intermediate", "genes.noninoc.rel.RDS")) #Genes relative abundance (non-inoculum assembly)
gff.noninoc <- readRDS(here("data", "intermediate", "gff.noninoc.RDS")) #Contig information for the non-inoculum only assembly
annot.noninoc <- readRDS(file = here("data","intermediate", "annot.noninoc.RDS"))
map <- readRDS(here("data", "intermediate", "map.RDS"))
palette(paletteer_d("ggthemes::Tableau_20"))

#Fix the lgt file, missing genes for some contigs
lgt.summary <-  lgt |>
  left_join(gff.noninoc, by = c("CONTIG_NAME"="contig_id"), keep =TRUE) |>
  group_by(contig_id) |>
  summarise(nbgenes=n())
lgt.summary$contig_id==lgt$CONTIG_NAME #Check order to make sure - all true

#Add "~" at the end of synteny that are too short - need two passes as some miss two "~"
lgt.fix <- lgt |>
  mutate(SYNTENY_2=case_when(
    CONTIG_NAME%in%lgt[lgt.summary$nbgenes!=lengths(str_split(SYNTENY, pattern = "")), "CONTIG_NAME"]
    ~ paste(SYNTENY, "~", sep=""),
    TRUE ~ paste(SYNTENY))
  ) |>
  mutate(SYNTENY_3=case_when(
    CONTIG_NAME%in%lgt[lgt.summary$nbgenes!=lengths(str_split(SYNTENY_2, pattern = "")), "CONTIG_NAME"]
    ~ paste(SYNTENY_2, "~", sep=""),
    TRUE ~ paste(SYNTENY_2))
  )

#check if fixed - should be 460
sum(lengths(str_split(lgt.fix$SYNTENY_3, "")))

#Plot
lgt.plot <-  lgt.fix |>
  left_join(gff.noninoc, by = c("CONTIG_NAME"="contig_id"), keep =TRUE) |>
  mutate(CLADE=unlist(str_split(lgt.fix$SYNTENY_3,""))) |>
  mutate(CLADE_2=case_when(
    CLADE=="A" ~ CLADE_A,
    CLADE=="B" ~ CLADE_B,
    TRUE ~ "NULL"
  ))
  
ggplot(lgt.plot, aes(xmin = start, xmax = end, y = contig_id, fill = CLADE, label=gene_id)) +
  geom_gene_arrow(aes(forward = strand == "+")) +
  geom_gene_label(align = "left")+
  #facet_wrap(~ contig_id, scales = "free", ncol = 2) +
  scale_fill_manual(values = palette()) 

saveRDS(lgt.plot, here("data","intermediate","lgt.plot.RDS"))

##LGT genes
lgt.genes <- lgt.plot |>
  filter(CLADE=="B") |>
  select(gene_id) |>
  left_join(annot.noninoc)

#Table for publication
lgt.genes.pub <- lgt.genes |>
  select(gene_id, product_name, kegg_entry, kegg_definition, kegg_module_desc, kegg_pathway_desc, tax_phylum, tax_order, tax_genus)|>
  mutate(kegg_pathway_desc = gsub("==.*$","", kegg_pathway_desc)) #Get rid of multiple pathways - keep first (imperfect)

saveRDS(lgt.genes.pub, here("data","intermediate","lgt.genes.pub.RDS"))

#Summary stats for text
lgt.genes.summary <- lgt.genes |>
  group_by(tax_phylum) |>
  summarise(count = n()) |>
  mutate(percent = count/sum(count)*100)

lgt.genes.summary.2 <- lgt.genes |>
  group_by(kegg_pathway_desc) |>
  summarise(count = n()) |>
  mutate(percent = count/sum(count)*100)

#Which ones are found in inoculum?-none!
#Find genes from Waafle that have match in blast file
lgt.genes.blast <- lgt.genes |> 
  left_join(blast_output, by = c("gene_id"="sseqid"), keep = TRUE)
sum(!is.na(lgt.genes.blast$qseqid)) #0

##Look for distribution of genes across treatments
lgt.genes.rel <-  lgt.genes |>
  left_join(mutate(genes.noninoc.rel, gene_id=row.names(genes.noninoc.rel)))

lgt.genes.anova <- lgt.genes.rel[,c(1,39:74)] |>
  column_to_rownames(var = "gene_id") |>
  t() |> 
  data.frame() |>
  rownames_to_column(var = "Sample")|>
  left_join(rownames_to_column(map, var = "Sample")) |>
  pivot_longer(cols = c(2:115), names_to = "gene_id", values_to = "RelAbund")
  
#Test for significance
#Normality assumption
lgt.genes.anova |> 
  group_by(gene_id) |>
  shapiro_test(RelAbund) #Many not ok

#Equality of variances assumption
lgt.genes.anova |> 
  group_by(gene_id) |>
  levene_test(RelAbund~Inoculum*SoilWaterContent) 

#Computing the statistical tests - Paired Wilcoxon
stat.test.lgt.SWC <- lgt.genes.anova |>
  group_by(gene_id) |>
  wilcox_test(RelAbund~SoilWaterContent, paired = TRUE)|>
  filter(p<0.05) |>
  left_join(lgt.genes.pub)
stat.test.lgt.SWC #28 significant
summary.stat.lgt.SWC <- lgt.genes.anova |>
  mutate(RelAbund=RelAbund*1e6)|>
  group_by(gene_id, SoilWaterContent) |>
  get_summary_stats(RelAbund, show = "mean")|>
  filter(gene_id%in%stat.test.lgt.SWC$gene_id)|>
  left_join(lgt.genes.pub)
summary.stat.lgt.SWC #14 are more abundant in LW
stat.test.lgt.inoc <- lgt.genes.anova |>
  group_by(gene_id) |>
  wilcox_test(RelAbund~Inoculum, paired = TRUE) |>
  filter(p.adj<0.05)|>
  left_join(lgt.genes.pub)
stat.test.lgt.inoc #2 significant, only one with ctrl
summary.stat.lgt.inoc <- lgt.genes.anova |>
  mutate(RelAbund=RelAbund*1e6)|>
  group_by(gene_id, Inoculum) |>
  get_summary_stats(RelAbund, show = "mean")|>
  filter(gene_id%in%stat.test.lgt.inoc$gene_id)|>
  left_join(lgt.genes.pub)
summary.stat.lgt.inoc #groES is more abundant in the intermittent inoculum
