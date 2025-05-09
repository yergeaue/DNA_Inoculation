##Waafle method to find lgt
lgt <- readRDS(here("data","intermediate", "lgt.RDS"))
blast_output <- readRDS(here("data", "intermediate", "blast_output.RDS")) #To link with inoculum
genes.inoc.rel <- genes.inoc.rel[,order(colnames(genes.inoc.rel))] #Contigs made from the inoculum
gff.inoc <- readRDS(here("data", "intermediate", "gff.inoc.RDS")) #Contig information for the inoculum only assembly
gff.noninoc <- readRDS(here("data", "intermediate", "gff.noninoc.RDS")) #Contig information for the non-inoculum only assembly
annot.inoc <- readRDS(here("data", "intermediate", "annot.inoc.RDS"))
annot.noninoc <- readRDS(file = here("data","intermediate", "annot.noninoc.RDS"))

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
  scale_fill_manual(values = color6[c(2,4,6)]) 


#Are some lgt contigs absent in controls? 
contigs.lgt.absent <- genes.inoc.rel[rowSums(genes.inoc.rel[,map.nosoil$Inoculum=="ctrl"])==0,  ]

#LGT contigs (non-inoculum assembly) also found in inoculum assembly by blast
lgt.inoc <- blast_output |> 
  left_join(gff.noninoc, by = c("sseqid"="gene_id")) |> #add contig info for subject genes
  filter(contig_id%in%lgt.fix$CONTIG_NAME) |> #Keep only genes potentially transferred
  left_join(gff.inoc, by = c("qseqid"="gene_id"), suffix = c("_s", "_q")) |> #add contig info for query genes
  mutate(contig.q.vs.s = interaction(contig_id_q,contig_id_s, sep = " vs ")) |> #Add column linking the two contigs
  left_join(select(annot.noninoc, gene_id, product_name, kegg_entry, kegg_definition, tax_phylum, tax_genus), by = c("sseqid"="gene_id")) |>
  left_join(select(annot.inoc, gene_id, product_name, kegg_entry, kegg_definition, tax_phylum, tax_genus), 
            by = c("qseqid"="gene_id"), suffix = c("_s", "_q"))



gene.i.15.inoc <- gff.inoc |>
  filter(contig_id %in% blast.i.15$contig_id_q) |> #Keep only contigs from blast file
  left_join(select(annot.inoc, gene_id, product_name, kegg_entry, kegg_definition, tax_phylum)) |> #Add annotations (selected var)
  left_join(select(blast.i.15, contig_id_q, contig.q.vs.s), by = c("contig_id"="contig_id_q")) |> #Add contig linkage info
  mutate(focal = gene_id%in%blast.i.15$qseqid) |>
  mutate(assembly = "inoculum")

gene.i.15.noninoc <- gff.noninoc |>
  filter(contig_id %in% blast.i.15$contig_id_s) |> #Keep only contigs from blast file
  left_join(select(annot.noninoc, gene_id, product_name, kegg_entry, kegg_definition, tax_phylum)) |> #Add annotations (selected var)
  left_join(select(blast.i.15, contig_id_s, contig.q.vs.s), by = c("contig_id"="contig_id_s")) |> #Add contig linkage info
  mutate(focal = gene_id%in%blast.i.15$sseqid) |>
  mutate(assembly = "non-inoculum")

gene.i.15.both <- rbind(gene.i.15.inoc, gene.i.15.noninoc)

ggplot(gene.i.15.both, aes(xmin = start, xmax = end, y = contig_id, fill = focal, label=product_name)) +
  geom_gene_arrow(aes(forward = strand == "+")) +
  geom_gene_label(align = "left")+
  facet_wrap(~ contig.q.vs.s, scales = "free", ncol = 2) +
  scale_fill_manual(values = color6[c(2,4)]) 

#Arrange lgt for plotting

