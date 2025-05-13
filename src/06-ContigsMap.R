##Visualise genes in the contigs
##Compare genes in Inoculum assembly vs. non-inoculum assembly (also contains controls - but selected genes should be absent from control)
##Blast the inoculum genes onto the non-inoculum genes (on the compute cluster: 
##Running from /project/6004719/projects/dna_recruitment_2nd_attempt/non_inoculums/gene_prediction
##blastn -query ../../inoculums_only/gene_prediction/Contigs_renamed.fna  -db Contigs_renamed.fna -outfmt 6 -num_alignments 1 >blast_output.txt)
##blast_output.txt copied in the raw data folder
##Find, visualize and compare the contigs

blast_output <- readRDS(here("data", "intermediate", "blast_output.RDS"))
gff.inoc <- readRDS(here("data", "intermediate", "gff.inoc.RDS")) #Contig information for the inoculum only assembly
gff.noninoc <- readRDS(here("data", "intermediate", "gff.noninoc.RDS")) #Contig information for the non-inoculum only assembly
annot.inoc <- readRDS(file = here("data","intermediate", "annot.inoc.RDS"))
annot.noninoc <- readRDS(file = here("data","intermediate", "annot.noninoc.RDS"))
i.15.table.0 <- readRDS(here("data","intermediate","i.15.table.0.RDS"))
i.50.table.0 <- readRDS(here("data","intermediate","i.50.table.0.RDS"))
ni.15.table.0 <- readRDS(here("data","intermediate","ni.15.table.0.RDS"))
ni.50.table.0 <- readRDS(here("data","intermediate","ni.50.table.0.RDS"))

#I inoculum, 15% SWHC
blast.i.15 <- blast_output |> 
  filter(qseqid%in%i.15.table.0$`Gene ID`) |> #Keep only genes potentially transferred
  left_join(gff.inoc, by = c("qseqid"="gene_id"), keep = TRUE) |> #add contig info for query genes
  left_join(gff.noninoc, by = c("sseqid"="gene_id"), keep = TRUE, suffix = c("_q", "_s")) |> #add contig info for subject genes
  mutate(contig.q.vs.s = interaction(contig_id_q,contig_id_s, sep = " vs ")) #Add column linking the two contigs

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

#I inoculum, 50% SWHC
blast.i.50 <- blast_output |> 
  filter(qseqid%in%i.50.table.0$`Gene ID`) |> #Keep only genes potentially transferred
  left_join(gff.inoc, by = c("qseqid"="gene_id"), keep = TRUE) |> #add contig info for query genes
  left_join(gff.noninoc, by = c("sseqid"="gene_id"), keep = TRUE, suffix = c("_q", "_s")) |> #add contig info for subject genes
  mutate(contig.q.vs.s = interaction(contig_id_q,contig_id_s, sep = " vs ")) #Add column linking the two contigs

gene.i.50.inoc <- gff.inoc |>
  filter(contig_id %in% blast.i.50$contig_id_q) |> #Keep only contigs from blast file
  left_join(select(annot.inoc, gene_id, product_name, kegg_entry, kegg_definition, tax_phylum)) |> #Add annotations (selected var)
  left_join(select(blast.i.50, contig_id_q, contig.q.vs.s), by = c("contig_id"="contig_id_q")) |> #Add contig linkage info
  mutate(focal = gene_id%in%blast.i.50$qseqid) |>
  mutate(assembly = "inoculum")

gene.i.50.noninoc <- gff.noninoc |>
  filter(contig_id %in% blast.i.50$contig_id_s) |> #Keep only contigs from blast file
  left_join(select(annot.noninoc, gene_id, product_name, kegg_entry, kegg_definition, tax_phylum)) |> #Add annotations (selected var)
  left_join(select(blast.i.50, contig_id_s, contig.q.vs.s), by = c("contig_id"="contig_id_s")) |> #Add contig linkage info
  mutate(focal = gene_id%in%blast.i.50$sseqid) |>
  mutate(assembly = "non-inoculum")

gene.i.50.both <- rbind(gene.i.50.inoc, gene.i.50.noninoc)

ggplot(gene.i.50.both, aes(xmin = start, xmax = end, y = contig_id, fill = focal, label=product_name)) +
  geom_gene_arrow(aes(forward = strand == "+")) +
  geom_gene_label(align = "left")+
  facet_wrap(~ contig.q.vs.s, scales = "free", ncol = 2) +
  scale_fill_manual(values = color6[c(2,4)]) 

#NI inoculum, 15% SWHC
blast.ni.15 <- blast_output |> 
  filter(qseqid%in%ni.15.table.0$`Gene ID`) |> #Keep only genes potentially transferred
  left_join(gff.inoc, by = c("qseqid"="gene_id"), keep = TRUE) |> #add contig info for query genes
  left_join(gff.noninoc, by = c("sseqid"="gene_id"), keep = TRUE, suffix = c("_q", "_s")) |> #add contig info for subject genes
  mutate(contig.q.vs.s = interaction(contig_id_q,contig_id_s, sep = " vs ")) #Add column linking the two contigs

gene.ni.15.inoc <- gff.inoc |>
  filter(contig_id %in% blast.ni.15$contig_id_q) |> #Keep only contigs from blast file
  left_join(select(annot.inoc, gene_id, product_name, kegg_entry, kegg_definition, tax_phylum)) |> #Add annotations (selected var)
  left_join(select(blast.ni.15, contig_id_q, contig.q.vs.s), by = c("contig_id"="contig_id_q")) |> #Add contig linkage info
  mutate(focal = gene_id%in%blast.ni.15$qseqid) |>
  mutate(assembly = "inoculum")

gene.ni.15.noninoc <- gff.noninoc |>
  filter(contig_id %in% blast.ni.15$contig_id_s) |> #Keep only contigs from blast file
  left_join(select(annot.noninoc, gene_id, product_name, kegg_entry, kegg_definition, tax_phylum)) |> #Add annotations (selected var)
  left_join(select(blast.ni.15, contig_id_s, contig.q.vs.s), by = c("contig_id"="contig_id_s")) |> #Add contig linkage info
  mutate(focal = gene_id%in%blast.ni.15$sseqid) |>
  mutate(assembly = "non-inoculum")

gene.ni.15.both <- rbind(gene.ni.15.inoc, gene.ni.15.noninoc)

ggplot(gene.ni.15.both, aes(xmin = start, xmax = end, y = contig_id, fill = focal, label=product_name)) +
  geom_gene_arrow(aes(forward = strand == "+")) +
  geom_gene_label(align = "left")+
  facet_wrap(~ contig.q.vs.s, scales = "free", ncol = 2) +
  scale_fill_manual(values = color6[c(2,4)]) 

#NI inoculum, 50% SWHC
blast.ni.50 <- blast_output |> 
  filter(qseqid%in%ni.50.table.0$`Gene ID`) |> #Keep only genes potentially transferred
  left_join(gff.inoc, by = c("qseqid"="gene_id"), keep = TRUE) |> #add contig info for query genes
  left_join(gff.noninoc, by = c("sseqid"="gene_id"), keep = TRUE, suffix = c("_q", "_s")) |> #add contig info for subject genes
  mutate(contig.q.vs.s = interaction(contig_id_q,contig_id_s, sep = " vs ")) #Add column linking the two contigs

gene.ni.50.inoc <- gff.inoc |>
  filter(contig_id %in% blast.ni.50$contig_id_q) |> #Keep only contigs from blast file
  left_join(select(annot.inoc, gene_id, product_name, kegg_entry, kegg_definition, tax_phylum)) |> #Add annotations (selected var)
  left_join(select(blast.ni.50, contig_id_q, contig.q.vs.s), by = c("contig_id"="contig_id_q")) |> #Add contig linkage info
  mutate(focal = gene_id%in%blast.ni.50$qseqid) |>
  mutate(assembly = "inoculum")

gene.ni.50.noninoc <- gff.noninoc |>
  filter(contig_id %in% blast.ni.50$contig_id_s) |> #Keep only contigs from blast file
  left_join(select(annot.noninoc, gene_id, product_name, kegg_entry, kegg_definition, tax_phylum)) |> #Add annotations (selected var)
  left_join(select(blast.ni.50, contig_id_s, contig.q.vs.s), by = c("contig_id"="contig_id_s")) |> #Add contig linkage info
  mutate(focal = gene_id%in%blast.ni.50$sseqid) |>
  mutate(assembly = "non-inoculum")

gene.ni.50.both <- rbind(gene.ni.50.inoc, gene.ni.50.noninoc)

ggplot(gene.ni.50.both, aes(xmin = start, xmax = end, y = contig_id, fill = focal,  label=product_name)) +
  geom_gene_arrow(aes(forward = strand == "+")) +
  geom_gene_label(align = "left")+
  facet_wrap(~ contig.q.vs.s, scales = "free", ncol = 2) +
  scale_fill_manual(values = color6[c(2,4)]) 
