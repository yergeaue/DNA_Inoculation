#Find contigs where genes have discordant taxonomy --> indicative of HGT, as in 
# https://academic.oup.com/ismecommun/article/4/1/ycae073/7671050

#Load data
gene.id <- readRDS(file = here("data","intermediate", "gene.id.RDS")) #Gene_id in same order as Contig_id_N
contig.id <- readRDS(file = here("data","intermediate", "contig.id.RDS")) #Contig_id_N
gene.tax <- readRDS(file = here("data","intermediate", "gene.tax.RDS")) #Taxonomy assigned to each genes. NCBI encoded
tax.names <- readRDS(file = here("data","intermediate", "tax.names.RDS")) #NCBI codes for taxonomy
annot.noninoc <- readRDS(file = here("data","intermediate", "annot.noninoc.RDS")) #Annotations for non-inoculum assembly

#Arrange data
contig.N.gene <- cbind(contig.id, gene.id)#Merge
colnames(contig.N.gene) <- c("contig_id_N", "gene_id") 
contig.N.gene.contig <- contig.N.gene |>  
  left_join(select(annot.noninoc, gene_id, contig_id)) #Add contig names, without the N
colnames(gene.tax) <- c("contig_id_N", "gene_N", "tax_string", "bit_score")
tax.names.less <- tax.names[tax.names$V7=="scientific name", c(1,3,5,7)] #Keep only relevant columns and lines

#Create object containing contigs that have more than one phylum assignation
gene.tax.trans <- gene.tax |>
  left_join(contig.N.gene.contig)|> #Add taxonomy strings
  separate(tax_string, c("Root", "Subroot", "Domain", "Superphylum", "Phylum"), 
                                             fill = "right") |> #Split taxonomy string 
  left_join(tax.names.less, join_by(Phylum == V1) ) |> #Taxonomy names corresponding to number in taxonomy strings
  filter(!is.na(V3)) |> #Remove NAs
  group_by(contig_id)  |> 
  reframe(NbPhylum = n_distinct(V3), gene_id = gene_id, Phylum = V3) |> #Count nb of phylum per contig
  filter(NbPhylum>1)

saveRDS(gene.tax.trans, file = here("data","intermediate", "gene.tax.trans.RDS"))
n_distinct(gene.tax.trans$contig_id)#How many suspicious contigs ? 257,853 - 4.66%
n_distinct(annot.noninoc$contig_id)#How many original? 5,532,321

#Compare contigs to the ones found in 06-ContigsMap.R
trans.intersect.i.15 <- intersect(gene.tax.trans$contig_id, gene.i.15.noninoc$contig_id) # two contigs: k91_3243898 and k91_3901876
trans.intersect.i.50 <- intersect(gene.tax.trans$contig_id, gene.i.50.noninoc$contig_id) # two contig : k91_3901876 and k91_11805324
trans.intersect.ni.15 <- intersect(gene.tax.trans$contig_id, gene.ni.15.noninoc$contig_id) # two contigs: k91_3243898 and k91_3901876
trans.intersect.ni.50 <- intersect(gene.tax.trans$contig_id, gene.ni.50.noninoc$contig_id) # three contigs :k91_3243898, k91_3901876 and k91_11805324

#Focus on these three contigs: k91_3243898, k91_3901876 and k91_11805324
focus.contigs <- gene.ni.50.noninoc |>
  filter(contig_id %in% c("k91_3243898", "k91_3901876", "k91_11805324")) |>
  left_join(gene.tax.trans, join_by(gene_id))

ggplot(focus.contigs, aes(xmin = start, xmax = end, y = contig_id.x, fill = Phylum, color = focal, label=product_name)) +
  geom_gene_arrow(aes(forward = strand == "+")) +
  geom_gene_label(align = "left")+
  facet_wrap(~ contig_id.x, scales = "free", ncol = 1) +
  scale_fill_manual(values = color6) +
  scale_color_manual(values = c("white", "black"))
