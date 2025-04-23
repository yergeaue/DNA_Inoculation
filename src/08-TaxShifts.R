#Find contigs where genes have discordant taxonomy --> indicative of HGT, as in 
# https://academic.oup.com/ismecommun/article/4/1/ycae073/7671050
#
#
#Load data
#Gene_id in same order as Contig_id_N
#Contig_id_N
contig.N.gene <- cbind(contig_id, gene_id)#Merge
colnames(contig.N.gene) <- c("contig_id_N", "gene_id")
contig.N.gene.contig <- contig.N.gene |>  #Add contig names
  left_join(select(annot.noninoc, gene_id, X.contig_id))
#Assigned taxonomy strings for each Contig_id_N
colnames(gene.tax) <- c("contig_id_N", "gene_N", "tax_string", "bit_score")
#Taxonomy names corresponding to number in taxonomy strings


#Find strings that are dissimilar for genes on the same contig
gene.tax.contig <- gene.tax |>
  left_join(contig.N.gene.contig)

gene.tax.contig.distinct <- gene.tax.contig |>
  distinct(interaction(X.contig_id, tax_string))


gene.tax.contig[gene.tax.contig$X.contig_id%in%gene.i.15.noninoc$contig_id,]
