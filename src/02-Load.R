#Palette
color6 <- c("#016180", "#4a8cb0", "#cde3f0","#EED9C4","#8c5c47", "#4c2c17")
palette(brewer.pal(n=6, name="BrBG"))

#Mapping file
map <- read.table(file=here("data", "raw", "mapping_file.tsv"), header=T, 
                    sep="\t", comment.char = "", row.names = 1)
row.names(map)[40] <- "T1-I1-5"#Bug with row name
rownames(map) <- gsub("-",".",rownames(map))
map <- map[order(row.names(map)),]
saveRDS(map, file = here("data","intermediate", "map.RDS"))

#Plant
plant <- read.table(file=here("data","raw", "plant.txt", sep=""), 
                    header=T, sep="\t",  row.names=1)
saveRDS(plant, file = here("data","intermediate", "plant.RDS"))

#Annotations
#All samples assembly
annot.all <- read.table(file=here("data","raw", "annotations.tsv"), header=T, 
                 sep="\t", row.names=2, comment.char = "", quote="")
annot.all <- annot.all[-4385263,]#remove gene_id_gene_id_000 (not in gene file)
annot.all<- annot.all[order(row.names(annot.all)),]#sort
annot.all$tax_genus <- str_match(annot.all$tax_genus, "[A-Z][a-z]*") #Fix the genus column
annot.all$tax_genus <- gsub("\\bN\\b", "NULL", annot.all$tax_genus) #Fix the genus column
tax.all <- annot.all[,c(1,26:32)]
COG.all <- annot.all[,c(1,16:19)]
saveRDS(annot.all, file = here("data","intermediate", "annot.all.RDS"))
saveRDS(tax.all, file = here("data","intermediate", "tax.all.RDS"))
saveRDS(COG.all, file = here("data","intermediate", "COG.all.RDS"))
#Innoculum assembly
annot.inoc <- read.table(file=here("data","raw", "annotations_inoc.tsv"), header=T, 
                         sep="\t", comment.char = "", quote="") #309078 obs of 32 variables
row.names(annot.inoc) <- annot.inoc$gene_id
annot.inoc <- annot.inoc[-123371,]#remove gene_id_gene_id_000 (not in gene file)
annot.inoc<- annot.inoc |> arrange (gene_id)#sort
annot.inoc$tax_genus <- str_match(annot.inoc$tax_genus, "[A-Z][a-z]*") #Fix the genus column
annot.inoc$tax_genus <- gsub("\\bN\\b", "NULL", annot.inoc$tax_genus) #Fix the genus column
annot.inoc <- annot.inoc[,-33]
saveRDS(annot.inoc, file = here("data","intermediate", "annot.inoc.RDS"))
#Non-inoculated assembly
annot.noninoc <- read.table(file=here("data","raw", "annotations_noninoc.tsv"), header=T, 
                         sep="\t", row.names=2, comment.char = "", quote="") #11095067 obs of 32 variables
annot.noninoc <- annot.noninoc[-123371,]#remove gene_id_gene_id_000 (not in gene file)
annot.noninoc<- annot.noninoc[order(row.names(annot.noninoc)),]#sort
annot.noninoc$tax_genus <- str_match(annot.noninoc$tax_genus, "[A-Z][a-z]*") #Fix the genus column
annot.noninoc$tax_genus <- gsub("\\bN\\b", "NULL", annot.noninoc$tax_genus) #Fix the genus column
saveRDS(annot.noninoc, file = here("data","intermediate", "annot.noninoc.RDS"))


#Genes
#All samples assembly
genes.all <- read.table(file=here("data", "raw", "merged_gene_abundance.tsv"), 
                        header=T, sep="\t", comment.char = "", row.names = 1)
genes.all <- genes.all[order(row.names(genes.all)),]#sort
colnames(genes.all)[4] <- "T1.I1.5"#Bug with row name
genes.all.rel <- data.frame(t(apply(genes.all, 1, "/", colSums(genes.all))))#Normalization
colSums(genes.all.rel)#Sanity check - ok all 1s
saveRDS(genes.all, file = here("data","intermediate", "genes.all.RDS"))
saveRDS(genes.all.rel, file = here("data","intermediate", "genes.all.rel.RDS"))
#Inoculum only assembly -- 2 files, one with the samples and another with the inoculum
genes.inoc <- read.table(file=here("data", "raw", "merged_gene_abundance_inoc.tsv"), 
                        header=T, sep="\t", comment.char = "", row.names = 1) #309078 obs. of 36 variables
genes.inoc.inoc <- read.table(file=here("data", "raw", "merged_gene_abundance_inoc_inoc.tsv"), 
                              header=T, sep="\t", comment.char = "", row.names = 1) #309078 obs. of 2 variables
genes.inoc <- genes.inoc[order(row.names(genes.inoc)),]#sort
genes.inoc.inoc <- genes.inoc.inoc[order(row.names(genes.inoc.inoc)),]#sort
sum(row.names(genes.inoc)==row.names(genes.inoc.inoc)) #309078
colnames(genes.inoc)[3] <- "T1.I1.5"#Bug with row name
genes.inoc <- cbind(genes.inoc.inoc, genes.inoc) #merge
genes.inoc.rel <- data.frame(t(apply(genes.inoc, 1, "/", colSums(genes.inoc))))#Normalization
colSums(genes.inoc.rel)#Sanity check - ok all 1s
saveRDS(genes.inoc, file = here("data","intermediate", "genes.inoc.RDS"))
saveRDS(genes.inoc.rel, file = here("data","intermediate", "genes.inoc.rel.RDS"))

#Genes position on contigs
#Inoculum only assembly
gff.inoc <- read.table(here("data","raw","Contigs_renamed_inoc.gff"), sep="\t") #309078 obs of 9 variables
gff.inoc[c('gene_id', 'Other')] <- str_split_fixed(gff.inoc$V9, ';', 2) #Get gene_id out of last column
gff.inoc <- gff.inoc[,c(1,4,5,7,10)]
colnames(gff.inoc) <- c("contig_id", "start", "end", "strand", "gene_id")
#Non inoculums assembly (also no soil)
gff.noninoc <- read.table(here("data","raw","Contigs_renamed_noninoc.gff"), sep="\t") #11095067 obs of 9 variables
gff.noninoc[c('gene_id', 'Other')] <- str_split_fixed(gff.noninoc$V9, ';', 2) #Get gene_id out of last column
gff.noninoc <- gff.noninoc[,c(1,4,5,7,10)]
colnames(gff.noninoc) <- c("contig_id", "start", "end", "strand", "gene_id")
#Save intermediate
saveRDS(gff.inoc, file = here("data","intermediate", "gff.inoc.RDS"))
saveRDS(gff.noninoc, file = here("data","intermediate", "gff.noninoc.RDS"))

#Blast output
blast_output <- read.table(here("data","raw","blast_output.txt"), sep = "\t")
colnames(blast_output) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                            "qstart", "qend", "sstart", "send", "evalue", "bitscore")
saveRDS(blast_output, file = here("data","intermediate","blast_output.RDS"))
