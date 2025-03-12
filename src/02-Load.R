#Palette
color6 <- c("#016180", "#4a8cb0", "#cde3f0","#EED9C4","#8c5c47", "#4c2c17")
palette(brewer.pal(n=6, name="BrBG"))

#Plant
plant <- read.table(file=here("data","raw", "plant.txt", sep=""), 
                    header=T, sep="\t",  row.names=1)
saveRDS(plant, file = here("data","intermediate", "plant.RDS"))

#Annotations
annot.all <- read.table(file=here("data","raw", "annotations.tsv"), header=T, 
                 sep="\t", row.names=2, comment.char = "", quote="")
annot.all <- annot.all[-4385263,]#remove gene_id_gene_id_000 (not in gene file)
annot.all<- annot.all[order(row.names(annot.all)),]#sort
tax.all <- annot.all[,c(1,26:32)]
COG.all <- annot.all[,c(1,16:19)]
saveRDS(annot.all, file = here("data","intermediate", "annot.all.RDS"))
saveRDS(tax.all, file = here("data","intermediate", "tax.all.RDS"))
saveRDS(COG.all, file = here("data","intermediate", "COG.all.RDS"))

#Genes
genes.all <- read.table(file=here("data", "raw", "merged_gene_abundance.tsv"), 
                        header=T, sep="\t", comment.char = "", row.names = 1)
genes.all <- genes.all[order(row.names(genes.all)),]#sort
genes.all.rel <- data.frame(t(apply(genes.all, 1, "/", colSums(genes.all))))#Normalization
colSums(genes.all.rel)#Sanity check - ok all 1s
saveRDS(genes.all, file = here("data","intermediate", "genes.all.RDS"))
saveRDS(genes.all.rel, file = here("data","intermediate", "genes.all.rel.RDS"))
