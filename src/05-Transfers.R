###Find genes potentially transferred
##Find genes
#Load data: all samples mapped on assembly of inoculum
genes.inoc.rel <- readRDS(here("data","intermediate", "genes.inoc.rel.RDS"))
map <- readRDS(here("data", "intermediate", "map.RDS"))
map.nosoil <- map[map$SoilWaterContent != "Soil",]
genes.inoc.rel <- genes.inoc.rel[,order(colnames(genes.inoc.rel))]
sum(rownames(map.nosoil)==colnames(genes.inoc.rel)) #38

#Select genes that are 
#1) present in inoculum (given, since assembly if from inoculum only)
#2)Absent in controls
genes.inoc.rel.absent <- genes.inoc.rel[rowSums(genes.inoc.rel[,map.nosoil$Inoculum=="ctrl"])==0,  ] #subset to absent in non inoculated
#3)Present in rhizosphere
genes.inoc.rel.absent <- genes.inoc.rel.absent[rowSums(genes.inoc.rel.absent[,c(3:14,21:32)])!=0,] #subset to present in the rhizosphere
#Keep only genes that are defined at the kegg_entry
annot.inoc.def <- annot.inoc[annot.inoc$kegg_entry != "NULL", ]
genes.inoc.rel.absent <- genes.inoc.rel.absent[row.names(genes.inoc.rel.absent) %in% row.names(annot.inoc.def),]

#4)Abundant at a safety margin (x) above relics
dil <- 850/(50*100/0.5*1500/1000) #Dilution factor - added 850 ug in 15,000ug
margin <- 10

#Irrigated inoculum, 15% SWHC
thresh.i.15 <- margin*dil*genes.inoc.rel.absent[,1] #threshold for each gene
count.i.15 <- rowSums(genes.inoc.rel.absent[,(map.nosoil$Inoculum =="i" & map.nosoil$SoilWaterContent == "15%")]>thresh.i.15)
sum(count.i.15>1)
i.15 <- data.frame(thresh.i.15, count.i.15, 
                   genes.inoc.rel.absent[,(map.nosoil$Inoculum =="i" & 
                                             (map.nosoil$SoilWaterContent == "15%" |
                                                map.nosoil$SoilWaterContent == "Inoculum"))])

i.15.annot <- annot.inoc.def[row.names(annot.inoc.def)%in%row.names(i.15),]
sum(row.names(i.15)==row.names(i.15.annot))#1060
i.15.final <- cbind(i.15, i.15.annot[,c(1,3,4,27,29,31)])#Add annotations
i.15.final <- i.15.final[i.15.final$count.i.15>=2,] #Keep only genes transferred 2 times and above
i.15.pub <- i.15.final[,c(2,11:15)] #Keep only relevant columns for publication
i.15.pub$tax_genus <- str_match(i.15.pub$tax_genus, "[A-Z][a-z]*") #Fix the genus column
i.15.pub$tax_genus <- gsub("\\bN\\b", "NULL", i.15.pub$tax_genus) #Fix the genus column
i.15.pub$gene_id <- row.names(i.15.final) #Put back the row names
i.15.pub$gene_id <- NULL
i.15.pub <- i.15.pub[order(i.15.pub$count.i.15, decreasing = TRUE),] #Sort by number of samples above threshold
colnames(i.15.pub) <- c("Count", "KEGG entry", "KEGG definition", "Phylum", "Order", "Genus")
write.table(i.15.pub, file = here("output", "tables", "Table2.txt"), sep = "\t")

#Ambient inoculum, 15% SWHC
thresh.ni.15 <- margin*dil*genes.inoc.rel.absent[,2] #threshold for each gene
count.ni.15 <- rowSums(genes.inoc.rel.absent[,(map.nosoil$Inoculum =="ni" & map.nosoil$SoilWaterContent == "15%")]>thresh.ni.15)
sum(count.ni.15>1)
ni.15 <- data.frame(thresh.ni.15, count.ni.15, 
                   genes.inoc.rel.absent[,(map.nosoil$Inoculum =="ni" & 
                                             (map.nosoil$SoilWaterContent == "15%" |
                                                map.nosoil$SoilWaterContent == "Inoculum"))])

ni.15.annot <- annot.inoc.def[row.names(annot.inoc.def)%in%row.names(ni.15),]
sum(row.names(ni.15)==row.names(ni.15.annot))#1060
ni.15.final <- cbind(ni.15, ni.15.annot[,c(1,3,4,27,29,31)])#Add annotations
ni.15.final <- ni.15.final[ni.15.final$count.ni.15>=2,] #Keep only genes transferred 2 times and above
ni.15.pub <- ni.15.final[,c(2,11:15)] #Keep only relevant columns for publication
ni.15.pub$tax_genus <- str_match(ni.15.pub$tax_genus, "[A-Z][a-z]*") #Fix the genus column
ni.15.pub$tax_genus <- gsub("\\bN\\b", "NULL", ni.15.pub$tax_genus) #Fix the genus column
ni.15.pub$gene_id <- row.names(ni.15.final) #Put back the row names
ni.15.pub$gene_id <- NULL
ni.15.pub <- ni.15.pub[order(ni.15.pub$count.ni.15, decreasing = TRUE),] #Sort by number of samples above threshold
colnames(ni.15.pub) <- c("Count", "KEGG entry", "KEGG definition", "Phylum", "Order", "Genus")
write.table(ni.15.pub, file = here("output", "tables", "Table3.txt"), sep = "\t")

#Irrigated inoculum, 50% SWHC
thresh.i.50 <- margin*dil*genes.inoc.rel.absent[,1] #threshold for each gene
count.i.50 <- rowSums(genes.inoc.rel.absent[,(map.nosoil$Inoculum =="i" & map.nosoil$SoilWaterContent == "50%")]>thresh.i.50)
sum(count.i.50>1)
i.50 <- data.frame(thresh.i.50, count.i.50, 
                   genes.inoc.rel.absent[,(map.nosoil$Inoculum =="i" & 
                                             (map.nosoil$SoilWaterContent == "50%" |
                                                map.nosoil$SoilWaterContent == "Inoculum"))])

i.50.annot <- annot.inoc.def[row.names(annot.inoc.def)%in%row.names(i.50),]
sum(row.names(i.50)==row.names(i.50.annot))#1060
i.50.final <- cbind(i.50, i.50.annot[,c(1,3,4,27,29,31)])#Add annotations
i.50.final <- i.50.final[i.50.final$count.i.50>=2,] #Keep only genes transferred 2 times and above
i.50.pub <- i.50.final[,c(2,11:15)] #Keep only relevant columns for publication
i.50.pub$tax_genus <- str_match(i.50.pub$tax_genus, "[A-Z][a-z]*") #Fix the genus column
i.50.pub$tax_genus <- gsub("\\bN\\b", "NULL", i.50.pub$tax_genus) #Fix the genus column
i.50.pub$gene_id <- row.names(i.50.final) #Put back the row names
i.50.pub$gene_id <- NULL
i.50.pub <- i.50.pub[order(i.50.pub$count.i.50, decreasing = TRUE),] #Sort by number of samples above threshold
colnames(i.50.pub) <- c("Count", "KEGG entry", "KEGG definition", "Phylum", "Order", "Genus")
write.table(i.50.pub, file = here("output", "tables", "Table4.txt"), sep = "\t")

#Ambient inoculum, 50% SWHC
thresh.ni.50 <- margin*dil*genes.inoc.rel.absent[,2] #threshold for each gene
count.ni.50 <- rowSums(genes.inoc.rel.absent[,(map.nosoil$Inoculum =="ni" & map.nosoil$SoilWaterContent == "50%")]>thresh.ni.50)
sum(count.ni.50>1)
ni.50 <- data.frame(thresh.ni.50, count.ni.50, 
                   genes.inoc.rel.absent[,(map.nosoil$Inoculum =="ni" & 
                                             (map.nosoil$SoilWaterContent == "50%" |
                                                map.nosoil$SoilWaterContent == "Inoculum"))])

ni.50.annot <- annot.inoc.def[row.names(annot.inoc.def)%in%row.names(ni.50),]
sum(row.names(ni.50)==row.names(ni.50.annot))#1060
ni.50.final <- cbind(ni.50, ni.50.annot[,c(1,3,4,27,29,31)])#Add annotations
ni.50.final <- ni.50.final[ni.50.final$count.ni.50>=2,] #Keep only genes transferred 2 times and above
ni.50.pub <- ni.50.final[,c(2,11:15)] #Keep only relevant columns for publication
ni.50.pub$tax_genus <- str_match(ni.50.pub$tax_genus, "[A-Z][a-z]*") #Fix the genus column
ni.50.pub$tax_genus <- gsub("\\bN\\b", "NULL", ni.50.pub$tax_genus) #Fix the genus column
ni.50.pub$gene_id <- row.names(ni.50.final) #Put back the row names
ni.50.pub$gene_id <- NULL
ni.50.pub <- ni.50.pub[order(ni.50.pub$count.ni.50, decreasing = TRUE),] #Sort by number of samples above threshold
colnames(ni.50.pub) <- c("Count", "KEGG entry", "KEGG definition", "Phylum", "Order", "Genus")
write.table(ni.50.pub, file = here("output", "tables", "Table5.txt"), sep = "\t")

##Visualise genes in the contigs
##Compare genes in Inoculum assembly vs. samples assembly (also contains controls - but selected genes should be absent from control - see above)
##Use Blast to link the genes to one another
##Find, visualize and compare the contigs

gff.inoc[gff.inoc$gene_id == "gene_id_196738",]
gff.inoc[gff.inoc$contig_id == "k91_342706",]


gff.noninoc[gff.noninoc$gene_id == "gene_id_4999252",]
gff.noninoc[gff.noninoc$contig_id == "k91_5127639",]


library(gggenes)
library(ggplot2)
library(data.table)

# Example data
genes_df <- data.frame(
  molecule = c("contig1", "contig1", "contig1", "contig2", "contig2"),
  gene = c("gene1", "gene2", "gene3", "gene4", "gene5"),
  start = c(1, 500, 1200, 1, 800),
  end = c(400, 900, 1800, 600, 1300),
  strand = c(1, -1, 1, 1, -1),
  type = c("typeA", "typeB", "typeA", "typeC", "typeB")
)

# Create the plot
ggplot(genes_df, aes(xmin = start, xmax = end, y = molecule, fill = type)) +
  geom_gene_arrow(aes(forward = strand == 1)) +
  facet_wrap(~ molecule, scales = "free") +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()


# Example with bed file:
setwd("/project/6004719/projects/GROW/MT_2024-10-06/bacteria/")
bed_file <- "./gene_prediction/Contigs_genes.bed"
df <- data.frame(fread(bed_file))
colnames(df) = c("contig_id", "start", "end", "gene_id")
df$strand = 1

df2 <- df[df$contig_id %in% c("k91_11481"),]
head(df2)
ggplot(df2, aes(xmin=start, xmax=end, y=contig_id)) +
  geom_gene_arrow(aes(forward = strand == 1), fill="blue") +
  facet_wrap(~ contig_id, scales = "free") +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()

