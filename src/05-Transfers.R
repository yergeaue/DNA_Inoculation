###Find genes potentially transferred
##Find genes
#Load data: all samples mapped on assembly of inoculum
genes.inoc.rel <- readRDS(here("data","intermediate", "genes.inoc.rel.RDS"))
map <- readRDS(here("data", "intermediate", "map.RDS"))
map.nosoil <- map[map$SoilWaterContent != "Soil",]
genes.inoc.rel <- genes.inoc.rel[,order(colnames(genes.inoc.rel))]
sum(rownames(map.nosoil)==colnames(genes.inoc.rel)) #38

#Find genes that are 1) present in inoculum (given, since assembly if from inoculum only)
#                   2) absent in the non-inoculated samples
#                   3) are present in the rhizosphere
#                   4) are abundant at a level that is unlikely due to relic (see paper)
genes.inoc.rel.absent <- genes.inoc.rel[rowSums(genes.inoc.rel[,map.nosoil$Inoculum=="ctrl"])==0,  ] #subset to absent in non inoculated
genes.inoc.rel.absent <- genes.inoc.rel.absent[rowSums(genes.inoc.rel.absent[,c(3:14,21:32)])!=0,] #subset to present in the rhizosphere
dil <- 850/(50*100/0.5*1500/1000) #Dilution factor - added 850 ug in 15,000ug
genes.inoc.rel.thresh.i <- genes.inoc.rel.absent[apply(genes.inoc.rel.absent[,c(3:14,21:32)], 1, max)
                                                 >100*dil*genes.inoc.rel.absent[,1],]
genes.inoc.rel.thresh.i <- genes.inoc.rel.thresh.i[genes.inoc.rel.thresh.i[,1]!=0,]
genes.inoc.rel.thresh.ni <- genes.inoc.rel.absent[apply(genes.inoc.rel.absent[,c(3:14,21:32)], 1, max)
                                                  >100*dil*genes.inoc.rel.absent[,2],]
genes.inoc.rel.thresh.ni <- genes.inoc.rel.thresh.ni[genes.inoc.rel.thresh.ni[,2]!=0,]

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

