#Palette
color6 <- c("#016180", "#4a8cb0", "#cde3f0","#EED9C4","#8c5c47", "#4c2c17")
palette(brewer.pal(n=6, name="BrBG"))
#Plant
plant <- read.table(file=here("data","raw", "plant.txt", sep=""), 
                    header=T, sep="\t",  row.names=1)
saveRDS(plant, file = here("data","intermediate", "plant.RDS"))


