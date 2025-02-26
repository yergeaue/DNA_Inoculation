#Plant
plant <- read.table(file=here("data","raw", "plant.txt", sep=""), 
                    header=T, sep="\t",  row.names=1)
saveRDS(plant, file = here("data","intermediate", "plant.RDS"))
