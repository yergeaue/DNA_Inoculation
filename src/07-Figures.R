
#Figure 1.
box.plant <- readRDS(file=here("data", "intermediate", "box.plant.RDS"))
ggsave(file= here("output", "figures", "Fig1.pdf"), box.plant, width = 7, units = "in")
ggsave(file= here("output", "figures", "Fig1.tiff"), box.plant, width = 7, units = "in", 
       dpi=600, compression = "lzw")
