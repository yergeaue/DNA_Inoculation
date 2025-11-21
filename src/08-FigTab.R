
#Figure 1.
box.plant <- readRDS(file=here("data", "intermediate", "box.plant.RDS"))
ggsave(file= here("output", "figures", "Fig1.pdf"), box.plant, width = 7, units = "in")
ggsave(file= here("output", "figures", "Fig1.tiff"), box.plant, width = 7, units = "in", 
       dpi=600, compression = "lzw")

#Figure 2.
fig2<-ggarrange(ordi.noninoc, stack.all, nrow = 2, labels = c("a","b"))
fig2
ggsave(file= here("output", "figures", "Fig2.pdf"), fig2, width = 7, height = 11, units = "in")
ggsave(file= here("output", "figures", "Fig2.tiff"), fig2, width = 7, height = 11, units = "in", 
       dpi=600, compression = "lzw")

#Figure 3.
fig3<-ggarrange(stack.trans.tax, stack.trans.fun, nrow = 2, labels = c("a","b"))
fig3
ggsave(file= here("output", "figures", "Fig3.pdf"), fig3, width = 7, height = 11, units = "in")
ggsave(file= here("output", "figures", "Fig3.tiff"), fig3, width = 7, height = 11, units = "in", 
       dpi=600, compression = "lzw")

#Table S1
tableS1 <- gt(lgt.genes.pub) |>
  tab_header(title = "Table S1", subtitle="Potentially LGT genes identified by WAAFLE")
gtsave(tableS1, filename = here("output","tables","tableS1.docx"))

#Table S2
tableS2 <- gt(summary.stat.lgt.SWC) |>
  tab_header(title = "Table S2", subtitle="Potentially LGT genes identified by WAAFLE that are significantly affected by soil water content")
gtsave(tableS2, filename = here("output","tables","tableS2.docx"))

#Table S3
tableS3 <- gt(trans.pub) |>
  tab_header(title = "Table S3", subtitle="Potentially LGT genes overrepresented in inoculated samples")
gtsave(tableS3, filename = here("output","tables","tableS3.docx"))
