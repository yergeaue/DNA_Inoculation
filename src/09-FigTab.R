
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

#Table 2
table2 <- gt(i.15.table) |>
  tab_header(title = "Table 2", subtitle="Potential naturally transferred genes from the irrigated inoculum under LW")
gtsave(table2, filename = here("output","tables","table2.docx"))

#Table 3
table3 <- gt(ni.15.table) |>
  tab_header(title = "Table 3", subtitle="Potential naturally transferred genes from the ambient inoculum under LW")
gtsave(table3, filename = here("output","tables","table3.docx"))

#Table 4
table4 <- gt(i.50.table) |>
  tab_header(title = "Table 4", subtitle="Potential naturally transferred genes from the irrigated inoculum under HW")
gtsave(table4, filename = here("output","tables","table4.docx"))

#Table 5
table5 <- gt(ni.50.table) |>
  tab_header(title = "Table 5", subtitle="Potential naturally transferred genes from the ambient inoculum under HW")
gtsave(table5, filename = here("output","tables","table5.docx"))

#Table S1
tableS1 <- gt(lgt.genes.pub) |>
  tab_header(title = "Table S1", subtitle="Potentially LGT genes identified by WAAFLE")
gtsave(tableS5, filename = here("output","tables","tableS1.docx"))

#Table S2
tableS2 <- gt(tax.trans.0.def) |>
  tab_header(title = "Table S2", subtitle="Potentially LGT genes overrepresented in inoculated samples")
gtsave(tableS2, filename = here("output","tables","tableS2.docx"))
