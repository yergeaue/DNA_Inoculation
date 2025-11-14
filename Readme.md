---
title: "Readme"
output: html_document
---

This R project was used to analyse the data of Giard-Laliberté et al. "Natural genetic transformation of the wheat rhizosphere microbial communities through DNA inoculations"

Data is available on Zenodo: https://doi.org/10.5281/zenodo.3755973. Data should be copied in the /data/raw folder.

Scripts are to be run in the following order:

**01-Packages.R** - Installing and loading packages necessary for the Project.

**02-Load.R** - Load all the raw data, pre-process and save as intermediary files.

**03-Plant.R** - Analysis of the plant data.

**04-Community.R** - Analysis of the microbial community: PCOA, permanova, stack barcharts

**05-Waafle-LGT.R** - Analysis of the WAAFLE output.

**06-Transfers.R** - Identification of potentially transferred genes

**07-TransTax.R** - Compare the relative abundance of putatively transferred genes vs. rhizosphere and inoculum

**08-FigTab.R** - Create multipannel figures and tables for publication.
