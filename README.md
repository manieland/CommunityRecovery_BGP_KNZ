
Open research: The data from this study are publically available in the Environmental Data Initiative repository and the National Center for Biotechnology Information. Data can be found in the following links:

[Plant community dataset](https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-knz.17.19)

[Raw 16S rRNA gene sequence data](https://www.ncbi.nlm.nih.gov/bioproject/)
- PRJNA577961
- PRJNA1478011

[Precipitation dataset](https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-knz.4.22)


## CommunityRecovery_BGP_KNZ

# Introduction
This repository contains R scripts and the 16S rRNA gene community R object that clean, analyze, and plot data from a five year study at the Konza Prairie Biological Station (KNZ). The purpose of this study was to assess how a tallgrass prairie ecosystem responds to the cessation of long-term (30 years) nitrogen enrichment. Specifically, we had two questions: (1) Does cessation allow recovery of plant and soil microbial communities?, and (2) Is this recovery affected by annual prescribed fires?

This study took place at the Belowground Plot Experiment. From 1986 to 2016, nitrogen was applied every year in the late spring (typically May). Starting in 2017, long-term enrichment was terminated at the subplot scale to study ecosystem recovery. The altered treatments at the subplot scale are referred to as the "Recovering" treatment. We continued applying nitrogen fertilizer at a smaller scale within recovering plots to continue assessing soil microbial community responses to enrichment. This treatment is referred to as "still-fertilized" (or simply "Fertilized" in R scripts).

![Diagram of the Belowground Plot Experiment at Konza Prairie Biological Station.](https://github.com/manieland/CommunityRecovery_BGP_KNZ/blob/main/BGPELayout.pdf)

Plant community composition data were initially collected every five years. Beginning in 2016, vegetation composition was collected every year. 16S rRNA gene community composition data were collected in June and August for each growing season during the study (2017-2021). 

# Workflow
All data to replicate analyses and figures are located in this repository. Each script is split for each community (plant and soil microbial communities), and the script to calculate the cumulative amount of precipitation within a single growing season is also included. For transparency, we included the script to transform the 16S rRNA gene data (the .qza file for those familiar with these type of data) into a phyloseq object. All other data can be downloaded directly using the code in the R scripts.

# Usage
Main analyses were conducted using R version 4.4.1 (2024-06-14) (R Core Team 2021) within the RStudio environment 2024.04.02+764 (Posit team 2024).

Attached package versions: DESeq2 1.44.0, SummarizedExperiment 1.34.0, Biobase 2.64.0, MatrixGenerics 1.16.0, matrixStats 1.5.0, GenomicRanges 1.56.2, GenomeInfoDb 1.40.1, IRanges 2.38.1, S4Vectors 0.42.1, BiocGenerics 0.50.0, phyloseq 1.48.0, tidyverse 2.0.0, TH.data 1.1-3, MASS 7.3-65, survival 3.8-3, mvtnorm 1.3-3, ggpubr 0.6.1, cowplot 1.2.0, emmeans 1.11.2, lmerTest 3.1-3, lme4 1.1-37, Matrix 1.7-3, vegan 2.7-1, permute 0.9-8, lubridate 1.9.4, forcats 1.0.0, stringr 1.5.1, dplyr 1.1.4, purrr 1.1.0, readr 2.1.5, tidyr 1.3.1, tibble 3.3.0, ggplot2 3.5.2 

Full session information (including dependency versions) is provided in sessionInfo.txt.

# Contact Information
For inquiries related to data and/or code, please contact Matthew Nieland at mnieland AT umass.edu

# Funding
This work was supported by the National Science Foundation - Division of Enivironmental Biology, grant DEB-2025849.
