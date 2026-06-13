# TITLE:        Belowground Plot Experiment - Initial soil microbial data set wrangling
# AUTHOR:       Matthew Nieland
# 
# DATA INPUT:   Features_2017_2021BGPEDataSamplesOnly.qza, 
#               MCCRecovery_AllBGPESamplesMappingFile.txt, 
#               Taxonomy_2017_2021BGPEDataSamplesOnly.qza
#
# DATA OUTPUT:  phyloseq.rds  
#
# PROJECT:      Belowground Plot Experiment (KNZ LTER) fertilization legacies
#
# DATE:         June 2026

# This script is used to create a phyloseq object (as .rds) that will be used for 
# subsequent microbial community analysis. 

# The following immediate script downloads necessary packages for this initial 
# setup.


# Add "BiocManager" and "remotes" packages. "BiocManager" package helps install 
# software from the Bioconductor project.

# "remotes" package will be used to install the "qiime2R" package, which will 
# process uploading .qza data a lot easier.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager") # I have v1.30.26

if (!require("remotes", quietly = TRUE))
  install.packages("remotes") # I have v2.5.0

# Download packages "phyloseq" from Bioconductor and "qiime2R" created by 
# Jordan Bisanz (https://github.com/jbisanz/qiime2R).

BiocManager::install("phyloseq")
remotes::install_github("jbisanz/qiime2R")





# Clear the environment in R workspace.
rm(list=ls())

# Load packages
library(phyloseq) # v0.99.6
library(qiime2R) # v1.48.0

# Set working directory 
setwd("~/Desktop/Research/BGPE/MCCRecovery") # directory path from my laptop


# The subsequent script consolidates feature table, metadata, and taxonomy of 
# ASVs into a phyloseq object.

# Import amplicon sequence variant (ASV) feature table
asv_table = qza_to_phyloseq(features="Data/Features_2017_2021BGPEDataSamplesOnly.qza")
ASV	=	otu_table(as.matrix(asv_table),	taxa_are_rows	=	TRUE)

# Import metadata
metadata	= read.table("Data/MCCRecovery_AllBGPESamplesMappingFile.txt",	
                      header = T, row.names=1) 
META = sample_data(metadata)

# Import taxonomy
TAX1 = read_qza("Data/Taxonomy_2017_2021BGPEDataSamplesOnly.qza")
# In its current state, taxonomy for each ASV is one long, continuous chain.
# It is easier to split taxonomy into groups (phylum, class, order, family,
# genus, and "species")
# TAX2 object has taxonomy parsed into separate columns
TAX2 = parse_taxonomy(TAX1$data)

head(TAX2) # check
TAX = tax_table(as.matrix(TAX2))


# Combine groups to make a phyloseq object
physeq = phyloseq(ASV, META, TAX)


# Save phyloseq as an individual object that can be reloaded later.
saveRDS(physeq, file = "phyloseq.rds")

