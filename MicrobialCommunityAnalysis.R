# TITLE:        Belowground Plot Experiment - Microbial community comp. analysis 
# AUTHOR:       Matthew Nieland
# 
# DATA INPUT:   phyloseq.rds
#
# OUTPUT:       TotalNumberofReads_ASVsSamples.tiff
#               Microbial_NMDS.tiff
#               BPG_MicrobialAlphaDiversity.csv
#               Alpha_MicrobialRecoveryLineGraph.tiff
#               MicrobialFunctionalGroup_LineGraph.tiff
#               MicrobialFunctionalGroup_Boxplot.tiff
#               BC_BN_SignificantASVlist.csv
#               BC_BR_SignificantASVlist.csv
#               UBC_UBN_SignificantASVlist.csv
#               UBC_UBR_SignificantASVlist.csv
#               Microbial_DifferentialAbundance.tiff
#
# PROJECT:      Belowground Plot Experiment (KNZ LTER) fertilization legacies
#
# DATE:         June 2026

# This script is used for the microbial community composition analysis of the 
# project. The output object from the InitialMicrobialAnalysisStep.R file 
# (phyloseq.rds) will be used for this script.

# After uploading the phyloseq object, the first step is to normalize the data. 
# Two normalization approaches are used, each appropriate for the analysis.
# The first normalization is through rarefaction for alpha diversity, and the 
# second is proportional transformation for beta diversity. 

# After normalizing data, the script then details the analysis, first with beta
# diversity and ending with alpha diversity (in order of results in manuscript).

rm(list=ls()) # clear environment

# Set working directory
setwd("~/Desktop/Research/BGPE/MCCRecovery") # change as necessary

# Install R packages if not already
BiocManager::install("DESeq2") # used for differential abundance analysis


# Load packages
library(tidyverse) # v2.0.0
library(phyloseq) # v0.99.6
library(DESeq2) # v1.44.0
library(vegan) # v2.7-1
library(lmerTest) # v3.1-3
library(emmeans) # v1.11.2
library(cowplot) # v1.2.0
library(ggpubr) # v0.6.1


#Bring in data
physeq <- readRDS("phyloseq.rds") #file can be found on github


# The first step of normalization via rarifying is to determine read depth per 
# sample. Because the number of reads are not constant across samples, we have
# to make a call on the number of sequence reads to pull from each sample. 
# However, the tradeoff is that we will lose some samples if they fall below the
# threshold value. 


# First, check to make sure whether there are any ASVs that are included in 
# physeq that are not found in any of the samples.
any(taxa_sums(physeq) == 0) # FALSE - good thing



# Next, check the number of sequence reads for the amplicon sequence variants 
# (ASVs) and samples.

# Code largely borrowed from Paul McMurdie's tutorial on phyloseq
# https://joey711.github.io/phyloseq/

readsumsdf = data.frame(nreads = sort(taxa_sums(physeq), TRUE), 
                        sorted = 1:ntaxa(physeq), 
                        type = "OTUs")

sum(readsumsdf$nreads) # Total number of sequences in library = 9279966

nrow(readsumsdf) # Total number of ASVs = 23486

readsumsdf = rbind(readsumsdf, 
                   data.frame(nreads = sort(sample_sums(physeq), TRUE),
                              sorted = 1:nsamples(physeq),
                              type = "Samples"))


title = "Total number of reads"

p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + 
  geom_bar(stat = "identity")

tiff(file="TotalNumberofReads_ASVsSamples.tiff", width = 7.5, height = 4, 
     pointsize = 1/300, units = 'in', res = 300)
p + ggtitle(title) +  
  scale_y_log10() + 
  annotation_logticks()+ 
  facet_wrap(~type, 1, scales = "free")
dev.off()

# The left panel shows the number of sequence reads for each ASV (OTU in figure),
# and the right panel shows the number of reads per sample

# The pattern in the left panel is expected, where few ASVs are a very abundant
# in terms of sequence reads, with the reads per ASV falling precipitously.
# This pattern is analogous to rank-abundance curves.

# The right panel shows that the number of reads is fairly consistent across 
# samples until the very end where they suddenly drop. This tells us that these
# particularly samples are likely poor in quality (likely poor PCR results). 

# Where the reads start to drop is the target gor sequence depth.

# The next set of script helps us to visualize the sequence reads easier.

data.frame(nreads = sort(sample_sums(physeq),TRUE), 
           sorted = 1:nsamples(physeq), type = "Samples") %>% 
  ggplot(aes(x=1, y=nreads))+
  geom_jitter(shape = 21, size = 3, fill="gray")+
  scale_y_log10()+
  labs(y="Sample sequence depth", x="Samples")+
  theme_bw()+
  theme(axis.text.y = element_text(colour="black", size=10),
        axis.text.x = element_blank(),
        axis.title.x= element_text(colour="black", size=10),
        panel.background = element_rect(fill = "white", color = "black", size = 1))
# In the figure from the above script, it becomes very clear that there are 15
# samples that fall below the cluster. This gives us a clue that the targeted 
# sequence read depth show be above 1e+04 (or 10,000 reads).


data.frame(nreads = sort(sample_sums(physeq),TRUE), 
           sorted = 1:nsamples(physeq), type = "Samples") %>% 
  arrange(nreads) %>% 
  ggplot(aes(x=1:nrow(.), y=nreads))+
  geom_line() +
  theme_bw()
# In this figure, it is easier to see the variability in sequence reads per sample.
# A few samples have much greater number of reads than the rest, and a few samples
# have fewer reads than the others.

# This area where it drops is the target zone in determining the read depth.

# Zoom in this area.
data.frame(nreads = sort(sample_sums(physeq),TRUE), 
           sorted = 1:nsamples(physeq), type = "Samples") %>% 
  arrange(nreads) %>% 
  ggplot(aes(x = 1:nrow(.), y = nreads))+
  geom_line()+
  geom_point(shape = 21, size = 3, fill="gray")+
  theme_bw()+
  coord_cartesian(xlim = c(0, 30),
                  ylim = c(0, 25000))
# Here, it's easier to see where this drop begins. The x-axis is sample order 
# arranged by increasing sequence depth. At sample 20, the "cliff" begins.


# Check the number of sequences in the smallest 25 samples.

data.frame(nreads = sort(sample_sums(physeq),TRUE), 
           sorted = 1:nsamples(physeq), type = "Samples") %>% 
  arrange(nreads) %>% 
  head(.,25)

# Based on this and the figures, the cut-off will be 19000


# The next series of code will create a nice looking rarefaction curve, which
# further illustrates that the sampling depth (19000) is capturing most diversity.

# Remember that DADA2 is being used to denoise data 
# DOI: 10.1038/nmeth.3869

# Apologies in advance, this is quite messy and might take time to run!
tab <- otu_table(physeq)
class(tab) <- "matrix" # as.matrix() will do nothing
## you get a warning here, but this is what we need to have
tab <- t(tab) # transpose observations to rows
out <- rarecurve(tab, step=100, lwd=2, ylab="ASV",  label=F)

names(out) = rownames(tab)

# Coerce data into "long" form.
protox <- mapply(FUN = function(x, y) {
  mydf <- as.data.frame(x)
  colnames(mydf) <- "value"
  mydf$SampleID <- y
  mydf$subsample <- attr(x, "Subsample")
  mydf
}, x = out, y = as.list(names(out)), SIMPLIFY = FALSE)

xy <- do.call(rbind, protox)
rownames(xy) <- NULL  
xy = data.frame(xy, 
                sample_data(physeq)[match(xy$SampleID, 
                                          rownames(sample_data(physeq))), ])

xy$Fertilization <- factor(xy$Fertilization, 
                           levels = c("Control","Recovering", "Fertilized"))

xy$Fire <- factor(xy$Fire, levels = c("Burned", "Unburned"))
levels(xy$Fire)[levels(xy$Fire)=="Burned"] <- "Burned prairie"
levels(xy$Fire)[levels(xy$Fire)=="Unburned"] <- "Unburned prairie"
xy$Fire <- factor(xy$Fire, levels = c("Burned prairie", "Unburned prairie"))


# Plot Rarefaction curve
ggplot(xy, aes(x = subsample, y = value, group = SampleID, color = Fertilization)) +
  scale_color_brewer(palette="Dark2")+
  geom_line() +
  xlim(NA,25000)+
  geom_vline(xintercept=19000, color= "black", linetype='dashed') + 
  labs(x="Sequenced Reads", y="ASVs Detected")+
  facet_grid(Year ~ Fire)+
  theme_bw()+
  theme(axis.text.y = element_text(colour="black", size=10),
        axis.text.x = element_text(colour="black", size=10),
        #axis.title.x=element_blank(),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        legend.title = element_blank())



##### Normalization #####

bgpePrune = prune_samples(sample_sums(physeq)>=19000, physeq)

# Check if there are any ASVs included in dataset but have no reads
any(taxa_sums(bgpePrune) == 0) # TRUE - get rid of these ASVs

# This line removes ASVs if there are no sequence reads
bgpePrune = prune_taxa(taxa_sums(bgpePrune) > 0, bgpePrune)

any(taxa_sums(bgpePrune) == 0) # FALSE - now continue


# New number of reads and ASVs with samples having >= 19000 reads
readsumsdf = data.frame(nreads = sort(taxa_sums(bgpePrune), TRUE), sorted = 1:ntaxa(bgpePrune), 
                        type = "OTUs")
sum(readsumsdf$nreads) # Total number of sequences in library = 9154830
# Loss of 125136 sequences reads

nrow(readsumsdf) # Total number of ASVs = 23371
# Loss of 115 ASVs

###### Rarefied #####
# Rarefied dataset to analyze alpha diversity 
bgpeAlpha = rarefy_even_depth(bgpePrune, replace=TRUE, 
                              rngseed = 20230413,
                              sample.size = 19000)

readsumsdf = data.frame(nreads = sort(taxa_sums(bgpeAlpha), TRUE), sorted = 1:ntaxa(bgpeAlpha), 
                        type = "OTUs")
sum(readsumsdf$nreads) # Total number of sequences in library = 4104000

nrow(readsumsdf) # Total number of ASVs = 21867

###### Proportional transformation #####
bgpeBeta = transform_sample_counts(bgpePrune, 
                                   function(x) 10000 * x/sum(x))

otu_table(bgpeBeta) # ASVs are rows, samples are columns. Common matrix.

readsumsdf = data.frame(nreads = sort(taxa_sums(bgpeBeta), TRUE), sorted = 1:ntaxa(bgpeBeta), 
                        type = "OTUs")
sum(readsumsdf$nreads) # Total number of sequences in library = 2160000

nrow(readsumsdf) # Total number of ASVs = 23371

sample_sums(bgpeBeta) # Each sample has 10,000 "reads" (remember, proportional transformation)


##### Preparing microbial dataset for beta diversity analysis #####

# This is another doozy of a script.
# - This is similar to the plant community analysis, where ASV abundance will be
# averaged by treatment. Now, with this in mind, the script does the following:

# 1) Pivots data to long format to average ASV abundance by field treatment and month
# 2) Returns data to "wide" format to matrix form.
# 3) Re-import averaged abundance data as a new phyloseq object.

# This takes a bit of time to run - fyi
# Step 1
bgpe_average <- bgpeBeta %>% 
  psmelt() %>% # This is equivalent to tidyr::pivot_longer with each row being the
  # abundance of a single ASV per sample
  dplyr::group_by(OTU, Month, Year, Fire, Fertilization) %>% 
  summarize(mean.abundance = mean(Abundance))

# Create a new column that consolidates columns into single metadata column
bgpe_average <- bgpe_average %>% 
  unite("TrtID", c("Fire", "Fertilization", "Month", "Year"), 
        sep="_", remove = T)

# Step 2
ASVmatrix <- bgpe_average %>%
  pivot_wider(names_from = OTU, 
              values_from = mean.abundance) %>% 
  tibble::column_to_rownames("TrtID")
# Rows are samples, columns are ASVs

# Step 3
ASV	=	phyloseq::otu_table(as.matrix(ASVmatrix),	taxa_are_rows	=	FALSE) # phyloseq object

# This code is a bit funky, but what it does is gather the different factors.
# Afterwards, the factors will be concatenated so it matches the sample names 
# in the new ASV matrix table.
bgpe_metaAverage <- bgpeBeta %>% 
  sample_data() %>% 
  as_tibble() %>% 
  add_column(number = 1) %>% 
  group_by(Fire, Fertilization, Month, as.factor(Year)) %>% 
  summarize(mean.number = mean(number))

# Unfortunately, we have a column named 'as.factor(Year)'. This is less than
# ideal, so we will change the name
bgpe_metaAverage$Year <- as.factor(bgpe_metaAverage$"as.factor(Year)")

# Now, create a new meta data for the new phyloseq object
MetaASV <- bgpe_metaAverage %>% 
  unite("TrtID", c("Fire", "Fertilization", "Month", "Year"), 
        sep = "_", remove = FALSE) %>% 
  dplyr::select(-c(mean.number, "as.factor(Year)")) %>% 
  as.data.frame() %>% 
  column_to_rownames("TrtID")

META_ASV = phyloseq::sample_data(MetaASV) # phyloseq object

# Extract taxonomy table from original physeq object
TAX <- tax_table(physeq)

# Combine groups to make a phyloseq object
physeq2 = phyloseq(ASV,	META_ASV, TAX)

# Check to make sure there are no lingering ASVs with no sequence reads
any(taxa_sums(physeq2) == 0) # FALSE - proceed with analysis

##### Non-metric multidimensional scaling (NMDS) #####

# Using "ordinate" function in phyloseq package for NMDS
set.seed(20230413)
microbial.NMDS = ordinate(physeq2, "NMDS", "bray") # stress = 0.186

plot(microbial.NMDS) # base R graph to quickly check plot
# Not super helpful of a visual, just a heads up

###### NMDS publication figure (Figure 3) #####

# Much of the script is tidying the data and to make sure the figure is
# generated correctly.

nmds_microbial <- microbial.NMDS$points %>%
  data.frame() %>% 
  tibble::rownames_to_column(var = "ID") %>% 
  separate(ID, into = c("Fire", "Fertilization",
                        "Month", "Year"), sep = "_") %>% 
  mutate(TrtCombo = case_when(
    Fire == "Burned" & Fertilization == "Control" ~ "BC",
    Fire == "Burned" & Fertilization == "Recovering" ~ "BR",
    Fire == "Burned" & Fertilization == "Fertilized" ~ "BN",
    Fire == "Unburned" & Fertilization == "Control" ~ "UBC",
    Fire == "Unburned" & Fertilization == "Recovering" ~ "UBR",
    Fire == "Unburned" & Fertilization == "Fertilized" ~ "UBN"
  ))

nmds_microbial$TrtCombo <- as.factor(nmds_microbial$TrtCombo)
nmds_microbial$TrtCombo <- factor(nmds_microbial$TrtCombo,
                                  levels = c("BC","BR", "BN", 
                                             "UBC", "UBR", "UBN"))


nmds_microbial$Fire <- as.factor(nmds_microbial$Fire)
levels(nmds_microbial$Fire)[levels(nmds_microbial$Fire)=="Unburned"] <- "Unburned prairie"
levels(nmds_microbial$Fire)[levels(nmds_microbial$Fire)=="Burned"] <- "Burned prairie"
nmds_microbial$Fire <- factor(nmds_microbial$Fire, 
                              levels = c("Burned prairie", "Unburned prairie"))

nmds_microbial$Fertilization <- as.factor(nmds_microbial$Fertilization)
nmds_microbial$Fertilization <- factor(nmds_microbial$Fertilization, 
                                 levels = c("Control", "Recovering", "Fertilized"))

nmds_microbial$Month <- as.factor(nmds_microbial$Month)
nmds_microbial$Month <- factor(nmds_microbial$Month, 
                         levels = c("June", "August"))


nmds_microbial$Year <- as.factor(nmds_microbial$Year)
nmds_microbial$Year <- factor(nmds_microbial$Year, 
                        levels = c("2017", "2018", "2019", "2020", "2021"))


nmds_microbial <- nmds_microbial %>% 
  unite("Grouping2",c("Fertilization","Month","Year"), 
        sep="", remove = FALSE) %>% 
  unite("GroupingTime",c("Month","Year"), 
        sep="", remove = FALSE) %>% 
  mutate(Order = case_when(
    GroupingTime == "June2017" ~ 1,
    GroupingTime == "August2017" ~ 2,
    GroupingTime == "June2018" ~ 3,
    GroupingTime == "August2018" ~ 4,
    GroupingTime == "June2019" ~ 5,
    GroupingTime == "August2019" ~ 6,
    GroupingTime == "June2020" ~ 7,
    GroupingTime == "August2020" ~ 8,
    GroupingTime == "June2021" ~ 9,
    GroupingTime == "August2021" ~ 10
  ))


nmds_microbial$Grouping2 <- as.factor(nmds_microbial$Grouping2)
nmds_microbial$Grouping2 <- factor(nmds_microbial$Grouping2, 
                                   levels = c(
                                     "ControlJune2017", "ControlAugust2017",
                                     "ControlJune2018", "ControlAugust2018",
                                     "ControlJune2019", "ControlAugust2019",
                                     "ControlJune2020", "ControlAugust2020",
                                     "ControlJune2021", "ControlAugust2021",
                                     "RecoveringJune2017", "RecoveringAugust2017",
                                     "RecoveringJune2018", "RecoveringAugust2018",
                                     "RecoveringJune2019", "RecoveringAugust2019",
                                     "RecoveringJune2020", "RecoveringAugust2020",
                                     "RecoveringJune2021", "RecoveringAugust2021",
                                     "FertilizedJune2017", "FertilizedAugust2017",
                                     "FertilizedJune2018", "FertilizedAugust2018",
                                     "FertilizedJune2019", "FertilizedAugust2019",
                                     "FertilizedJune2020", "FertilizedAugust2020",
                                     "FertilizedJune2021", "FertilizedAugust2021"
                                   )
)

nmds_microbial$GroupingTime <- as.factor(nmds_microbial$GroupingTime)

groupingtime_list <- c("June2017","August2017",
                       "June2018","August2018",
                       "June2019","August2019",
                       "June2020","August2020",
                       "June2021","August2021")

nmds_microbial$GroupingTime <- factor(nmds_microbial$GroupingTime, 
                                      levels = groupingtime_list)


Grouping2Color = c("#EDF5F1", "#DBECE5", "#C9E1D7", "#B7D7CA", "#A4CDBC",
                   "#93C4AF", "#80B9A2", "#6FAF94","#5CA686", "#4B9C79",
                   
                   "#FAF0E9", "#F4E0D4", "#EFD1BE", "#EAC2A9", "#E4B393",
                   "#DFA37E", "#DA9468", "#D58553", "#D0763D", "#CA6627",
                   
                   "#F1F0F7", "#E3E2EF", "#D5D4E7", "#C8C6DF", "#ABA9CE",
                   "#ABA9CE", "#9D9BC7", "#908DBE", "#827EB6", "#7470AE")

MCCRecoveryshape <- c(21,22,24,21,22,24)

tiff(file="Microbial_NMDS.tiff",width = 7, height = 4, pointsize = 1/300, units = 'in', res =
       300)
ggplot(data = nmds_microbial[order(nmds_microbial$Order),], 
       aes(x=MDS1, y=MDS2, group = TrtCombo, 
           color=Fertilization, fill=Grouping2, shape=Fertilization))+
  geom_path(aes(colour = Fertilization), linewidth =0.75)+
  geom_point(aes(fill = Grouping2), size = 4, color = "black")+
  scale_shape_manual(values=MCCRecoveryshape)+
  scale_color_brewer(palette="Dark2")+
  scale_fill_manual(values=Grouping2Color)+
  labs(x= "NMDS 1", y="NMDS 2")+
  facet_grid(~Fire)+
  theme_bw()+
  theme(#text = element_text(family = "sans", color = "black"),
    axis.text.y = element_text(colour="black", size=12),
    axis.text.x = element_text(colour="black", size=12),
    strip.text.x = element_text(size = 12, colour = "black"),
    panel.background = element_rect(fill = "white", color = "black", size = 1),
    legend.title = element_blank(),
    legend.position = "none")
dev.off()

##### Permutational analysis of variation (PERMANOVA) #####

# Using adonis2 function for PERMANOVA

df = as(sample_data(physeq2), "data.frame")

# Next three lines makes sure that environmental data are factors
df$Fire <- as.factor(df$Fire)
df$Fertilization <- as.factor(df$Fertilization)
df$Year <- as.factor(df$Year)

# Compute Bray-Curtis dissimilarity matrix
bray_distance.microbial = phyloseq::distance(physeq2, method = "bray")

# Because of the permutational nature of analysis, it is necessary to set a seed
# for reproducibility 
set.seed(20230413)
Microbial_adonis_v1 = adonis2(bray_distance.microbial ~ Fire*Fertilization*Year,
                              by="terms", df)

Microbial_adonis_v1


##### Alpha diversity for microbial community analysis #####

# Make sure to follow this code in order!

Alpha <- bgpeAlpha %>% 
  estimate_richness(split = TRUE, measures = c("Observed", "Shannon")) %>% 
  mutate(Evenness = Shannon / log(Observed)) # alpha metrics

# The Alpha object currently has the row names with "X" in front of it
# because they're numbers, so we will change this to make it easier for future.

row.names(Alpha) <- rownames(sample_data(bgpeAlpha)) # change rownames

# This will change the dataframe into a tibble, which will be easier to add
# meta data to the alpha diversity data set.
Alpha <- Alpha %>% 
  rownames_to_column(var = "ID") %>% # add column
  as_tibble() # change to tibble

# Gather metadata
AlphaMeta <- bgpeAlpha %>% 
  sample_data() %>% 
  data.frame()

row.names(AlphaMeta) <- rownames(sample_data(bgpeAlpha)) # change rownames

AlphaMeta <- AlphaMeta %>% 
  rownames_to_column(var = "ID") %>% 
  as_tibble()

ExperimentalPlotSetup <- data.frame(BlockNumber = c("B1", "B1", "B1", # In manuscript,
                                                    "B2", "B2", "B2", # this is the 
                                                    "B3", "B3", "B3", # plot
                                                    "B4", "B4", "B4",
                                                    "B5", "B5", "B5",
                                                    "B6", "B6", "B6",
                                                    "B7", "B7", "B7",
                                                    "B8", "B8", "B8"),
                                    TrtCombo = c("UBR", "UBC", "UB+N",
                                                 "BC", "BR", "B+N",
                                                 "UBR", "UBC", "UB+N",
                                                 "BC", "BR", "B+N",
                                                 "UBC", "UBR", "UB+N",
                                                 "BR", "BC", "B+N",
                                                 "BR", "BC", "B+N",
                                                 "UBC", "UBR", "UB+N"),
                                    OldPlot = c("P01", "P02", "P01", # b1 - subplot
                                                "P13", "P16", "P16", # b2
                                                "P18", "P20", "P18", # b3
                                                "P25", "P28", "P28", # b4
                                                "P38", "P40", "P40", # b5
                                                "P41", "P42", "P41", # b6
                                                "P54", "P56", "P54", # b7
                                                "P57", "P60", "P60"), # b8
                                    PlotNumber = c("P01", "P02", "P01F", # b1 - small plot
                                                   "P13", "P16", "P16F", # b2 for +N
                                                   "P18", "P20", "P18F", # b3
                                                   "P25", "P28", "P28F", # b4
                                                   "P38", "P40", "P40F", # b5
                                                   "P41", "P42", "P41F", # b6
                                                   "P54", "P56", "P54F", # b7
                                                   "P57", "P60", "P60F"))

# Merge metadata with experiment design
AlphaMeta <- left_join(AlphaMeta, 
                       ExperimentalPlotSetup, 
                       by = c("TrtCombo", "BlockNumber"))

# Join alpha diversity with updated metadata
AlphaDataframe <- left_join(Alpha, AlphaMeta, by = "ID")

# Fot future purposes to work with data out of R, if desired
write.csv(AlphaDataframe, "BPG_MicrobialAlphaDiversity.csv")

# Adding day of year (DOY) to account for repeated sampling within a year
SampleDates <- data.frame(Year = c(2017, 2017,
                                   2018, 2018,
                                   2019, 2019,
                                   2020, 2020,
                                   2021, 2021),
                          Month = c("June","August",
                                    "June","August",
                                    "June","August",
                                    "June","August",
                                    "June","August"),
                          DOY = c(160, 230,
                                  157, 225,
                                  156, 232,
                                  163, 225,
                                  161, 223))


AlphaDataframe <- AlphaDataframe %>% 
  left_join(., SampleDates, by = c("Year", "Month")) %>% 
  mutate(doy.gmc = DOY - mean(DOY)) # day of year scaled to zero

# Before running linear mixed models, make sure fixed effects are factors
AlphaDataframe$Fire <- as.factor(AlphaDataframe$Fire)
AlphaDataframe$Fertilization <- as.factor(AlphaDataframe$Fertilization)
AlphaDataframe$Year <- as.factor(AlphaDataframe$Year)
AlphaDataframe$BlockNumber <- as.factor(AlphaDataframe$BlockNumber)
AlphaDataframe$OldPlot <- as.factor(AlphaDataframe$OldPlot)
AlphaDataframe$PlotNumber <- as.factor(AlphaDataframe$PlotNumber)


# making sure factor levels are in order
AlphaDataframe$Fire <- factor(AlphaDataframe$Fire, 
                              levels = c("Burned", "Unburned"))
AlphaDataframe$Fertilization <- factor(AlphaDataframe$Fertilization, 
                                       levels = c("Control", "Recovering", "Fertilized"))
AlphaDataframe$Year <- factor(AlphaDataframe$Year,
                              levels = c("2017", "2018", 
                                         "2019", "2020", 
                                         "2021"))
AlphaDataframe$BlockNumber <- factor(AlphaDataframe$BlockNumber, 
                                     levels = c("B1", "B2", "B3", "B4",
                                                "B5", "B6", "B7", "B8"))
AlphaDataframe$OldPlot <- factor(AlphaDataframe$OldPlot,
                                 levels = c("P01", "P02", # b1 - subplot
                                             "P13", "P16", # b2
                                             "P18", "P20", # b3
                                             "P25", "P28", # b4
                                             "P38", "P40", # b5
                                             "P41", "P42", # b6
                                             "P54", "P56", # b7
                                             "P57", "P60")) # b8
AlphaDataframe$PlotNumber <- factor(AlphaDataframe$PlotNumber,
                                    levels = c("P01", "P02", "P01F", # b1 - small plot for +N
                                               "P13", "P16", "P16F", # b2 
                                               "P18", "P20", "P18F", # b3
                                               "P25", "P28", "P28F", # b4
                                               "P38", "P40", "P40F", # b5
                                               "P41", "P42", "P41F", # b6
                                               "P54", "P56", "P54F", # b7
                                               "P57", "P60", "P60F")) # b8

###### Richness (Observed) #####
hist(AlphaDataframe$Observed) # Normally distributed

lm.Richness <- lmer(Observed ~ Fire*Fertilization*Year + (1|BlockNumber) + 
                      (1|doy.gmc) + (1|OldPlot:PlotNumber),
                    data = AlphaDataframe, contrasts=contr.sum, REML = TRUE,
                    na.action = na.exclude)

plot(lm.Richness) # appears homoscedastic
qqnorm(resid(lm.Richness)) 
qqline(resid(lm.Richness)) # QQ-plot looks OK
hist(resid(lm.Richness)) # residuals follow normal distribution, small right-tail

anova(lm.Richness)

emm <- emmeans(lm.Richness, pairwise~Fertilization|Year, adjust = "tukey")
multcomp::cld(emm, Letters = letters, reversed = T)
# use multcomp:: because it masks the select function in dplyr




###### Shannon #####
hist(AlphaDataframe$Shannon) # Follows normal distribution

lm.Shannon <- lmer(Shannon ~ Fire*Fertilization*Year + (1|BlockNumber) + 
                     (1|doy.gmc) + (1|OldPlot:PlotNumber),
                   data = AlphaDataframe, contrasts=contr.sum, REML = TRUE,
                   na.action = na.exclude)

plot(lm.Shannon) # appears homoscedastic
qqnorm(resid(lm.Shannon)) 
qqline(resid(lm.Shannon)) # QQ-plot does not look great
hist(resid(lm.Shannon)) # residuals follow normal distribution

anova(lm.Shannon)

emm <- emmeans(lm.Shannon, pairwise~Fertilization|Year, adjust = "tukey")
multcomp::cld(emm, Letters = letters, reversed = T)

###### Evenness #####
hist(AlphaDataframe$Evenness) # Follows normal distribution

lm.Evenness <- lmer(Evenness ~ Fire*Fertilization*Year + (1|BlockNumber) + 
                     (1|doy.gmc) + (1|OldPlot:PlotNumber),
                   data = AlphaDataframe, contrasts=contr.sum, REML = TRUE,
                   na.action = na.exclude)

plot(lm.Evenness) # appears homoscedastic
qqnorm(resid(lm.Evenness)) 
qqline(resid(lm.Evenness)) # QQ-plot does not look great
hist(resid(lm.Evenness)) # residuals follow normal distribution

anova(lm.Evenness)

emm <- emmeans(lm.Evenness, pairwise~Fertilization|Year, adjust = "tukey")
multcomp::cld(emm, Letters = letters, reversed = T)



###### Microbial alpha diversity line graph (Supplemental Figure 4) #####
df.AlphaGraph <- AlphaDataframe %>%
  gather(key = Metric, value = Value, Observed, Shannon, Evenness) %>% 
  group_by(Metric, Fire, Fertilization, Year) %>%
  summarize(mean.metric = mean(Value, na.rm = T), 
            se.metric = plotrix::std.error(Value)) %>% 
  mutate(Treatment = case_when(
    Fire == "Burned" & Fertilization == "Control" ~ "BC",
    Fire == "Burned" & Fertilization == "Recovering" ~ "BR",
    Fire == "Burned" & Fertilization == "Fertilized" ~ "BN",
    Fire == "Unburned" & Fertilization == "Control" ~ "UBC",
    Fire == "Unburned" & Fertilization == "Recovering" ~ "UBR",
    Fire == "Unburned" & Fertilization == "Fertilized" ~ "UBN"
  ))

# Re-making sure data are formatted correctly.
df.AlphaGraph$Treatment <- as.factor(df.AlphaGraph$Treatment)
df.AlphaGraph$Treatment <- factor(df.AlphaGraph$Treatment, 
                                  levels = c("BC", "BR", "BN", 
                                             "UBC", "UBR", "UBN"))

df.AlphaGraph$Metric <- as.factor(df.AlphaGraph$Metric)
df.AlphaGraph$Metric <- factor(df.AlphaGraph$Metric, 
                               levels = c("Observed", 
                                          "Shannon", 
                                          "Evenness"))


df.AlphaGraph$Fire <- as.factor(df.AlphaGraph$Fire)
levels(df.AlphaGraph$Fire)[levels(df.AlphaGraph$Fire)=="Unburned"] <- "Unburned prairie"
levels(df.AlphaGraph$Fire)[levels(df.AlphaGraph$Fire)=="Burned"] <- "Burned prairie"
df.AlphaGraph$Fire <- factor(df.AlphaGraph$Fire, 
                             levels = c("Burned prairie",
                                        "Unburned prairie"))

df.AlphaGraph$Fertilization <- as.factor(df.AlphaGraph$Fertilization)
df.AlphaGraph$Fertilization <- factor(df.AlphaGraph$Fertilization, 
                                      levels = c("Control", 
                                                 "Recovering", 
                                                 "Fertilized"))

df.AlphaGraph$Year <- as.factor(df.AlphaGraph$Year)
df.AlphaGraph$Year <- factor(df.AlphaGraph$Year,
                             levels = c("2017", "2018",
                                        "2019", "2020",
                                        "2021"))


MCCRecoverycolor <- c("#1b9e77", "#d95f02", "#7570b3",
                      "#1b9e77", "#d95f02", "#7570b3")
MCCRecoveryshape <- c(21, 22, 24,
                      21, 22, 24)

# Highlight lines from tiff() to dev.off() to recreate figure.
tiff(file="Alpha_MicrobialRecoveryLineGraph.tiff", width = 6.5, height = 7.5, 
     pointsize = 1/300, units = 'in', res = 300)
ggplot(df.AlphaGraph, aes(x = Year, y = mean.metric, group = Treatment))+
  geom_line(data = df.AlphaGraph[!is.na(df.AlphaGraph$mean.metric),], 
            aes(colour = Treatment), size =0.75)+
  geom_errorbar(aes(ymin = mean.metric-se.metric, 
                    ymax = mean.metric+se.metric), 
                colour = "black", width = .1)+
  geom_point(aes(fill = Treatment, shape = Treatment), 
             width = 0.15, size = 4, color = "black")+
  scale_fill_manual(values = MCCRecoverycolor)+
  scale_colour_manual(values = MCCRecoverycolor)+
  scale_shape_manual(values = MCCRecoveryshape)+
  theme_bw()+
  labs(x = "Year")+
  facet_grid(Metric ~ Fire, scales = "free")+
  theme_bw()+
  theme(axis.text.y = element_text(colour="black", size=10),
        axis.text.x = element_text(colour="black", size=10),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        legend.title = element_blank(),
        legend.position = "none")
dev.off()
# Figure S4 was further annotated in Microsoft Powerpoint


##### Functional group analysis #####

# The next set of script is a bit much, but it essentially classifies taxa
# as different functional groups, if the meet criteria, and then merges
# different objects into one data set.

# We will use the non-rarified dataset.

# Run this code in order

# First, copiotrophs and oligotrophs

ProteobacteriaClass <- c("Alphaproteobacteria", 
                         "Gammaproteobacteria") # Proteobacteria classes

# Because the Alphaproteobacteria and Gammaproteobacteria are classes,
# we will have to separate them as a different object
CopioOligoRelativeAbundance_noProteo <- bgpeBeta %>% 
  subset_taxa(Phylum %in% c("Acidobacteriota", "Actinobacteriota",
                            "Planctomycetota", "Verrucomicrobiota", 
                            "Crenarchaeota", "Myxococcota"))

CopioOligoRelativeAbundance_Proteo <- bgpeBeta %>% 
  subset_taxa(Class %in% c("Alphaproteobacteria", "Gammaproteobacteria"))

CopioOligoRelativeAbundance <- merge_phyloseq(CopioOligoRelativeAbundance_noProteo, 
                                              CopioOligoRelativeAbundance_Proteo)

# Now that the Alphaproteobacteria and Gammaproteobacteria are included,
# and are the only Proteobacteria in the dataset, we can refer just use
# Proteobacteria in the dataset.
CopiotrophList <- c("Crenarchaeota", "Actinobacteriota", 
                    "Proteobacteria") 

OligotrophList <- c("Acidobacteriota", "Planctomycetota", 
                    "Myxococcota", "Verrucomicrobiota")

Copio <- CopioOligoRelativeAbundance %>% 
  subset_taxa(Phylum %in% CopiotrophList) # 5227 taxa

Oligo <- CopioOligoRelativeAbundance %>% 
  subset_taxa(Phylum %in% OligotrophList) # 8001 taxa


# Next, ammonia oxidizing archaea (AOA) and bacteria (AOB)
AOA <- bgpeBeta %>% 
  subset_taxa(Order == "Nitrososphaerales")


AOBList <- c("Nitrosospira", "Nitrosomonas",
             "Nitrospira","Nitrosopumilus",
             "Nitrososphaera", "Nitrosopumilus",
             "Kuenenia", "Nitrococcus") 


AOB <- bgpeBeta %>% 
  subset_taxa(Genus %in% AOBList)


# Now, sum the reads and calculate relative abundance
Copio <- data.frame(Copio %>%
                       sample_sums(.)/10000*100) %>% # Relative abundance (%)
  dplyr::rename(Copiotrophs = 'Copio.....sample_sums....10000...100') %>% 
  rownames_to_column(., var = "SampleID")

Oligo <- data.frame(Oligo %>%
                      sample_sums(.)/10000*100) %>% # Relative abundance (%)
  dplyr::rename(Oligotrophs = 'Oligo.....sample_sums....10000...100') %>% 
  rownames_to_column(., var = "SampleID")

AOA <- data.frame(AOA %>%
                    sample_sums(.)/10000*100) %>% # Relative abundance (%)
  dplyr::rename(AOA = 'AOA.....sample_sums....10000...100') %>% 
  rownames_to_column(., var = "SampleID")

AOB <- data.frame(AOB %>%
                    sample_sums(.)/10000*100) %>% # Relative abundance (%)
  dplyr::rename(AOB = 'AOB.....sample_sums....10000...100') %>% 
  rownames_to_column(., var = "SampleID")

# Meta data 
BetaMeta <- bgpePrune %>%
  sample_data() %>% 
  data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID")

# Join objects
# First joining all relative abundance of functional groups
# Then adding meta data for statistics and graphing purposes

df.FunctionalGroup <- Copio %>% 
  left_join(., Oligo, by = "SampleID") %>% 
  left_join(., AOA, by = "SampleID") %>% 
  left_join(., AOB, by = "SampleID") %>% 
  left_join(., BetaMeta, by = "SampleID") %>% 
  left_join(., ExperimentalPlotSetup, 
            by = c("TrtCombo", "BlockNumber")) %>% 
  dplyr::select(SampleID, Fire, Fertilization, Month, Year,
                TrtCombo, BlockNumber, OldPlot, PlotNumber,
                Copiotrophs, Oligotrophs, AOA, AOB) %>% 
  left_join(., SampleDates, by = c("Year", "Month")) %>% 
  mutate(doy.gmc = DOY - mean(DOY))


# Before running linear mixed models, make sure fixed effects are factors
df.FunctionalGroup$Fire <- as.factor(df.FunctionalGroup$Fire)
df.FunctionalGroup$Fertilization <- as.factor(df.FunctionalGroup$Fertilization)
df.FunctionalGroup$Year <- as.factor(df.FunctionalGroup$Year)
df.FunctionalGroup$BlockNumber <- as.factor(df.FunctionalGroup$BlockNumber)
df.FunctionalGroup$OldPlot <- as.factor(df.FunctionalGroup$OldPlot)
df.FunctionalGroup$PlotNumber <- as.factor(df.FunctionalGroup$PlotNumber)


# making sure factor levels are in order
df.FunctionalGroup$Fire <- factor(df.FunctionalGroup$Fire, 
                              levels = c("Burned", "Unburned"))
df.FunctionalGroup$Fertilization <- factor(df.FunctionalGroup$Fertilization, 
                                       levels = c("Control", "Recovering", "Fertilized"))
df.FunctionalGroup$Year <- factor(df.FunctionalGroup$Year,
                              levels = c("2017", "2018", 
                                         "2019", "2020", 
                                         "2021"))
df.FunctionalGroup$BlockNumber <- factor(df.FunctionalGroup$BlockNumber, 
                                     levels = c("B1", "B2", "B3", "B4",
                                                "B5", "B6", "B7", "B8"))
df.FunctionalGroup$OldPlot <- factor(df.FunctionalGroup$OldPlot,
                                 levels = c("P01", "P02", # b1 - subplot
                                            "P13", "P16", # b2
                                            "P18", "P20", # b3
                                            "P25", "P28", # b4
                                            "P38", "P40", # b5
                                            "P41", "P42", # b6
                                            "P54", "P56", # b7
                                            "P57", "P60")) # b8
df.FunctionalGroup$PlotNumber <- factor(df.FunctionalGroup$PlotNumber,
                                    levels = c("P01", "P02", "P01F", # b1 - small plot for +N
                                               "P13", "P16", "P16F", # b2 
                                               "P18", "P20", "P18F", # b3
                                               "P25", "P28", "P28F", # b4
                                               "P38", "P40", "P40F", # b5
                                               "P41", "P42", "P41F", # b6
                                               "P54", "P56", "P54F", # b7
                                               "P57", "P60", "P60F")) # b8


###### Copiotrophs #####
hist(df.FunctionalGroup$Copiotrophs) # somewhat normally distributed, right-tail

lm.Copio <- lmer(Copiotrophs ~ Fire*Fertilization*Year + (1|BlockNumber) + 
                   (1|doy.gmc) + (1|OldPlot:PlotNumber),
                 data = df.FunctionalGroup, contrasts=contr.sum, REML = TRUE,
                 na.action = na.exclude)

plot(lm.Copio) # appears homoscedastic, but trend on far right side
qqnorm(resid(lm.Copio)) 
qqline(resid(lm.Copio)) # QQ-plot looks OK
hist(resid(lm.Copio)) # long right-tail


# Try log-transformation
lm.Copio <- lmer(log(Copiotrophs+0.1) ~ Fire*Fertilization*Year + (1|BlockNumber) + 
                   (1|doy.gmc) + (1|OldPlot:PlotNumber),
                 data = df.FunctionalGroup, contrasts=contr.sum, REML = TRUE,
                 na.action = na.exclude)

plot(lm.Copio) # appears homoscedastic, pattern on far right side disappeared
qqnorm(resid(lm.Copio)) 
qqline(resid(lm.Copio)) # QQ-plot looks OK
hist(resid(lm.Copio)) # residuals normally distributed

anova(lm.Copio)

emm <- emmeans(lm.Copio, pairwise ~ Fire:Fertilization, adjust = "tukey")
multcomp::cld(emm, Letters = letters, reversed = F)


emm <- emmeans(lm.Copio, pairwise ~ Fertilization|Year, adjust = "tukey")
multcomp::cld(emm, Letters = letters, reversed = T)


###### Oligotrophs #####
hist(df.FunctionalGroup$Oligotrophs) # somewhat normally distributed, left-tail

lm.Oligo <- lmer(Oligotrophs ~ Fire*Fertilization*Year + (1|BlockNumber) + 
                   (1|doy.gmc) + (1|OldPlot:PlotNumber),
                 data = df.FunctionalGroup, contrasts=contr.sum, REML = TRUE,
                 na.action = na.exclude)

plot(lm.Oligo) # appears homoscedastic
qqnorm(resid(lm.Oligo)) 
qqline(resid(lm.Oligo)) # QQ-plot looks OK
hist(resid(lm.Oligo)) # appears normally distributed

anova(lm.Oligo)

emm <- emmeans(lm.Oligo, pairwise ~ Fertilization|Year, adjust = "tukey")
multcomp::cld(emm, Letters = letters, reversed = F)


###### AOA #####
hist(df.FunctionalGroup$AOA) # not normally distributed, right tail

lm.AOA <- lmer(AOA ~ Fire*Fertilization*Year + (1|BlockNumber) + 
                   (1|doy.gmc) + (1|OldPlot:PlotNumber),
               data = df.FunctionalGroup, contrasts=contr.sum, REML = TRUE,
               na.action = na.exclude)

plot(lm.AOA) # mostly homoscedastic; seems to spread out on the right
qqnorm(resid(lm.AOA)) 
qqline(resid(lm.AOA)) # QQ-plot looks OK
hist(resid(lm.AOA)) # appears normally distributed

anova(lm.AOA)

emm <- emmeans(lm.AOA, pairwise ~ Fertilization|Year, adjust = "tukey")
multcomp::cld(emm, Letters = letters, reversed = T)


###### AOB #####
hist(df.FunctionalGroup$AOB) # not normally distributed, right tail

lm.AOB <- lmer(AOB ~ Fire*Fertilization*Year + (1|BlockNumber) + 
                 (1|doy.gmc) + (1|OldPlot:PlotNumber),
               data = df.FunctionalGroup, contrasts=contr.sum, REML = TRUE,
               na.action = na.exclude)

plot(lm.AOB) # sort of homoscedastic; a bit constrained on the left side, but overall fine
qqnorm(resid(lm.AOB)) 
qqline(resid(lm.AOB)) # QQ-plot looks OK
hist(resid(lm.AOB)) # appears normally distributed

anova(lm.AOB)

emm <- emmeans(lm.AOB, pairwise ~ Fertilization|Year, adjust = "tukey")
multcomp::cld(emm, Letters = letters, reversed = T)


###### Microbial functional group publication figure (Figure 4) #####
df.FunctionalGroupGraph <- df.FunctionalGroup %>%
  pivot_longer(
    cols = c(Copiotrophs, Oligotrophs, AOB, AOA),
    names_to = "Lifehistory", 
    values_to = "RelativeAbundance"
  ) %>% 
  group_by(Lifehistory, Fire, Fertilization, Year) %>%
  summarize(mean.RelativeAbundance = mean(RelativeAbundance, na.rm = T), 
            se.RelativeAbundance = plotrix::std.error(RelativeAbundance)) %>% 
  mutate(Treatment = case_when(
    Fire == "Burned" & Fertilization == "Control" ~ "BC",
    Fire == "Burned" & Fertilization == "Recovering" ~ "BR",
    Fire == "Burned" & Fertilization == "Fertilized" ~ "BN",
    Fire == "Unburned" & Fertilization == "Control" ~ "UBC",
    Fire == "Unburned" & Fertilization == "Recovering" ~ "UBR",
    Fire == "Unburned" & Fertilization == "Fertilized" ~ "UBN"
  ))

# Re-making sure data are formatted correctly.
df.FunctionalGroupGraph$Treatment <- as.factor(df.FunctionalGroupGraph$Treatment)
df.FunctionalGroupGraph$Treatment <- factor(df.FunctionalGroupGraph$Treatment, 
                                  levels = c("BC", "BR", "BN", 
                                             "UBC", "UBR", "UBN"))

df.FunctionalGroupGraph$Lifehistory <- as.factor(df.FunctionalGroupGraph$Lifehistory)
df.FunctionalGroupGraph$Lifehistory <- factor(df.FunctionalGroupGraph$Lifehistory, 
                               levels = c("Copiotrophs", 
                                          "Oligotrophs", 
                                          "AOA",
                                          "AOB"))


df.FunctionalGroupGraph$Fire <- as.factor(df.FunctionalGroupGraph$Fire)
levels(df.FunctionalGroupGraph$Fire)[levels(df.FunctionalGroupGraph$Fire)=="Unburned"] <- "Unburned prairie"
levels(df.FunctionalGroupGraph$Fire)[levels(df.FunctionalGroupGraph$Fire)=="Burned"] <- "Burned prairie"
df.FunctionalGroupGraph$Fire <- factor(df.FunctionalGroupGraph$Fire, 
                             levels = c("Burned prairie",
                                        "Unburned prairie"))

df.FunctionalGroupGraph$Fertilization <- as.factor(df.FunctionalGroupGraph$Fertilization)
df.FunctionalGroupGraph$Fertilization <- factor(df.FunctionalGroupGraph$Fertilization, 
                                      levels = c("Control", 
                                                 "Recovering", 
                                                 "Fertilized"))

df.FunctionalGroupGraph$Year <- as.factor(df.FunctionalGroupGraph$Year)
df.FunctionalGroupGraph$Year <- factor(df.FunctionalGroupGraph$Year,
                             levels = c("2017", "2018",
                                        "2019", "2020",
                                        "2021"))


MCCRecoverycolor <- c("#1b9e77", "#d95f02", "#7570b3",
                      "#1b9e77", "#d95f02", "#7570b3")
MCCRecoveryshape <- c(21, 22, 24,
                      21, 22, 24)

# Highlight lines from tiff() to dev.off() to recreate figure.
tiff(file="MicrobialFunctionalGroup_LineGraph.tiff", width = 5, height = 8, 
     pointsize = 1/300, units = 'in', res = 300)
ggplot(df.FunctionalGroupGraph, aes(x = Year, y = mean.RelativeAbundance, 
                                    group = Treatment))+
  geom_errorbar(aes(ymin = mean.RelativeAbundance-se.RelativeAbundance,
                    ymax = mean.RelativeAbundance+se.RelativeAbundance),
                colour="black", width=.1)+
  geom_line(data = df.FunctionalGroupGraph, aes(colour = Treatment), size = 0.65)+
  geom_point(aes(fill = Treatment, shape = Treatment), 
             width = 0.15, size = 4, color = "black")+
  scale_fill_manual(values = MCCRecoverycolor)+
  scale_colour_manual(values = MCCRecoverycolor)+
  scale_shape_manual(values = MCCRecoveryshape)+
  theme_bw()+
  labs(x = "Year",
       y = "Relative abundance (%)")+
  facet_grid(Lifehistory ~ Fire, scales = "free")+
  theme_bw()+
  theme(axis.text.y = element_text(colour="black", size=10),
        axis.text.x = element_text(colour="black", size=10),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        legend.title = element_blank(),
        legend.position = "none")
dev.off()
# Figure 4 was further annotated in Microsoft Powerpoint

# This figure include the legend
ggplot(df.FunctionalGroupGraph, aes(x = Year, y = mean.RelativeAbundance, 
                                               group = Treatment))+
  geom_errorbar(aes(ymin = mean.RelativeAbundance-se.RelativeAbundance,
                    ymax = mean.RelativeAbundance+se.RelativeAbundance),
                colour="black", width=.1)+
  geom_line(data = df.FunctionalGroupGraph, aes(colour = Treatment), size = 0.65)+
  geom_point(aes(fill = Treatment, shape = Treatment), 
             width = 0.15, size = 4, color = "black")+
  scale_fill_manual(values = MCCRecoverycolor)+
  scale_colour_manual(values = MCCRecoverycolor)+
  scale_shape_manual(values = MCCRecoveryshape)+
  labs(x = "Year",
       y = "Relative abundance (%)")+
  facet_grid(Lifehistory ~ Fire, scales = "free")+
  theme_bw()+
  theme(axis.text.y = element_text(colour="black", size=10),
        axis.text.x = element_text(colour="black", size=10),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        legend.title = element_blank())

# THIS SCRIPT IS PURELY TO EXTRACT A LEGEND FOR PUBLICATION
minimal <- df.FunctionalGroupGraph %>% 
  group_by(Fertilization, Year) %>% 
  summarize(mean.RelativeAbundance = mean(mean.RelativeAbundance, na.rm = T)) %>% 
  ggplot(., aes(x = Year, y = mean.RelativeAbundance, 
                                      group = Fertilization))+
  geom_line(aes(colour = Fertilization), size = 0.65)+
  geom_point(aes(fill = Fertilization, shape = Fertilization), 
             width = 0.15, size = 4, color = "black")+
  scale_fill_manual(values = MCCRecoverycolor)+
  scale_colour_manual(values = MCCRecoverycolor)+
  scale_shape_manual(values = MCCRecoveryshape)+
  labs(x = "Year",
       y = "Relative abundance (%)")+
  theme_bw()+
  theme(axis.text.y = element_text(colour="black", size=10),
        axis.text.x = element_text(colour="black", size=10),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        legend.title = element_blank())

leg <- get_legend(minimal)

tiff("FertilizationTreatment_Legend.tiff", 
     width = 4, height = 4, pointsize = 1/300, units = 'in', res = 300)
as_ggplot(leg)
dev.off()


###### Microbial functional group boxplot (Supplemental Figure 5) #####
df.FunctionalGroupGraph <- df.FunctionalGroup %>%
  pivot_longer(
    cols = c(Copiotrophs, Oligotrophs, AOB, AOA),
    names_to = "Lifehistory", 
    values_to = "RelativeAbundance"
  )

# Re-making sure data are formatted correctly.
df.FunctionalGroupGraph$Lifehistory <- as.factor(df.FunctionalGroupGraph$Lifehistory)
df.FunctionalGroupGraph$Lifehistory <- factor(df.FunctionalGroupGraph$Lifehistory, 
                                              levels = c("Copiotrophs", 
                                                         "Oligotrophs", 
                                                         "AOA",
                                                         "AOB"))


df.FunctionalGroupGraph$Fire <- as.factor(df.FunctionalGroupGraph$Fire)
levels(df.FunctionalGroupGraph$Fire)[levels(df.FunctionalGroupGraph$Fire)=="Unburned"] <- "Unburned prairie"
levels(df.FunctionalGroupGraph$Fire)[levels(df.FunctionalGroupGraph$Fire)=="Burned"] <- "Burned prairie"
df.FunctionalGroupGraph$Fire <- factor(df.FunctionalGroupGraph$Fire, 
                                       levels = c("Burned prairie",
                                                  "Unburned prairie"))

df.FunctionalGroupGraph$Fertilization <- as.factor(df.FunctionalGroupGraph$Fertilization)
df.FunctionalGroupGraph$Fertilization <- factor(df.FunctionalGroupGraph$Fertilization, 
                                                levels = c("Control", 
                                                           "Recovering", 
                                                           "Fertilized"))

MCCRecoverycolor <- c("#1b9e77", "#d95f02", "#7570b3")
MCCRecoveryshape <- c(21, 22, 24)

# Highlight lines from tiff() to dev.off() to recreate figure.
tiff(file="MicrobialFunctionalGroup_Boxplot.tiff", width = 6, height = 6, 
     pointsize = 1/300, units = 'in', res = 300)
ggplot(df.FunctionalGroupGraph, aes(x = Fire, y = RelativeAbundance, 
                                    color = Fertilization))+
  geom_boxplot(linewidth = 1, outlier.shape = NA, fill = "white")+
  geom_point(aes(fill = Fertilization, shape = Fertilization), 
             size = 2, position = position_jitterdodge(), alpha = 1/3)+
  scale_colour_manual(values = MCCRecoverycolor)+
  scale_fill_manual(values = MCCRecoverycolor)+
  scale_shape_manual(values = MCCRecoveryshape)+
  labs(y = "Relative Abundance (%)")+
  facet_wrap(~Lifehistory, scales = "free")+
  theme_bw()+
  theme(axis.text.y = element_text(colour="black", size=10),
        axis.text.x = element_text(colour="black", size=10),
        axis.title.x=element_blank(),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        legend.title = element_blank(),
        legend.position = "none")
dev.off()
# Figure S5 was further annotated in Microsoft Powerpoint

minimal <- ggplot(df.FunctionalGroupGraph, aes(x = Fire, y = RelativeAbundance, 
                                               color = Fertilization))+
  geom_boxplot(linewidth = 1, outlier.shape = NA, fill = "white")+
  geom_point(aes(fill = Fertilization, shape = Fertilization), 
             size = 2, position = position_jitterdodge(), alpha = 1/3)+
  scale_colour_manual(values = MCCRecoverycolor)+
  scale_fill_manual(values = MCCRecoverycolor)+
  scale_shape_manual(values = MCCRecoveryshape)+
  labs(y = "Relative Abundance (%)")+
  facet_wrap(~Lifehistory, scales = "free")+
  theme_bw()+
  theme(axis.text.y = element_text(colour="black", size=10),
        axis.text.x = element_text(colour="black", size=10),
        axis.title.x=element_blank(),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        legend.title = element_blank())

leg <- get_legend(minimal)

tiff("MicrobialFunctionalGroup_Boxplot_Legend.tiff", width = 4, height = 4, pointsize = 1/300, units = 'in', res = 300)
as_ggplot(leg)
dev.off()


##### Differential abundance analysis (DESeq2) #####


# The general workflow is the following:
# 1) Create a temporary phyloseq object (bgpeBeta_temp) with treatments of interest
# 2) Perform the differential abundance (DA) analysis with DESeq2
# 3) Select DA ASV taxa at p < 0.01
# 4) Create a CSV file to create a record of the treatment comparison.

# Step 1) DA between BC and B+N
bgpeBeta_temp <- bgpeBeta %>% 
  subset_samples(TrtCombo %in% c("BC", "B+N"))

# Step 2) Perform DESeq2 analysis
burned_DA = phyloseq_to_deseq2(bgpeBeta_temp, ~ Fertilization)
burned_DA = DESeq(burned_DA, test="Wald", fitType="parametric")

# Step 3) Identify differentially abundanct taxa at P < 0.01
res = results(burned_DA, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), 
               as(tax_table(bgpeBeta_temp)[rownames(sigtab), ], "matrix"))

nrow(sigtab) #152 ASVs

# Step 4) Create new object and save as .csv
bc_bn <- sigtab

bc_bn %>% 
  arrange(desc(log2FoldChange)) %>% 
  data.frame() %>% 
  write.csv("BC_BN_SignificantASVlist.csv")



# Step 1) DA between BC and BR
bgpeBeta_temp <- bgpeBeta %>% 
  subset_samples(TrtCombo %in% c("BC", "BR"))

# Step 2)
burned_DA = phyloseq_to_deseq2(bgpeBeta_temp, ~ Fertilization)
burned_DA = DESeq(burned_DA, test="Wald", fitType="parametric")

# Step 3)
res = results(burned_DA, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(bgpeBeta_temp)[rownames(sigtab), ], "matrix"))

nrow(sigtab) #233 ASVs

# Step 4)
bc_br <- sigtab

bc_br %>% 
  arrange(desc(log2FoldChange)) %>% 
  data.frame() %>% 
  write.csv("BC_BR_SignificantASVlist.csv")



# Step 1) DA analysis for UBC and UB+N
bgpeBeta_temp <- bgpeBeta %>% 
  subset_samples(TrtCombo %in% c("UBC", "UB+N"))

# Step 2)
unburned_DA = phyloseq_to_deseq2(bgpeBeta_temp, ~ Fertilization)
unburned_DA = DESeq(unburned_DA, test="Wald", fitType="parametric")

# Step 3)
res = results(unburned_DA, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(bgpeBeta_temp)[rownames(sigtab), ], "matrix"))

nrow(sigtab) #55 ASVs

# Step 4)
ubc_ubn <- sigtab

ubc_ubn %>% 
  arrange(desc(log2FoldChange)) %>% 
  data.frame() %>% 
  write.csv("UBC_UBN_SignificantASVlist.csv")


# Step 1) DA analysis for UBC and UBR
bgpeBeta_temp <- bgpeBeta %>% 
  subset_samples(TrtCombo %in% c("UBC", "UBR"))

# Step 2)
unburned_DA = phyloseq_to_deseq2(bgpeBeta_temp, ~ Fertilization)
unburned_DA = DESeq(unburned_DA, test="Wald", fitType="parametric")

# Step 3)
res = results(unburned_DA, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(bgpeBeta_temp)[rownames(sigtab), ], "matrix"))

nrow(sigtab) #114 ASVs

# Step 4)
ubc_ubr <- sigtab

ubc_ubr %>% 
  arrange(desc(log2FoldChange)) %>% 
  data.frame() %>% 
  write.csv("UBC_UBR_SignificantASVlist.csv")



# Statistics for whether log-fold changes in abundance are different between
# fire and enrichment treatments.

# First, separate positive and negative log-fold changes in abundance.
bc_bn.lfc.positive <- bc_bn %>% 
  dplyr::filter(log2FoldChange > 0) %>% 
  dplyr::select(log2FoldChange) %>% 
  mutate(Fire = "Burned",
         Fert = "Fertilized")

bc_br.lfc.positive <- bc_br %>% 
  dplyr::filter(log2FoldChange > 0) %>% 
  dplyr::select(log2FoldChange) %>% 
  mutate(Fire = "Burned",
         Fert = "Recovering")

ubc_ubn.lfc.positive <- ubc_ubn %>% 
  dplyr::filter(log2FoldChange > 0) %>% 
  dplyr::select(log2FoldChange) %>% 
  mutate(Fire = "Unburned",
         Fert = "Fertilized")

ubc_ubr.lfc.positive <- ubc_ubr %>% 
  dplyr::filter(log2FoldChange > 0) %>% 
  dplyr::select(log2FoldChange) %>% 
  mutate(Fire = "Unburned",
         Fert = "Recovering")

# Combine positive log2FoldChange changes
rbind(bc_bn.lfc.positive,
      bc_br.lfc.positive,
      ubc_ubn.lfc.positive,
      ubc_ubr.lfc.positive) %>% 
  ggplot(., aes(x = Fire, y = log2FoldChange))+
  geom_jitter(shape = 21, size = 2, color="black")+
  theme_bw()
# There is a clear outlier (log2FoldChange > 20) in the unburned prairie
# Let's remove it.

df.lfc.positive <- rbind(bc_bn.lfc.positive,
                         bc_br.lfc.positive,
                         ubc_ubn.lfc.positive,
                         ubc_ubr.lfc.positive) %>% 
  dplyr::filter(log2FoldChange < 10) %>% 
  mutate(Change = "Positive")


# Next, repeat the same steps but for negative changes
bc_bn.lfc.negative <- bc_bn %>% 
  dplyr::filter(log2FoldChange < 0) %>% 
  dplyr::select(log2FoldChange) %>% 
  mutate(Fire = "Burned",
         Fert = "Fertilized")

bc_br.lfc.negative <- bc_br %>% 
  dplyr::filter(log2FoldChange < 0) %>% 
  dplyr::select(log2FoldChange) %>% 
  mutate(Fire = "Burned",
         Fert = "Recovering")


ubc_ubn.lfc.negative <- ubc_ubn %>% 
  dplyr::filter(log2FoldChange < 0) %>% 
  dplyr::select(log2FoldChange) %>% 
  mutate(Fire = "Unburned",
         Fert = "Fertilized")

ubc_ubr.lfc.negative <- ubc_ubr %>% 
  dplyr::filter(log2FoldChange < 0) %>% 
  dplyr::select(log2FoldChange) %>% 
  mutate(Fire = "Unburned",
         Fert = "Recovering")


# Combine negative log2FoldChange changes
rbind(bc_bn.lfc.negative,
      bc_br.lfc.negative,
      ubc_ubn.lfc.negative,
      ubc_ubr.lfc.negative) %>% 
  ggplot(., aes(x = Fire, y = log2FoldChange))+
  geom_jitter(shape = 21, size = 2, color="black")+
  theme_bw()
# No outliers here; no need to remove data.

df.lfc.negative <- rbind(bc_bn.lfc.negative,
                         bc_br.lfc.negative,
                         ubc_ubn.lfc.negative,
                         ubc_ubr.lfc.negative) %>% 
  mutate(Change = "Negative")


###### Positive differential abundance #####
lm.lfc.positive <- lm(log2FoldChange ~ Fire*Fert, 
                      data = df.lfc.positive)

plot(resid(lm.lfc.positive))
abline(h = 0) # appears homoscedastic

qqnorm(resid(lm.lfc.positive)) 
qqline(resid(lm.lfc.positive)) # QQ-plot look OK
hist(resid(lm.lfc.positive)) # appears normally distributed

anova(lm.lfc.positive)

# How many ASVs are differentially abundant?
df.lfc.positive %>% 
  group_by(Fire, Fert) %>% 
  summarize(count = n())
# Remember that one datum was excluded, so UB+N is actually 42


# The following script will be used to fill in Table S2
# Burned recovering - copiotrophs
bc_br.lfc.positive %>% 
  tibble::rownames_to_column("ASV") %>% 
  pull(ASV) %>% 
  prune_taxa(., CopioOligoRelativeAbundance) %>% 
  subset_taxa(Phylum %in% CopiotrophList) # 30 ASVs
30/23371*100 # 0.1283642

# Burned recovering - oligotrophs
bc_br.lfc.positive %>% 
  tibble::rownames_to_column("ASV") %>% 
  pull(ASV) %>% 
  prune_taxa(., CopioOligoRelativeAbundance) %>% 
  subset_taxa(Phylum %in% OligotrophList) # 58 ASVs
58/23371*100 # 0.2481708

# Burned recovering - neither
106-58-30 # 18
18/23371*100 # 0.07701853


# Burned fertilized - copiotrophs
bc_bn.lfc.positive %>% 
  tibble::rownames_to_column("ASV") %>% 
  pull(ASV) %>% 
  prune_taxa(., CopioOligoRelativeAbundance) %>% 
  subset_taxa(Phylum %in% CopiotrophList) # 33 ASVs
33/23371*100 # 0.1412006

# Burned fertilized - oligotrophs
bc_bn.lfc.positive %>% 
  tibble::rownames_to_column("ASV") %>% 
  pull(ASV) %>% 
  prune_taxa(., CopioOligoRelativeAbundance) %>% 
  subset_taxa(Phylum %in% OligotrophList) # 25 ASVs
25/23371*100 # 0.1069702

# Burned fertilized - neither
75-33-25 # 17
17/23371*100 # 0.07273972


# Unburned recovering - copiotrophs
ubc_ubr.lfc.positive %>% 
  tibble::rownames_to_column("ASV") %>% 
  pull(ASV) %>% 
  prune_taxa(., CopioOligoRelativeAbundance) %>% 
  subset_taxa(Phylum %in% CopiotrophList) # 33 ASVs
33/23371*100 # 0.1412006

# Unburned recovering - oligotrophs
ubc_ubr.lfc.positive %>% 
  tibble::rownames_to_column("ASV") %>% 
  pull(ASV) %>% 
  prune_taxa(., CopioOligoRelativeAbundance) %>% 
  subset_taxa(Phylum %in% OligotrophList) # 29 ASVs
29/23371*100 # 0.1240854

# Unburned recovering - neither
71-33-29 # 9
9/23371*100 # 0.03850926


# Unburned fertilized - copiotrophs
ubc_ubn.lfc.positive %>% 
  tibble::rownames_to_column("ASV") %>% 
  pull(ASV) %>% 
  prune_taxa(., CopioOligoRelativeAbundance) %>% 
  subset_taxa(Phylum %in% CopiotrophList) # 25 ASVs
25/23371*100 # 0.1069702

# Unburned fertilized - oligotrophs
ubc_ubn.lfc.positive %>% 
  tibble::rownames_to_column("ASV") %>% 
  pull(ASV) %>% 
  prune_taxa(., CopioOligoRelativeAbundance) %>% 
  subset_taxa(Phylum %in% OligotrophList) # 16 ASVs
16/23371*100 # 0.06846091

# Unburned fertilized - neither
# UB+N is actually 42
42-25-16 # 1
1/23371*100 # 0.004278807





###### Negative differential abundance #####
lm.lfc.negative <- lm(log2FoldChange ~ Fire*Fert, 
                      data = df.lfc.negative)

plot(resid(lm.lfc.negative))
abline(h = 0) # appears homoscedastic

qqnorm(resid(lm.lfc.negative)) 
qqline(resid(lm.lfc.negative)) # QQ-plot look OK
hist(resid(lm.lfc.negative)) # appears normally distributed

anova(lm.lfc.negative)

# How many ASVs are differentially abundant?
df.lfc.negative %>% 
  group_by(Fire, Fert) %>% 
  summarize(count = n())



# How many ASVs are differentially abundant?
df.lfc.negative %>% 
  group_by(Fire, Fert) %>% 
  summarize(count = n())


# The following script will be used to fill in Table S2
# Burned recovering - copiotrophs
bc_br.lfc.negative %>% 
  tibble::rownames_to_column("ASV") %>% 
  pull(ASV) %>% 
  prune_taxa(., CopioOligoRelativeAbundance) %>% 
  subset_taxa(Phylum %in% CopiotrophList) # 47 ASVs
47/23371*100 # 0.2011039

# Burned recovering - oligotrophs
bc_br.lfc.negative %>% 
  tibble::rownames_to_column("ASV") %>% 
  pull(ASV) %>% 
  prune_taxa(., CopioOligoRelativeAbundance) %>% 
  subset_taxa(Phylum %in% OligotrophList) # 56 ASVs
56/23371*100 # 0.2396132

# Burned recovering - neither
127-47-56 # 24
24/23371*100 # 0.1026914


# Burned fertilized - copiotrophs
bc_bn.lfc.negative %>% 
  tibble::rownames_to_column("ASV") %>% 
  pull(ASV) %>% 
  prune_taxa(., CopioOligoRelativeAbundance) %>% 
  subset_taxa(Phylum %in% CopiotrophList) # 18 ASVs
18/23371*100 # 0.07701853

# Burned fertilized - oligotrophs
bc_bn.lfc.negative %>% 
  tibble::rownames_to_column("ASV") %>% 
  pull(ASV) %>% 
  prune_taxa(., CopioOligoRelativeAbundance) %>% 
  subset_taxa(Phylum %in% OligotrophList) # 39 ASVs
39/23371*100 # 0.1069702

# Burned fertilized - neither
77-39-18 # 20
20/23371*100 # 0.08557614


# Unburned recovering - copiotrophs
ubc_ubr.lfc.negative %>% 
  tibble::rownames_to_column("ASV") %>% 
  pull(ASV) %>% 
  prune_taxa(., CopioOligoRelativeAbundance) %>% 
  subset_taxa(Phylum %in% CopiotrophList) # 8 ASVs
8/23371*100 # 0.03423046

# Unburned recovering - oligotrophs
ubc_ubr.lfc.negative %>% 
  tibble::rownames_to_column("ASV") %>% 
  pull(ASV) %>% 
  prune_taxa(., CopioOligoRelativeAbundance) %>% 
  subset_taxa(Phylum %in% OligotrophList) # 25 ASVs
25/23371*100 # 0.1069702

# Unburned recovering - neither
43-8-25 # 10
10/23371*100 # 0.04278807


# Unburned fertilized - copiotrophs
ubc_ubn.lfc.negative %>% 
  tibble::rownames_to_column("ASV") %>% 
  pull(ASV) %>% 
  prune_taxa(., CopioOligoRelativeAbundance) %>% 
  subset_taxa(Phylum %in% CopiotrophList) # 1 ASVs
1/23371*100 # 0.004278807

# Unburned fertilized - oligotrophs
ubc_ubn.lfc.negative %>% 
  tibble::rownames_to_column("ASV") %>% 
  pull(ASV) %>% 
  prune_taxa(., CopioOligoRelativeAbundance) %>% 
  subset_taxa(Phylum %in% OligotrophList) # 9 ASVs
9/23371*100 # 0.03850926

# Unburned fertilized - neither
13-1-9 # 3
3/23371*100 # 0.01283642



###### Differential abundance boxplot (Supplemental Figure 6) #####
df.DifferentialAbundanceGraph <- rbind(df.lfc.positive, df.lfc.negative)

# Re-making sure data are formatted correctly.
df.DifferentialAbundanceGraph$Change <- as.factor(df.DifferentialAbundanceGraph$Change)
df.DifferentialAbundanceGraph$Change <- factor(df.DifferentialAbundanceGraph$Change, 
                                               levels = c("Positive",
                                                          "Negative"))


df.DifferentialAbundanceGraph$Fire <- as.factor(df.DifferentialAbundanceGraph$Fire)
levels(df.DifferentialAbundanceGraph$Fire)[levels(df.DifferentialAbundanceGraph$Fire)=="Unburned"] <- "Unburned prairie"
levels(df.DifferentialAbundanceGraph$Fire)[levels(df.DifferentialAbundanceGraph$Fire)=="Burned"] <- "Burned prairie"
df.DifferentialAbundanceGraph$Fire <- factor(df.DifferentialAbundanceGraph$Fire, 
                                       levels = c("Burned prairie",
                                                  "Unburned prairie"))

df.DifferentialAbundanceGraph$Fert <- as.factor(df.DifferentialAbundanceGraph$Fert)
df.DifferentialAbundanceGraph$Fert <- factor(df.DifferentialAbundanceGraph$Fert,
                                             levels = c("Recovering",
                                                        "Fertilized"))

MCCRecoverycolor <- c("#d95f02", "#7570b3")
MCCRecoveryshape <- c(22, 24)

LFC_label = expression(paste("Log"["2"]," Fold Change"))

# Highlight lines from tiff() to dev.off() to recreate figure.
tiff(file="Microbial_DifferentialAbundance.tiff", width = 7, height = 4, 
     pointsize = 1/300, units = 'in', res = 300)
ggplot(df.DifferentialAbundanceGraph, aes(x = Fire, y = log2FoldChange, 
                                    color = Fert))+
  geom_boxplot(linewidth = 1, outlier.shape = NA, fill = "white")+
  geom_point(aes(fill = Fert, shape = Fert), 
             size = 2, position = position_jitterdodge(), alpha = 1/3)+
  scale_colour_manual(values = MCCRecoverycolor)+
  scale_fill_manual(values = MCCRecoverycolor)+
  scale_shape_manual(values = MCCRecoveryshape)+
  labs(y = LFC_label)+
  facet_wrap(~Change, scales = "free")+
  theme_bw()+
  theme(axis.text.y = element_text(colour="black", size=10),
        axis.text.x = element_text(colour="black", size=10),
        axis.title.x=element_blank(),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        legend.title = element_blank(),
        legend.position = "none")
dev.off()
# Figure S45 was further annotated in Microsoft Powerpoint

# Script is to extract a legend
minimal <- ggplot(df.DifferentialAbundanceGraph, aes(x = Fire, y = log2FoldChange, 
                                          color = Fert))+
  geom_point(aes(fill = Fert, shape = Fert), 
             size = 2, position = position_jitterdodge())+
  scale_colour_manual(values = MCCRecoverycolor)+
  scale_fill_manual(values = MCCRecoverycolor)+
  scale_shape_manual(values = MCCRecoveryshape)+
  labs(y = LFC_label)+
  facet_wrap(~Change, scales = "free")+
  theme_bw()+
  theme(axis.text.y = element_text(colour="black", size=10),
        axis.text.x = element_text(colour="black", size=10),
        axis.title.x=element_blank(),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        legend.title = element_blank(),
        legend.position = "none")


leg <- get_legend(minimal)

tiff("MicrobialDifferntialAbundance_Legend.tiff", width = 4, height = 4, pointsize = 1/300, units = 'in', res = 300)
as_ggplot(leg)
dev.off()


## END

