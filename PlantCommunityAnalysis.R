# TITLE:        Belowground Plot Experiment - Plant community comp. analysis 
# AUTHOR:       Matthew Nieland
# 
# DATA INPUT:   BGPVC dataset from EDI, 
#               PPS01 dataset from EDI
#
# OUTPUT:       Plant_NMDS_YearColor_ALL.tiff (Figure 1)
#               Alpha_PlantRecovery.tiff (Figure 2a-c)
#               Alpha_PlantRecoveryLineGraph.tiff (Figure S2)
#               PlantTotalFunctionalTypeCover_Boxplot.tiff (Figure 2d-g)
#               BurnedPrairieSIMPERspecies.tiff (Figure S3)
#               UnburnedPrairieSIMPERspecies.tiff (Figure S3)
#
# PROJECT:      Belowground Plot Experiment (KNZ LTER) fertilization legacies
#
# DATE:         June 2026

# This script is used for the plant community composition analysis of the 
# project. It first imports the data from the Environmental Data Initiative (EDI),
# and then runs analysis, starting beta diversity and then alpha diversity (in 
# order of results in manuscript).

rm(list=ls()) # clear environment

# Set working directory
setwd("~/Desktop/Research/BGPE/MCCRecovery") # change as necessary



# The script to import data from EDI was copy-and-pasted.

# Package ID: knb-lter-knz.17.14 Cataloging System:https://pasta.edirepository.org.
# Data set title: BGPVC Plant species composition in the Belowground Plot Experiment at Konza Prairie.
# Data set creator:  John Blair -  
# Metadata Provider:    - Konza LTER 
# Contact:  Konza LTER -    - knzlter@ksu.edu
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-knz/17/14/410e032a0651ce990c8c497be62c68f7" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "Datacode",     
                 "RecYear",     
                 "RecMonth",     
                 "RecDay",     
                 "Repsite",     
                 "SpeciesCode",     
                 "Ab_genus",     
                 "Ab_species",     
                 "Plot",     
                 "CoverClass",     
                 "Pid",     
                 "Comments"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$Datacode)!="factor") dt1$Datacode<- as.factor(dt1$Datacode)
if (class(dt1$RecYear)=="factor") dt1$RecYear <-as.numeric(levels(dt1$RecYear))[as.integer(dt1$RecYear) ]               
if (class(dt1$RecYear)=="character") dt1$RecYear <-as.numeric(dt1$RecYear)
if (class(dt1$RecMonth)=="factor") dt1$RecMonth <-as.numeric(levels(dt1$RecMonth))[as.integer(dt1$RecMonth) ]               
if (class(dt1$RecMonth)=="character") dt1$RecMonth <-as.numeric(dt1$RecMonth)
if (class(dt1$RecDay)=="factor") dt1$RecDay <-as.numeric(levels(dt1$RecDay))[as.integer(dt1$RecDay) ]               
if (class(dt1$RecDay)=="character") dt1$RecDay <-as.numeric(dt1$RecDay)
if (class(dt1$Repsite)!="factor") dt1$Repsite<- as.factor(dt1$Repsite)
if (class(dt1$SpeciesCode)=="factor") dt1$SpeciesCode <-as.numeric(levels(dt1$SpeciesCode))[as.integer(dt1$SpeciesCode) ]               
if (class(dt1$SpeciesCode)=="character") dt1$SpeciesCode <-as.numeric(dt1$SpeciesCode)
if (class(dt1$Ab_genus)!="factor") dt1$Ab_genus<- as.factor(dt1$Ab_genus)
if (class(dt1$Ab_species)!="factor") dt1$Ab_species<- as.factor(dt1$Ab_species)
if (class(dt1$Plot)=="factor") dt1$Plot <-as.numeric(levels(dt1$Plot))[as.integer(dt1$Plot) ]               
if (class(dt1$Plot)=="character") dt1$Plot <-as.numeric(dt1$Plot)
if (class(dt1$CoverClass)!="factor") dt1$CoverClass<- as.factor(dt1$CoverClass)
if (class(dt1$Pid)=="factor") dt1$Pid <-as.numeric(levels(dt1$Pid))[as.integer(dt1$Pid) ]               
if (class(dt1$Pid)=="character") dt1$Pid <-as.numeric(dt1$Pid)
if (class(dt1$Comments)!="factor") dt1$Comments<- as.factor(dt1$Comments)

# Convert Missing Values to NA for non-dates

dt1$Repsite <- as.factor(ifelse((trimws(as.character(dt1$Repsite))==trimws("blank")),NA,as.character(dt1$Repsite)))
dt1$SpeciesCode <- ifelse((trimws(as.character(dt1$SpeciesCode))==trimws("blank")),NA,dt1$SpeciesCode)               
suppressWarnings(dt1$SpeciesCode <- ifelse(!is.na(as.numeric("blank")) & (trimws(as.character(dt1$SpeciesCode))==as.character(as.numeric("blank"))),NA,dt1$SpeciesCode))

sessionInfo()

# Here is the structure of the input data frame:
str(dt1)                            
attach(dt1)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(Datacode)
summary(RecYear)
summary(RecMonth)
summary(RecDay)
summary(Repsite)
summary(SpeciesCode)
summary(Ab_genus)
summary(Ab_species)
summary(Plot)
summary(CoverClass)
summary(Pid)
summary(Comments) 
# Get more details on character variables

summary(as.factor(dt1$Datacode)) 
summary(as.factor(dt1$Repsite)) 
summary(as.factor(dt1$Ab_genus)) 
summary(as.factor(dt1$Ab_species)) 
summary(as.factor(dt1$CoverClass)) 
summary(as.factor(dt1$Comments))
detach(dt1)

stat(dt1$CoverClass)


# Save data as .csv file, so the above script doesn't have to be repeated.
write.csv(dt1, "BGP_PlantComp.csv")


# Next, download plant functional group dataset from KNZ LTER

# Package ID: knb-lter-knz.134.4 Cataloging System:https://pasta.edirepository.org.
# Data set title: PPS01 Konza prairie plant species list.
# Data set creator:  Jesse Nippert -  
# Data set creator:  John Blair -  
# Data set creator:  Jeffrey Taylor -  
# Metadata Provider:    - Konza LTER 
# Contact:    -  Konza LTER  - 
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-knz/134/4/4abeadee64638e46bd5088f2fa7d832e" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt.functional <-read.csv(infile1,header=F 
                         ,skip=1
                         ,sep=","  
                         ,quot='"' 
                         , col.names=c(
                           "updatedyear",     
                           "code",     
                           "gen",     
                           "spec",     
                           "genus",     
                           "species",     
                           "family",     
                           "lifespan",     
                           "growthform",     
                           "origin",     
                           "photo",     
                           "Comments"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt.functional$updatedyear)=="factor") dt.functional$updatedyear <-as.numeric(levels(dt.functional$updatedyear))[as.integer(dt.functional$updatedyear) ]               
if (class(dt.functional$updatedyear)=="character") dt.functional$updatedyear <-as.numeric(dt.functional$updatedyear)
if (class(dt.functional$code)=="factor") dt.functional$code <-as.numeric(levels(dt.functional$code))[as.integer(dt.functional$code) ]               
if (class(dt.functional$code)=="character") dt.functional$code <-as.numeric(dt.functional$code)
if (class(dt.functional$gen)!="factor") dt.functional$gen<- as.factor(dt.functional$gen)
if (class(dt.functional$spec)!="factor") dt.functional$spec<- as.factor(dt.functional$spec)
if (class(dt.functional$genus)!="factor") dt.functional$genus<- as.factor(dt.functional$genus)
if (class(dt.functional$species)!="factor") dt.functional$species<- as.factor(dt.functional$species)
if (class(dt.functional$family)!="factor") dt.functional$family<- as.factor(dt.functional$family)
if (class(dt.functional$lifespan)!="factor") dt.functional$lifespan<- as.factor(dt.functional$lifespan)
if (class(dt.functional$growthform)!="factor") dt.functional$growthform<- as.factor(dt.functional$growthform)
if (class(dt.functional$origin)!="factor") dt.functional$origin<- as.factor(dt.functional$origin)
if (class(dt.functional$photo)!="factor") dt.functional$photo<- as.factor(dt.functional$photo)
if (class(dt.functional$Comments)!="factor") dt.functional$Comments<- as.factor(dt.functional$Comments)

# Convert Missing Values to NA for non-dates

dt.functional$lifespan <- as.factor(ifelse((trimws(as.character(dt.functional$lifespan))==trimws(".")),NA,as.character(dt.functional$lifespan)))
dt.functional$growthform <- as.factor(ifelse((trimws(as.character(dt.functional$growthform))==trimws(".")),NA,as.character(dt.functional$growthform)))
dt.functional$origin <- as.factor(ifelse((trimws(as.character(dt.functional$origin))==trimws(".")),NA,as.character(dt.functional$origin)))
dt.functional$photo <- as.factor(ifelse((trimws(as.character(dt.functional$photo))==trimws(".")),NA,as.character(dt.functional$photo)))


# Here is the structure of the input data frame:
str(dt.functional)                            
attach(dt.functional)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(updatedyear)
summary(code)
summary(gen)
summary(spec)
summary(genus)
summary(species)
summary(family)
summary(lifespan)
summary(growthform)
summary(origin)
summary(photo)
summary(Comments) 
# Get more details on character variables

summary(as.factor(dt.functional$gen)) 
summary(as.factor(dt.functional$spec)) 
summary(as.factor(dt.functional$genus)) 
summary(as.factor(dt.functional$species)) 
summary(as.factor(dt.functional$family)) 
summary(as.factor(dt.functional$lifespan)) 
summary(as.factor(dt.functional$growthform)) 
summary(as.factor(dt.functional$origin)) 
summary(as.factor(dt.functional$photo)) 
summary(as.factor(dt.functional$Comments))
detach(dt.functional)     



# Install R packages if not already
if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse") # this an aggregate of packages within the 
# tidyverse language, which will be used throughout the analysis

if (!require("vegan", quietly = TRUE))
  install.packages("vegan") # package for community ecology

if (!require("lme4", quietly = TRUE))
  install.packages("lme4") # package for linear mixed effect models

if (!require("lmerTest", quietly = TRUE))
  install.packages("lmerTest") # package similar to lme4, but include p-values

if (!require("emmeans", quietly = TRUE))
  install.packages("emmeans") # package for post-hoc tests

if (!require("cowplot", quietly = TRUE))
  install.packages("cowplot") # package for figures

if (!require("ggpubr", quietly = TRUE))
  install.packages("ggpubr") # package for figures

if (!require("multcomp", quietly = TRUE))
  install.packages("multcomp") # used for pairwise comparisons

if (!require("plotrix", quietly = TRUE))
  install.packages("plotrix") # used to calculate std error for graphs



# Load packages
library(tidyverse) # v2.0.0
library(vegan) # v2.7-1
library(lmerTest) # v3.1-3
library(emmeans) # v1.11.2
library(cowplot) # v1.2.0
library(ggpubr) # v0.6.1



dt.functional %>% 
  filter(code == 201)

# Line of code below creates an object that consolidates plant genus and species,
# and identifies plant species' growth form, lifespan, and photosynthetic pathway
functionaltype.list <- dt.functional %>% 
  unite(GenusSpecies, c(gen,spec), sep="_") %>% 
  mutate(PlantFunctionalGroup = case_when(growthform == "g" & lifespan == "p" & photo == "c4" ~ "C4perennialgrass",
                                          growthform == "g" & lifespan == "p" & photo == "c3" ~ "C3perennialgrass",
                                          growthform == "g" & lifespan == "a" ~ "Annualgrass",
                                          growthform == "f" & lifespan == "a" ~ "Annualforb",
                                          family == "fabaceae" ~ "Legume",
                                          growthform == "f" & lifespan == "p" & family != "fabaceae" ~ "NonNfixingperennialforb",
                                          growthform == "w" ~ "Woody")) %>% 
  dplyr::select(GenusSpecies, PlantFunctionalGroup)



# The next series of script calculates the abundance of plant species by
# using the midpoint of cover based on a modified Daubenmire scale (CoverClass),
# and then classifying species based on life history traits


# First, list plots of interest and their field treatments 
# i.e., burned vs. unburned; control vs. recovering


PlotsToKeep <- c(1, 2, 13, 16, # The field experiment has other treatments that 
                 18, 20, 25, 28, # are not to be included in analysis.
                 38, 40, 41, 42,
                 54, 56, 57, 60)

# To see the field experiment's layout, please refer to the following link: 
# https://github.com/manieland/CommunityRecovery_BGP_KNZ/tree/main [Accessed June 8, 2026]

BurnedPlots <- c(13, 16, 25, 28,
                 41, 42, 54, 56)

UnburnedPlots <- c(1, 2, 18, 20,
                   38, 40, 57, 60)

ControlPlots <- c(2, 13, 20, 25,
                  38, 42, 56, 57)

RecoveringPlots <- c(1, 16, 18, 28,
                     40, 41, 54, 60)

##### Non-metric multidimensional scaling (NMDS) for all years #####

# This script creates a new object from the species composition dataset, with
# species abundance averaged across field treatments (see Collins et al. 2021)
# doi: 10.1111/ele.13676

# This is doozy of a script, but here's what is going on:
# 1) Removes data we don't want (i.e., exclude 2022 data) and cleans plant names
# 2) Calculates midpoint of cover in Daubenmire scale
# 3) Spreads data to "wide" format, and adds 0 to cover when species weren't observed
# 4) Pivots back to "long" format to average species cover by field treatment for each year
# 5) Returns data back to "wide" format to matrix form necessary for analysis

PlantSpeciesCover_allYears <- dt1 %>% 
  # Step 1
  filter(!RecYear == 2022) %>% # dataset from 1989 to 2022
  filter(Plot %in% PlotsToKeep) %>% 
  dplyr::select(!(c(Datacode, RecMonth, RecDay, Pid, Comments, SpeciesCode))) %>% 
  unite(GenusSpecies, c(Ab_genus,Ab_species), sep="_") %>% 
  
  # Step 2
  mutate(PercentCover = case_when(
    CoverClass == "1" ~ 0.5,
    CoverClass == "2" ~ 3.0,
    CoverClass == "3" ~ 15.0,
    CoverClass == "4" ~ 37.5,
    CoverClass == "5" ~ 62.5,
    CoverClass == "6" ~ 85.0,
    CoverClass == "7" ~ 97.5)) %>% 
  filter(!GenusSpecies == "_") %>% 
  dplyr::select(!(CoverClass)) %>% 
  distinct() %>% # silphi_integ (Silphium integrifolium) is duplicated in 1989 collections
  
  # Step 3
  pivot_wider(names_from = GenusSpecies, 
              values_from = PercentCover) %>%
  replace(is.na(.), 0) %>% 
  
  # Step 4
  pivot_longer(
    cols = !c(RecYear, Repsite, Plot),
    names_to = "GenusSpecies", 
    values_to = "PercentCover"
  ) %>%  # return back to "long" format
  mutate(Fire = case_when( # Add fire and enrichment treatments to object
    Plot %in% BurnedPlots ~ "Burned",
    Plot %in% UnburnedPlots ~ "Unburned"),
    Fert = case_when(
      Plot %in% ControlPlots ~ "Control",
      Plot %in% RecoveringPlots ~ "Recovering")) %>% 
  group_by(RecYear, Fire, Fert, GenusSpecies) %>% 
  summarize(PercentCover = mean(PercentCover, na.rm = T)) %>% 
  
  # Step 5
  pivot_wider(names_from = GenusSpecies, 
              values_from = PercentCover) %>%
  unite(SampleID,c(RecYear,Fire, Fert), sep="_", remove = FALSE) %>% 
  tibble::column_to_rownames("SampleID")

# Next, pull out environmental data
comm.env <- PlantSpeciesCover_allYears %>%
  dplyr::select(RecYear,Fire,Fert)

# And species cover data
comm <- PlantSpeciesCover_allYears %>%
  dplyr::select(-c(RecYear,Fire,Fert))

# Next three lines makes sure that environmental data are factors
comm.env$Fire <- as.factor(comm.env$Fire)
comm.env$Fert <- as.factor(comm.env$Fert)
comm.env$RecYear <- as.factor(comm.env$RecYear)

# Now run the nmds using "metaMDS" function in vegan package
set.seed(20230413)
plant.NMDS <- metaMDS(comm,
                      distance = "bray",
                      k = 2,
                      maxit = 999,
                      trymax = 500,
                      wascores = TRUE)

plot(plant.NMDS) # base R graph to quickly check plot


###### NMDS publication figure (Figure 1) #####
nmds_points <- plant.NMDS$points %>%
  data.frame() %>% 
  rownames_to_column(var = "ID") %>% # Make row names be a column
  as_tibble()


comm.envMeta <- comm.env %>% 
  rownames_to_column(var = "ID") %>% # Make row names be a column
  as_tibble() # convert to a tibble

# Join NMDS points with meta data
nmds_plant <- left_join(comm.envMeta,nmds_points, by = "ID")

# Make sure factors are coded as factors
nmds_plant$Fire <- as.factor(nmds_plant$Fire)
nmds_plant$Fert <- as.factor(nmds_plant$Fert)
nmds_plant$RecYear <- as.factor(nmds_plant$RecYear)

# Because N enrichment was ceased, this change will be reflected in the NMDS figure
nmds_plant <- nmds_plant %>% 
  mutate(Fert = case_when(
    Fert == "Control" ~ "Control",
    Fert == "Recovering" & RecYear %in% c("1989", "1994", "1999",
                                          "2005", "2010", "2015",
                                          "2016") ~ "Fertilized", # Years with fertilization
    Fert == "Recovering" & RecYear %in% c("2017", "2018", "2019",
                                          "2020", "2021") ~ "Recovering" # Post-fertilization
  ),
  TrtCombo = case_when(
    Fire == "Burned" & Fert == "Control" ~ "BC",
    Fire == "Burned" & Fert == "Recovering" ~ "BR",
    Fire == "Burned" & Fert == "Fertilized" ~ "BN",
    Fire == "Unburned" & Fert == "Control" ~ "UBC",
    Fire == "Unburned" & Fert == "Recovering" ~ "UBR",
    Fire == "Unburned" & Fert == "Fertilized" ~ "UBN"
  ))

nmds_plant$TrtCombo <- factor(nmds_plant$TrtCombo, levels = c("BC","BR", "BN", 
                                                          "UBC", "UBR", "UBN"))

# Script from L431 to L477 is cleaning the data and preparing for the NMDS plot
nmds_plant <- nmds_plant %>% 
  unite("Grouping",c("Fert","RecYear"), 
        sep="", remove = FALSE)

nmds_plant$Grouping <- as.factor(nmds_plant$Grouping)


# Changes "Fire" levels to add "prairie".
levels(nmds_plant$Fire)[levels(nmds_plant$Fire)=="Unburned"] <- "Unburned prairie"
levels(nmds_plant$Fire)[levels(nmds_plant$Fire)=="Burned"] <- "Burned prairie"
nmds_plant$Fire <- factor(nmds_plant$Fire , levels = c("Burned prairie", 
                                                       "Unburned prairie"))


nmds_plant$Fert <- factor(nmds_plant$Fert , levels = c("Control",
                                                   "Recovering",
                                                   "Fertilized"))

# Making sure points connect on the order of collection (based on year)
nmds_plant <- nmds_plant %>% 
  mutate(Order = case_when(
    RecYear == "1989" ~ "1",
    RecYear == "1994" ~ "2",
    RecYear == "1999" ~ "3",
    RecYear == "2005" ~ "4",
    RecYear == "2010" ~ "5",
    RecYear == "2015" ~ "6",
    RecYear == "2016" ~ "7",
    RecYear == "2017" ~ "8",
    RecYear == "2018" ~ "9",
    RecYear == "2019" ~ "10",
    RecYear == "2020" ~ "11",
    RecYear == "2021" ~ "12"))

nmds_plant$Order <- as.numeric(nmds_plant$Order)


# Lighter shades with earlier data, with saturated color representing more 
# recent data
GroupingColor = c("#DBECE5", "#DBECE5", "#DBECE5", "#B7D7CA","#B7D7CA", 
                  "#93C4AF","#93C4AF", "#93C4AF", "#6FAF94", "#6FAF94", 
                  "#4B9C79", "#4B9C79", "#E3E2EF", "#E3E2EF","#E3E2EF", 
                  "#C8C6DF","#ABA9CE", "#ABA9CE", "#908DBE", "#F4E0D4",
                  "#EAC2A9","#DFA37E","#D58553","#CA6627")


Fertshape <- c(21, 22, 24)

# NMDS figure made with ggplot, saved as a .tiff file
# To recreate figure with exact dimension, highlight all of L489 to L507
# and then run script.

tiff(file="Plant_NMDS_YearColor_ALL.tiff", width = 4, height = 4, 
     pointsize = 1/300, units = 'in', res = 300)
ggplot(data = nmds_plant[order(nmds_plant$Order),], 
       aes(x = MDS1, y = MDS2, group = TrtCombo, 
           color = Fert, fill = Grouping, shape = Fert))+
  geom_path(aes(colour = Fert), linewidth = 0.75)+
  geom_point(aes(fill = Grouping), size = 4, color = "black")+
  scale_shape_manual(values = Fertshape)+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_manual(values=GroupingColor)+
  labs(x = "NMDS 1", y="NMDS 2")+
  theme_bw()+
  theme(axis.text.y = element_text(colour="black", size=12),
        axis.text.x = element_text(colour="black", size=12),
        strip.text.x = element_text(size = 12, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", 
                                        size = 1),
        legend.title = element_blank(),
        legend.position = "none") # To check if colors match treatments, remove
                                  # legend.position = "none"
                                  # Note: check line colors  for fertilization trt.
dev.off()
# Figure 1 was further annotated in Microsoft Powerpoint

# This figure includes legend that shows fire and fertilization treatments
ggplot(data = nmds_plant[order(nmds_plant$Order),], 
       aes(x = MDS1, y = MDS2, group = TrtCombo, 
           color = Fert, fill = Grouping, shape = Fire))+
  geom_point(aes(fill = Grouping), size = 4)+
  labs(x = "NMDS 1", y="NMDS 2")+
  theme_bw()+
  theme(axis.text.y = element_text(colour="black", size=12),
        axis.text.x = element_text(colour="black", size=12),
        strip.text.x = element_text(size = 12, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", 
                                        size = 1),
        legend.title = element_blank())



##### Permutational analysis of variation (PERMANOVA) for 2017-2021 #####
# Only 2017-2021 were analyzed through PERMANOVA, because we were interested if
# prescribed burns supported species composition recovery from long-term N enrichment.

PlantSpeciesCover_2017_2021 <- dt1 %>% 
  filter(between(RecYear, 2017, 2021), # Keeps data only between 2017-2021
         Plot %in% PlotsToKeep) %>% 
  dplyr::select(!(c(Datacode, RecMonth, RecDay, Pid, Comments, SpeciesCode))) %>% 
  unite(GenusSpecies, c(Ab_genus,Ab_species), sep="_") %>% 
  
  # Step 2
  mutate(PercentCover = case_when(
    CoverClass == "1" ~ 0.5,
    CoverClass == "2" ~ 3.0,
    CoverClass == "3" ~ 15.0,
    CoverClass == "4" ~ 37.5,
    CoverClass == "5" ~ 62.5,
    CoverClass == "6" ~ 85.0,
    CoverClass == "7" ~ 97.5)) %>% 
  filter(!GenusSpecies == "_") %>% 
  dplyr::select(!(CoverClass)) %>% 
  distinct() %>% # silphi_integ (Silphium integrifolium) is duplicated in 1989 collections
  
  # Step 3
  pivot_wider(names_from = GenusSpecies, 
              values_from = PercentCover) %>%
  replace(is.na(.), 0) %>% 
  
  # Step 4
  pivot_longer(
    cols = !c(RecYear, Repsite, Plot),
    names_to = "GenusSpecies", 
    values_to = "PercentCover"
  ) %>%  # return back to "long" format
  mutate(Fire = case_when( # Add fire and enrichment treatments to object
    Plot %in% BurnedPlots ~ "Burned",
    Plot %in% UnburnedPlots ~ "Unburned"),
    Fert = case_when(
      Plot %in% ControlPlots ~ "Control",
      Plot %in% RecoveringPlots ~ "Recovering")) %>% 
  group_by(RecYear, Fire, Fert, GenusSpecies) %>% 
  summarize(PercentCover = mean(PercentCover, na.rm = T)) %>% 
  
  # Step 5
  pivot_wider(names_from = GenusSpecies, 
              values_from = PercentCover) %>%
  unite(SampleID,c(RecYear,Fire, Fert), sep="_", remove = FALSE) %>% 
  tibble::column_to_rownames("SampleID")


# Next, pull out environmental data
comm.env_2017_2021 <- PlantSpeciesCover_2017_2021 %>%
  dplyr::select(RecYear,Fire,Fert)

# And species cover data
comm_2017_2021 <- PlantSpeciesCover_2017_2021 %>%
  dplyr::select(-c(RecYear,Fire,Fert))

# Next three lines makes sure that environmental data are factors
comm.env_2017_2021$Fire <- as.factor(comm.env_2017_2021$Fire)
comm.env_2017_2021$Fert <- as.factor(comm.env_2017_2021$Fert)
comm.env_2017_2021$RecYear <- as.factor(comm.env_2017_2021$RecYear)


# Compute Bray-Curtis dissimilarity matrix
bray_distance.plant <- vegdist(comm_2017_2021, method = "bray")

# Because of the permutational nature of analysis, it is necessary to set a seed
# for reproducibility 
set.seed(20230413)
Plant_adonis_v1 = adonis2(bray_distance.plant ~ Fire*Fert*RecYear,
                       by = "terms",
                       comm.env_2017_2021, permutations = 999)

Plant_adonis_v1 # Perfect fit; all variance explained

# Exclude three way interaction
set.seed(20230413)
Plant_adonis_v2 = adonis2(bray_distance.plant ~ Fire*Fert + RecYear*Fire + RecYear*Fert,
                       by = "terms",
                       comm.env_2017_2021, permutations = 999)

Plant_adonis_v2




##### Alpha diversity for plant community analysis #####

# Recreate plant species matrix, but with subplots kept separate

PlantSpeciesMatrix_2017_2021 <- dt1 %>% 
  filter(between(RecYear, 2017, 2021), # Keeps data only between 2017-2021
         Plot %in% PlotsToKeep) %>% 
  dplyr::select(!(c(Datacode, RecMonth, RecDay, Pid, Comments, SpeciesCode))) %>% 
  unite(GenusSpecies, c(Ab_genus,Ab_species), sep="_") %>% 
  
  # Step 2
  mutate(PercentCover = case_when(
    CoverClass == "1" ~ 0.5,
    CoverClass == "2" ~ 3.0,
    CoverClass == "3" ~ 15.0,
    CoverClass == "4" ~ 37.5,
    CoverClass == "5" ~ 62.5,
    CoverClass == "6" ~ 85.0,
    CoverClass == "7" ~ 97.5)) %>% 
  filter(!GenusSpecies == "_") %>% 
  dplyr::select(!(CoverClass)) %>% 
  distinct() %>% # silphi_integ (Silphium integrifolium) is duplicated in 1989 collections
  
  # Step 3
  pivot_wider(names_from = GenusSpecies, 
              values_from = PercentCover) %>%
  replace(is.na(.), 0) %>% 
  
  # Step 4
  pivot_longer(
    cols = !c(RecYear, Repsite, Plot),
    names_to = "GenusSpecies", 
    values_to = "PercentCover"
  ) %>%
  dplyr::select(RecYear, Repsite, GenusSpecies, Plot, PercentCover) %>% 
  
  # Step 5 (Subplots are now separate)
  group_by(RecYear, GenusSpecies, Plot) %>% # subplot coded as "plot"
  summarize(PercentCover = mean(PercentCover)) %>% 
  pivot_wider(names_from = GenusSpecies, 
              values_from = PercentCover) %>% 
  unite(SampleID, c(RecYear, Plot), sep="_", remove = TRUE) %>% 
  tibble::column_to_rownames("SampleID")



# Create meta data to align subplot (Plot in this script) with field treatment
# Easier to do this way in my opinion
PlantSpeciesMatrix_meta <- PlotsToKeep %>% 
  as_tibble() %>% 
  rename(Plot = value) %>% 
  mutate(Fire = case_when(
    Plot %in% BurnedPlots ~ "Burned",
    Plot %in% UnburnedPlots ~ "Unburned"
  ),
  Fert = case_when(
    Plot %in% ControlPlots ~ "Control",
    Plot %in% RecoveringPlots ~ "Recovering"
  ))

# Make sure PlantSpeciesMatrix_meta$Plot is a character to match format as the
# object PlantSpeciesMatrix_2017_2021.
PlantSpeciesMatrix_meta$Plot <- as.character(PlantSpeciesMatrix_meta$Plot)

shannon <- PlantSpeciesMatrix_2017_2021 %>% 
  mutate(Shannon = vegan::diversity(., index = "shannon")) %>% 
  data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  dplyr::select(ID, Shannon)

richness <- PlantSpeciesMatrix_2017_2021 %>% 
  mutate(Richness = vegan::specnumber(.)) %>% 
  data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  dplyr::select(ID, Richness)

Alpha_Plant <- shannon %>% 
  left_join(., richness, by = "ID") %>% 
  mutate(Evenness = Shannon/log(Richness)) %>% 
  separate(ID, into = c("RecYear", "Plot"), sep = "_") %>% 
  left_join(., PlantSpeciesMatrix_meta, by = "Plot") %>% 
  mutate(BlockNumber = case_when( # adding experimental plot (coded "block")
    Plot == c("1", "2") ~ "B1",
    Plot == c("13", "16") ~ "B2",
    Plot == c("18", "20") ~ "B3",
    Plot == c("25", "28") ~ "B4",
    Plot == c("38", "40") ~ "B5",
    Plot == c("41", "42") ~ "B6",
    Plot == c("54", "56") ~ "B7",
    Plot == c("57", "60") ~ "B8"
  ))


# Before running linear mixed models, make sure fixed effects are factors
Alpha_Plant$RecYear <- as.factor(Alpha_Plant$RecYear)
Alpha_Plant$Fire <- as.factor(Alpha_Plant$Fire)
Alpha_Plant$Fert <- as.factor(Alpha_Plant$Fert)
Alpha_Plant$BlockNumber <- as.factor(Alpha_Plant$BlockNumber)


# making sure factor levels are in order
Alpha_Plant$Fire <- factor(Alpha_Plant$Fire, levels = c("Burned", "Unburned"))
Alpha_Plant$Fert <- factor(Alpha_Plant$Fert, levels = c("Control", "Recovering"))
Alpha_Plant$RecYear <- factor(Alpha_Plant$RecYear, levels = c("2017", "2018",
                                                              "2019", "2020",
                                                              "2021"))

Alpha_Plant$BlockNumber <- factor(Alpha_Plant$BlockNumber, levels = c("B1",
                                                                      "B2",
                                                                      "B3",
                                                                      "B4",
                                                                      "B5",
                                                                      "B6",
                                                                      "B7",
                                                                      "B8"))

###### Richness #####
hist(Alpha_Plant$Richness) # Normally distributed

lm.Richness <- lmer(Richness ~ Fire*Fert*RecYear + (1|BlockNumber),
                    data = Alpha_Plant, contrasts=contr.sum, REML = TRUE,
                    na.action = na.exclude)

plot(lm.Richness) # appears homoscedastic
qqnorm(resid(lm.Richness)) 
qqline(resid(lm.Richness)) # QQ-plot looks OK
hist(resid(lm.Richness)) # residuals follow normal distribution

anova(lm.Richness)

emm <- emmeans(lm.Richness, pairwise~Fire:Fert, adjust = "tukey")
multcomp::cld(emm, Letters = letters, reversed = T)
# use multcomp:: because it masks the select function in dplyr

emm <- emmeans(lm.Richness, pairwise~RecYear, adjust = "tukey")
multcomp::cld(emm, Letters = letters, reversed = T)

###### Shannon #####
hist(Alpha_Plant$Shannon) # Less normally distributed than richness

lm.Shannon <- lmer(Shannon ~ Fire*Fert*RecYear + (1|BlockNumber),
                   data = Alpha_Plant, contrasts=contr.sum, REML = TRUE,
                   na.action = na.exclude)

plot(lm.Shannon) # appears homoscedastic
qqnorm(resid(lm.Shannon)) 
qqline(resid(lm.Shannon)) # QQ-plot does not look great
hist(resid(lm.Shannon)) # appears bimodal

# Try log-transformation
lm.Shannon <- lmer(log(Shannon+0.1) ~ Fire*Fert*RecYear + (1|BlockNumber),
                   data = Alpha_Plant, contrasts=contr.sum, REML = TRUE,
                   na.action = na.exclude)

plot(lm.Shannon) # appears homoscedastic
qqnorm(resid(lm.Shannon)) 
qqline(resid(lm.Shannon)) # QQ-plot looks better (residuals follow along line)
hist(resid(lm.Shannon)) # residuals appear normally distributed

anova(lm.Shannon)

emm <- emmeans(lm.Shannon, pairwise~Fire:Fert, adjust = "tukey")
multcomp::cld(emm, Letters = letters, reversed = T)


###### Evenness #####
hist(Alpha_Plant$Evenness) # Right-skewed

lm.Evenness <- lmer(Evenness ~ Fire*Fert*RecYear + (1|BlockNumber),
                    data = Alpha_Plant, contrasts=contr.sum, REML = TRUE,
                    na.action = na.exclude)

plot(lm.Evenness) # homoscedastic
qqnorm(resid(lm.Evenness)) 
qqline(resid(lm.Evenness)) # QQ-plot looks OK
hist(resid(lm.Evenness)) # residuals are normally distributed

anova(lm.Evenness)

emm <- emmeans(lm.Evenness, pairwise~Fire:Fert, adjust = "tukey")
multcomp::cld(emm, Letters = letters, reversed = T)


###### Plant alpha diversity publication figure (Figure 2a-c) #####

# The following script formats data to smoothly make ggplot figure

df.AlphaGraph <- Alpha_Plant %>% 
  gather(key = Metric, value = Value, Shannon, Richness, Evenness)

# Make sure that factors are factors
df.AlphaGraph$Treatment <- as.factor(df.AlphaGraph$Treatment)
df.AlphaGraph$Treatment <- factor(df.AlphaGraph$Treatment, 
                                  levels = c("BC", "BR", "UBC", "UBR"))

df.AlphaGraph$Metric <- as.factor(df.AlphaGraph$Metric)
df.AlphaGraph$Metric <- factor(df.AlphaGraph$Metric,
                               levels = c("Richness", "Shannon", "Evenness"))

df.AlphaGraph$Fire <- as.factor(df.AlphaGraph$Fire)
levels(df.AlphaGraph$Fire)[levels(df.AlphaGraph$Fire)=="Unburned"] <- "Unburned prairie"
levels(df.AlphaGraph$Fire)[levels(df.AlphaGraph$Fire)=="Burned"] <- "Burned prairie"

df.AlphaGraph$Fire <- factor(df.AlphaGraph$Fire, 
                             levels = c("Burned prairie", "Unburned prairie"))

df.AlphaGraph$Fert <- as.factor(df.AlphaGraph$Fert)
df.AlphaGraph$Fert <- factor(df.AlphaGraph$Fert, 
                             levels = c("Control", "Recovering"))


Plantcolor <- c("#1b9e77", "#d95f02")
Plantshape <- c(16, 15)

# Figure made with ggplot, saved as a .tiff file

tiff(file="PlantTotalFunctionalTypeCover_Boxplot.tiff",width = 7, height = 2.5, 
     pointsize = 1/300, units = 'in', res = 300)
ggplot(df.AlphaGraph, aes(x = Fire, y = Value, color = Fert, shape = Fert))+
  geom_boxplot(linewidth = 1, outlier.shape = NA, fill = "white")+
  geom_point(aes(fill = Fert), size = 2, position = position_jitterdodge(), alpha = 1/3)+
  scale_colour_manual(values = Plantcolor)+
  scale_fill_manual(values = Plantcolor)+
  scale_shape_manual(values = Plantshape)+
  theme_bw()+
  facet_wrap(~Metric, scales = "free")+
  theme_bw()+
  theme(axis.text.y = element_text(colour = "black", size=10),
        axis.text.x = element_text(colour = "black", size=10),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        legend.title = element_blank(),
        legend.position = "none")
dev.off()
# Figure 2a-c was further annotated in Microsoft Powerpoint

###### Plant alpha diversity line graph (Supplemental Figure 2) #####
df.AlphaGraph

df.AlphaGraph.v2 <- Alpha_Plant %>% # New object for line graph
  gather(key = Metric, value = Value, Richness, Shannon, Evenness) %>% 
  group_by(Metric, Fire, Fert, RecYear) %>%
  summarize(mean.metric = mean(Value, na.rm = T), 
            se.metric = plotrix::std.error(Value, na.rm = T)) %>% 
  mutate(Treatment = case_when(
    Fire == "Burned" & Fert == "Control" ~ "BC",
    Fire == "Burned" & Fert == "Recovering" ~ "BR",
    Fire == "Unburned" & Fert == "Control" ~ "UBC",
    Fire == "Unburned" & Fert == "Recovering" ~ "UBR"
  ))

# Re-making sure data are formatted correctly.
df.AlphaGraph.v2$Treatment <- as.factor(df.AlphaGraph.v2$Treatment)
df.AlphaGraph.v2$Treatment <- factor(df.AlphaGraph.v2$Treatment, 
                                     levels = c("BC","BR", "UBC", "UBR"))

df.AlphaGraph.v2$Metric <- as.factor(df.AlphaGraph.v2$Metric)
df.AlphaGraph.v2$Metric <- factor(df.AlphaGraph.v2$Metric, 
                                  levels = c("Richness","Shannon", "Evenness"))

df.AlphaGraph.v2$Fire <- as.factor(df.AlphaGraph.v2$Fire)
levels(df.AlphaGraph.v2$Fire)[levels(df.AlphaGraph.v2$Fire)=="Unburned"] <- "Unburned prairie"
levels(df.AlphaGraph.v2$Fire)[levels(df.AlphaGraph.v2$Fire)=="Burned"] <- "Burned prairie"

df.AlphaGraph.v2$Fire <- factor(df.AlphaGraph.v2$Fire, 
                                levels = c("Burned prairie", "Unburned prairie"))

df.AlphaGraph.v2$Fert <- as.factor(df.AlphaGraph.v2$Fert)
df.AlphaGraph.v2$Fert <- factor(df.AlphaGraph.v2$Fert, 
                                levels = c("Control", "Recovering"))


PlantRecoverycolor <- c("#1b9e77", "#d95f02", "#1b9e77", "#d95f02")
PlantRecoveryshape <- c(21, 22, 21, 22)

tiff(file="Alpha_PlantRecoveryLineGraph.tiff",width = 6, height = 7.5, 
     pointsize = 1/300, units = 'in', res = 300)
ggplot(df.AlphaGraph.v2, aes(x = RecYear, y = mean.metric, group = Treatment))+
  geom_line(data = df.AlphaGraph.v2[!is.na(df.AlphaGraph.v2$mean.metric),], 
            aes(colour = Treatment), size = 0.75)+
  geom_errorbar(aes(ymin = mean.metric-se.metric, ymax = mean.metric+se.metric),
                colour = "black", width = 0.1)+
  geom_point(aes(fill = Treatment, shape = Treatment), width = 0.15, size = 4, 
             color = "black")+
  scale_fill_manual(values = PlantRecoverycolor)+
  scale_colour_manual(values = PlantRecoverycolor)+
  scale_shape_manual(values = PlantRecoveryshape)+
  theme_bw()+
  labs(x = "Year")+
  facet_grid(Metric ~ Fire, scales = "free")+
  #scale_y_continuous(trans='log10', labels = scales::comma)+
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
# Figure S2 was further annotated in Microsoft Powerpoint


##### Plant functional group analysis #####

# Create a new object that will be used for both plant functional group analysis
# and SIMPER analysis.


PlantSpeciesCover_2017_2021 <- dt1 %>% 
  filter(between(RecYear, 2017, 2021),
         Plot %in% PlotsToKeep) %>% 
  dplyr::select(!(c(Datacode, RecMonth, RecDay, Pid, Comments, SpeciesCode))) %>% 
  unite(GenusSpecies, c(Ab_genus,Ab_species), sep="_") %>% 
  mutate(PercentCover = case_when(
    CoverClass == "1" ~ 0.5, # Midpoint of cover in Daubenmire scale
    CoverClass == "2" ~ 3.0,
    CoverClass == "3" ~ 15.0,
    CoverClass == "4" ~ 37.5,
    CoverClass == "5" ~ 62.5,
    CoverClass == "6" ~ 85.0,
    CoverClass == "7" ~ 97.5)) %>% 
  
  filter(!GenusSpecies == "_") %>% 
  dplyr::select(!(CoverClass)) %>% 
  distinct() %>% # silphi_integ (Silphium integrifolium) is duplicated in 1989 collections
  pivot_wider(names_from = GenusSpecies, 
              values_from = PercentCover) %>% # Spread dataset so 0 can be added
  # to species not found in a subplot within a year
  replace(is.na(.), 0) %>%
  
  pivot_longer(
    cols = !c(RecYear, Repsite, Plot),
    names_to = "GenusSpecies", 
    values_to = "PercentCover"
    ) %>% # return back to long format

  group_by(RecYear, GenusSpecies, Plot) %>% 
  summarize(PercentCover = mean(PercentCover)) %>% 
  left_join(., functionaltype.list, by = "GenusSpecies") %>% 
  mutate(BlockNumber = case_when( # adding experimental plot (coded "block")
    Plot == c("1", "2") ~ "B1",
    Plot == c("13", "16") ~ "B2",
    Plot == c("18", "20") ~ "B3",
    Plot == c("25", "28") ~ "B4",
    Plot == c("38", "40") ~ "B5",
    Plot == c("41", "42") ~ "B6",
    Plot == c("54", "56") ~ "B7",
    Plot == c("57", "60") ~ "B8"
  ))

PlantSpeciesCover_2017_2021$Plot <- as.character(PlantSpeciesCover_2017_2021$Plot)

# Object for SIMPER analysis figures
PlantSpeciesCover_2017_2021 <- PlantSpeciesCover_2017_2021 %>% 
  left_join(., PlantSpeciesMatrix_meta, by = "Plot")

# Object for plant functional group analysis
df.FunctionalCover <- PlantSpeciesCover_2017_2021 %>% 
  group_by(RecYear, Plot, BlockNumber, Fire, Fert, PlantFunctionalGroup) %>% 
  summarize(TotalCover = sum(PercentCover)) %>% 
  tidyr::drop_na() %>% 
  dplyr::filter(PlantFunctionalGroup %in% c("C4perennialgrass", 
                                            "Legume",
                                            "NonNfixingperennialforb",
                                            "Woody"))
# Make sure fixed effects are factors
df.FunctionalCover$RecYear <- as.factor(df.FunctionalCover$RecYear)
df.FunctionalCover$Fire <- as.factor(df.FunctionalCover$Fire)
df.FunctionalCover$Fert <- as.factor(df.FunctionalCover$Fert)

df.FunctionalCover$RecYear <- factor(df.FunctionalCover$RecYear,
                                     levels = c("2017", "2018", "2019",
                                                "2020", "2021")
                                     )

df.FunctionalCover$Fire <- factor(df.FunctionalCover$Fire,
                                  levels = c("Burned", "Unburned")
                                  )

df.FunctionalCover$Fert <- factor(df.FunctionalCover$Fert,
                                     levels = c("Control", "Recovering")
                                  )

###### C4 grass #####

# First filter data so only C4 grass is included
df.c4grass <- df.FunctionalCover %>% 
  filter(PlantFunctionalGroup == "C4perennialgrass")

hist(df.c4grass$TotalCover) # Bimodal, as expected

lm.c4grass <- lmer(TotalCover ~ Fire*Fert*RecYear + (1|BlockNumber),
                    data = df.c4grass, contrasts=contr.sum, REML = TRUE,
                    na.action = na.exclude)

plot(lm.c4grass) # appears somewhat homoscedastic, pattern with lower values
qqnorm(resid(lm.c4grass)) 
qqline(resid(lm.c4grass)) # QQ-plot looks OK
hist(resid(lm.c4grass)) # residuals nearly follow normal distribution

anova(lm.c4grass)


###### Legume #####

df.legume <- df.FunctionalCover %>% 
  filter(PlantFunctionalGroup == "Legume")

hist(df.legume$TotalCover) # Left-skewed

lm.legume <- lmer(TotalCover ~ Fire*Fert*RecYear + (1|BlockNumber),
                  data = df.legume, contrasts=contr.sum, REML = TRUE,
                  na.action = na.exclude)

plot(lm.legume) # appears heteroscedastic, apparent pattern with lower values
qqnorm(resid(lm.legume)) 
qqline(resid(lm.legume)) # QQ-plot looks OK
hist(resid(lm.legume)) # residuals nearly follow normal distribution

# Try log-transformation
lm.legume <- lmer(log(TotalCover+0.1)~Fire*Fert*RecYear + (1|BlockNumber) ,
                  data = df.legume, contrasts=contr.sum, REML = TRUE,
                  na.action = na.exclude)

plot(lm.legume) # appears less heteroscedastic
qqnorm(resid(lm.legume)) 
qqline(resid(lm.legume)) # QQ-plot looks fine
hist(resid(lm.legume)) # appears a little more normal distribution


# Go with log-transformation
anova(lm.legume)


emm <- emmeans(lm.legume, pairwise~Fire:Fert, adjust = "tukey")
multcomp::cld(emm, Letters = letters, reversed = T)

###### Non-N-fixing forb #####

df.nonNfixingperennialforb <- df.FunctionalCover %>% 
  filter(PlantFunctionalGroup == "NonNfixingperennialforb")

hist(df.nonNfixingperennialforb$TotalCover) # Left-skewed

lm.nonNfixingperennialforb <- lmer(TotalCover ~ Fire*Fert*RecYear + (1|BlockNumber),
                                   data = df.nonNfixingperennialforb, contrasts=contr.sum, REML = TRUE,
                                   na.action = na.exclude)

plot(lm.nonNfixingperennialforb) # appears heteroscedastic, spreads out to right
qqnorm(resid(lm.nonNfixingperennialforb)) 
qqline(resid(lm.nonNfixingperennialforb)) # QQ-plot looks fine
hist(resid(lm.nonNfixingperennialforb)) # normally distributed, with small right-tail

# Try log-transformation
lm.nonNfixingperennialforb <- lmer(log(TotalCover+0.1)~Fire*Fert*RecYear + (1|BlockNumber),
                                   data = df.nonNfixingperennialforb, contrasts=contr.sum, REML = TRUE,
                                   na.action = na.exclude)

plot(lm.nonNfixingperennialforb) # improvement; appears homoscedastic
qqnorm(resid(lm.nonNfixingperennialforb)) 
qqline(resid(lm.nonNfixingperennialforb)) # QQ-plot looks OK
hist(resid(lm.nonNfixingperennialforb)) # little less normal distribution, but still fine


# Go with log-transformation
anova(lm.nonNfixingperennialforb)


###### Woody #####

df.woody <- df.FunctionalCover %>% 
  filter(PlantFunctionalGroup == "Woody")

hist(df.woody$TotalCover) # Left-skewed

lm.woody <- lmer(TotalCover ~ Fire*Fert*RecYear + (1|BlockNumber),
                 data = df.woody, contrasts=contr.sum, REML = TRUE,
                 na.action = na.exclude)

plot(lm.woody) # appears homoscedastic
qqnorm(resid(lm.woody)) 
qqline(resid(lm.woody)) # QQ-plot does not look OK
hist(resid(lm.woody)) # normally distributed=

# Try log-transformation
lm.woody <- lmer(log(TotalCover+0.1)~Fire*Fert*RecYear + (1|BlockNumber),
                 data = df.woody, contrasts=contr.sum, REML = TRUE,
                 na.action = na.exclude)

plot(lm.woody) # still appears homoscedastic
qqnorm(resid(lm.woody)) 
qqline(resid(lm.woody)) # QQ-plot looks better
hist(resid(lm.woody)) # normal distribution


# Go with log-transformation
anova(lm.woody)


emm <- emmeans(lm.woody, pairwise~Fire:Fert, adjust = "tukey")
multcomp::cld(emm, Letters = letters, reversed = T)



###### Plant functional group publication figure (Figure 2d-g) #####
df.FunctionalGroupGraph <- df.FunctionalCover

# Make sure factors are coded as factors for ggplot
df.FunctionalGroupGraph$PlantFunctionalGroup <- as.factor(df.FunctionalGroupGraph$PlantFunctionalGroup)
levels(df.FunctionalGroupGraph$PlantFunctionalGroup)[levels(df.FunctionalGroupGraph$PlantFunctionalGroup)=="C4perennialgrass"] <- "C4 grass"
levels(df.FunctionalGroupGraph$PlantFunctionalGroup)[levels(df.FunctionalGroupGraph$PlantFunctionalGroup)=="NonNfixingperennialforb"] <- "Non-N-fixing perennial forb"

df.FunctionalGroupGraph$PlantFunctionalGroup <- factor(df.FunctionalGroupGraph$PlantFunctionalGroup, 
                                          levels = c("C4 grass","Non-N-fixing perennial forb", 
                                                     "Legume", "Woody"))

df.FunctionalGroupGraph$Fire <- as.factor(df.FunctionalGroupGraph$Fire)
levels(df.FunctionalGroupGraph$Fire)[levels(df.FunctionalGroupGraph$Fire)=="Unburned"] <- "Unburned prairie"
levels(df.FunctionalGroupGraph$Fire)[levels(df.FunctionalGroupGraph$Fire)=="Burned"] <- "Burned prairie"
df.FunctionalGroupGraph$Fire <- factor(df.FunctionalGroupGraph$Fire, 
                                       levels = c("Burned prairie", "Unburned prairie"))


df.FunctionalGroupGraph$Fert <- as.factor(df.FunctionalGroupGraph$Fert)
df.FunctionalGroupGraph$Fert <- factor(df.FunctionalGroupGraph$Fert, 
                                       levels = c("Control", "Recovering"))


Plantcolor <- c("#1b9e77", "#d95f02")
Plantshape <- c(16, 15)

tiff(file="PlantTotalFunctionalTypeCover_Boxplot.tiff", 
     width = 5.5, height = 5, pointsize = 1/300, units = 'in', res = 300)
ggplot(df.FunctionalGroupGraph, aes(x = Fire, y = TotalCover, color = Fert, shape = Fert))+
  geom_boxplot(linewidth = 1, outlier.shape = NA, fill = "white")+
  geom_point(aes(fill = Fert), size = 2, position = position_jitterdodge(), alpha = 1/3)+
  scale_colour_manual(values=Plantcolor)+
  scale_fill_manual(values=Plantcolor)+
  scale_shape_manual(values=Plantshape)+
  theme_bw()+
  labs(y = "Total cover (%)")+
  facet_wrap(~PlantFunctionalGroup, scales = "free")+
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

# Figure 2d-g was further annotated in Microsoft Powerpoint


# To get legend for figure, use following code
minimal <- ggplot(df.FunctionalGroupGraph, aes(x = Fire, y = TotalCover, color = Fert, shape = Fert))+
  geom_point(aes(fill = Fert), size = 2, position = position_jitterdodge())+
  scale_colour_manual(values=Plantcolor)+
  scale_fill_manual(values=Plantcolor)+
  scale_shape_manual(values=Plantshape)+
  theme_bw()+
  labs(y = "Total cover (%)")+
  facet_wrap(~PlantFunctionalGroup, scales = "free")+
  theme_bw()+
  theme(axis.text.y = element_text(colour="black", size=10),
        axis.text.x = element_text(colour="black", size=10),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        legend.title = element_blank())


leg <- get_legend(minimal)

tiff("PlantTotalFunctionalTypeCover_Boxplot_Legend.tiff", 
     width = 4, height = 4, pointsize = 1/300, units = 'in', res = 300)
as_ggplot(leg)
dev.off()


##### SIMPER analysis #####

# Use the same matrix for the plant alpha diversity analysis 
# i.e., PlantSpeciesMatrix_2017_2021

df.simper <- PlantSpeciesMatrix_2017_2021 %>% 
  data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  separate(ID, into = c("RecYear", "Plot"), sep = "_") %>% 
  left_join(., PlantSpeciesMatrix_meta, by = "Plot") %>% 
  mutate(Fire = case_when(
    Plot %in% BurnedPlots ~ "Burned",
    Plot %in% UnburnedPlots ~ "Unburned"
  ),
  Fert = case_when(
    Plot %in% ControlPlots ~ "Control",
    Plot %in% RecoveringPlots ~ "Recovering"
  ))


###### Burned prairie #####
# Pull out environmental data
BurnedPlant.env <- df.simper %>%
  filter(Fire == "Burned") %>% 
  dplyr::select(RecYear, Plot, Fire, Fert)

# Plant cover data
BurnedPlant.comm <- df.simper %>%
  filter(Fire == "Burned") %>% 
  dplyr::select(-c(RecYear, Plot, Fire, Fert))

simper.BurnedPrairie <- simper(BurnedPlant.comm, BurnedPlant.env$Fert)

summary(simper.BurnedPrairie)


###### Burned prairie SIMPER line graph (Supplemental Figure 3) #####

# Script helps make sure ggplot is plotted correctly

simper.species.burned <- PlantSpeciesCover_2017_2021 %>% 
  filter(Fire == "Burned",
         GenusSpecies %in% c("panicu_virga", "androp_gerar",
                             "schiza_scopa", "solida_misso",
                             "ambros_psilo", "teucri_canad",
                             "lesped_capit", "symphy_erico",
                             "asclep_syria", "sorgha_nutan")) %>% 
  group_by(Fire, Fert, RecYear, GenusSpecies) %>% 
  summarize(mean.cover = mean(PercentCover),
            se.cover = plotrix::std.error(PercentCover))

simper.species.burned$RecYear <- as.factor(simper.species.burned$RecYear)
simper.species.burned$Fire <- as.factor(simper.species.burned$Fire)
simper.species.burned$Fert <- as.factor(simper.species.burned$Fert)
simper.species.burned$GenusSpecies <- as.factor(simper.species.burned$GenusSpecies)

levels(simper.species.burned$Fire)[levels(simper.species.burned$Fire)=="Burned"] <- "Burned prairie"

simper.species.burned$Fert <- factor(simper.species.burned$Fert, 
                                     levels = c("Control", "Recovering"))

simper.species.burned$RecYear <- factor(simper.species.burned$RecYear, 
                                        levels = c("2017", "2018",
                                                   "2019", "2020", "2021"))

simper.species.burned$GenusSpecies <- factor(simper.species.burned$GenusSpecies, 
                                             levels = c("panicu_virga", "androp_gerar",
                                                        "schiza_scopa","solida_misso",
                                                        "ambros_psilo","teucri_canad",
                                                        "lesped_capit","symphy_erico",
                                                        "asclep_syria","sorgha_nutan"))

species.color = c("#AB231F", "#FCEF96",
                  "#F9DE62", "#DACD75",
                  "#61A88A", "#B4DCCD",
                  "#93BAAB", "#618778",
                  "#A0BA63", "#6D8734")

tiff(file="BurnedPrairieSIMPERspecies.tiff",width = 5, height = 3, pointsize = 1/300, units = 'in', res =
       300)
ggplot(simper.species.burned, aes(x = RecYear, y = mean.cover, group = GenusSpecies))+
  geom_errorbar(aes(ymin=mean.cover-se.cover,ymax=mean.cover+se.cover,
                    color = GenusSpecies),width=.1)+
  geom_line(data=simper.species.burned, aes(colour = GenusSpecies), size =0.65)+
  geom_point(aes(fill = GenusSpecies), shape = 21, width = 0.15, 
             size = 4, color = "black")+
  scale_fill_manual(values=species.color)+
  scale_colour_manual(values=species.color)+
  theme_bw()+
  labs(x = "Year",
       y = "Cover (%)")+
  facet_grid(Fire ~ Fert, scales = "free")+
  theme_bw()+
  theme(axis.text.y = element_text(colour="black", size=10),
        axis.text.x = element_text(colour="black", size=10),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        legend.title = element_blank(),
        legend.position = "none")
dev.off()
# Figure S3 was further annotated in Microsoft Powerpoint

# Legend for the burned prairie SIMPER figure
minimal <- ggplot(simper.species.burned, aes(x = RecYear, y = mean.cover, group = GenusSpecies))+
  geom_errorbar(aes(ymin=mean.cover-se.cover,ymax=mean.cover+se.cover,
                    color = GenusSpecies),width=.1)+
  geom_line(data=simper.species.burned, aes(colour = GenusSpecies), size =0.65)+
  geom_point(aes(fill = GenusSpecies), shape = 21, width = 0.15, size = 4, color = "black")+
  scale_fill_manual(values=species.color)+
  scale_colour_manual(values=species.color)+
  theme_bw()+
  labs(x = "Year",
       y = "Cover (%)")+
  facet_grid(Fire ~ Fert, scales = "free")+

  theme_bw()+
  theme(axis.text.y = element_text(colour="black", size=10),
        axis.text.x = element_text(colour="black", size=10),
        #axis.title.x=element_blank(),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        legend.title = element_blank())

leg <- get_legend(minimal)

tiff("BurnedPrairieSIMPERspecies_Legend.tiff", width = 4, height = 4, pointsize = 1/300, units = 'in', res = 300)
as_ggplot(leg)
dev.off()



###### Unburned prairie #####
# Pull out environmental data
UnburnedPlant.env <- df.simper %>%
  filter(Fire == "Unburned") %>% 
  dplyr::select(RecYear, Plot, Fire, Fert)

# Plant cover data
UnburnedPlant.comm <- df.simper %>%
  filter(Fire == "Unburned") %>% 
  dplyr::select(-c(RecYear, Plot, Fire, Fert))

simper.UnburnedPrairie <- simper(UnburnedPlant.comm, UnburnedPlant.env$Fert)

summary(simper.UnburnedPrairie)


###### Unburned prairie SIMPER line graph (Supplemental Figure 3) #####
simper.species.unburned <- PlantSpeciesCover_2017_2021 %>% 
  filter(Fire == "Unburned",
         GenusSpecies %in% c("junipe_virgi", "rubus_occid",
                             "gledit_triac", "androp_gerar")) %>% 
  group_by(Fire, Fert, RecYear, GenusSpecies) %>% 
  summarize(mean.cover = mean(PercentCover),
            se.cover = plotrix::std.error(PercentCover))


simper.species.unburned$RecYear <- as.factor(simper.species.unburned$RecYear)
simper.species.unburned$Fire <- as.factor(simper.species.unburned$Fire)
simper.species.unburned$Fert <- as.factor(simper.species.unburned$Fert)
simper.species.unburned$GenusSpecies <- as.factor(simper.species.unburned$GenusSpecies)

levels(simper.species.unburned$Fire)[levels(simper.species.unburned$Fire)=="Unburned"] <- "Unburned prairie"
simper.species.unburned$Fert <- factor(simper.species.unburned$Fert, 
                                       levels = c("Control", "Recovering"))
simper.species.unburned$RecYear <- factor(simper.species.unburned$RecYear, 
                                          levels = c("2017", "2018", "2019", 
                                                     "2020", "2021"))
simper.species.unburned$GenusSpecies <- factor(simper.species.unburned$GenusSpecies, 
                                               levels = c("junipe_virgi", 
                                                          "rubus_occid",
                                                          "gledit_triac",
                                                          "androp_gerar"))

species.color = c("#AB231F", "#F9DE62",
                  "#61A88A", "#A0BA63")

tiff(file="UnburnedPrairieSIMPERspecies.tiff",width = 5, height = 3, pointsize = 1/300, units = 'in', res =
       300)
ggplot(simper.species.unburned, aes(x = RecYear, y = mean.cover, group = GenusSpecies))+
  geom_errorbar(aes(ymin=mean.cover-se.cover,ymax=mean.cover+se.cover,
                    color = GenusSpecies),width=.1)+
  geom_line(data=simper.species.unburned, aes(colour = GenusSpecies), size =0.65)+
  geom_point(aes(fill = GenusSpecies), shape = 21, width = 0.15, size = 4, color = "black")+
  scale_fill_manual(values=species.color)+
  scale_colour_manual(values=species.color)+
  theme_bw()+
  labs(x = "Year",
       y = "Cover (%)")+
  facet_grid(Fire ~ Fert, scales = "free")+
  theme_bw()+
  theme(axis.text.y = element_text(colour="black", size=10),
        axis.text.x = element_text(colour="black", size=10),
        #axis.title.y=element_blank(),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        legend.title = element_blank(),
        legend.position = "none")
dev.off()


# Legend for unburned prairie SIMPER figure
minimal <- ggplot(simper.species.unburned, aes(x = RecYear, y = mean.cover, group = GenusSpecies))+
  geom_errorbar(aes(ymin=mean.cover-se.cover,ymax=mean.cover+se.cover,
                    color = GenusSpecies),width=.1)+
  geom_line(data=simper.species.unburned, aes(colour = GenusSpecies), size =0.65)+
  geom_point(aes(fill = GenusSpecies), shape = 21, width = 0.15, size = 4, color = "black")+
  scale_fill_manual(values=species.color)+
  scale_colour_manual(values=species.color)+
  theme_bw()+
  labs(x = "Year",
       y = "Cover (%)")+
  facet_grid(Fire ~ Fert, scales = "free")+
  #scale_y_continuous(trans='log10')+
  #scale_y_continuous(limits = c(0, 75), expand = c(0,0))+
  theme_bw()+
  theme(axis.text.y = element_text(colour="black", size=10),
        axis.text.x = element_text(colour="black", size=10),
        #axis.title.x=element_blank(),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        legend.title = element_blank())

# Figure S3 was further annotated in Microsoft Powerpoint

leg <- get_legend(minimal)

tiff("UnburnedPrairieSIMPERspecies_Legend.tiff", width = 4, height = 4, pointsize = 1/300, units = 'in', res = 300)
as_ggplot(leg)
dev.off()


## END
