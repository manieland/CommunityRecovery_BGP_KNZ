# TITLE:        Belowground Plot Experiment - Precipitation 
# AUTHOR:       Matthew Nieland
# 
# DATA INPUT:   APT01 dataset from EDI
#
# OUTPUT:       2017-2021CumulativePrecipitation.tif
#
# PROJECT:      Belowground Plot Experiment (KNZ LTER) fertilization legacies
#
# DATE:         June 2026

# This script is used to show the variable of precipitation during the study.

rm(list=ls()) # clear environment

# Set working directory
setwd("~/Desktop/Research/BGPE/MCCRecovery") # change as necessary


# The script to import data from EDI was copy-and-pasted.

# Package ID: knb-lter-knz.4.17 Cataloging System:https://pasta.edirepository.org.
# Data set title: APT01 Daily precipitation amounts measured at multiple sites across konza prairie.
# Data set creator:  Jesse Nippert -  
# Metadata Provider:    - Konza LTER 
# Contact:  Konza LTER -    - knzlter@ksu.edu
# Contact:  Jeffrey Taylor -    - jht@ksu.edu
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-knz/4/17/434849250e532ae2029287269e677770" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "DataCode",     
                 "RecType",     
                 "RecDate",     
                 "watershed",     
                 "ppt",     
                 "Comments"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$DataCode)!="factor") dt1$DataCode<- as.factor(dt1$DataCode)
if (class(dt1$RecType)=="factor") dt1$RecType <-as.numeric(levels(dt1$RecType))[as.integer(dt1$RecType) ]               
if (class(dt1$RecType)=="character") dt1$RecType <-as.numeric(dt1$RecType)                                   
# attempting to convert dt1$RecDate dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1RecDate<-as.Date(dt1$RecDate,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1RecDate) == length(tmp1RecDate[!is.na(tmp1RecDate)])){dt1$RecDate <- tmp1RecDate } else {print("Date conversion failed for dt1$RecDate. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1RecDate) 
if (class(dt1$watershed)!="factor") dt1$watershed<- as.factor(dt1$watershed)
if (class(dt1$ppt)=="factor") dt1$ppt <-as.numeric(levels(dt1$ppt))[as.integer(dt1$ppt) ]               
if (class(dt1$ppt)=="character") dt1$ppt <-as.numeric(dt1$ppt)
if (class(dt1$Comments)!="factor") dt1$Comments<- as.factor(dt1$Comments)

# Convert Missing Values to NA for non-dates

dt1$ppt <- ifelse((trimws(as.character(dt1$ppt))==trimws(".")),NA,dt1$ppt)               
suppressWarnings(dt1$ppt <- ifelse(!is.na(as.numeric(".")) & (trimws(as.character(dt1$ppt))==as.character(as.numeric("."))),NA,dt1$ppt))


# Here is the structure of the input data frame:
str(dt1)                            
attach(dt1)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(DataCode)
summary(RecType)
summary(RecDate)
summary(watershed)
summary(ppt)
summary(Comments) 
# Get more details on character variables

summary(as.factor(dt1$DataCode)) 
summary(as.factor(dt1$watershed)) 
summary(as.factor(dt1$Comments))
detach(dt1)  


# open packages into library
library(tidyverse) # v2.0.0


###### Precipitation figure (Figure S1) #####
# Follow code in order!


# how is watershed defined in the dataframe?
class(dt1$watershed) # factors

# get column names of the dataframe
colnames(dt1)

# create new dataframe to remove gauges from other watersheds
dt2 <- subset(dt1,watershed=='HQ',
              select=c(RecDate,ppt,watershed))
dim(dt2)

# define date precipitation was collected as a 'date' type
dt2$RecDate <- as.Date(dt2$RecDate,format="%m/%d/%Y")

# create a new dataframe to define years of interest to calculate average precipitation
dt3 <- with(dt2,dt2[RecDate >= "1986-01-01", ])
tail(dt3) # 2021 is the last year

# calculate mean and standard deviation
# if you don't have the lubridate package, install it
dt3 %>% 
  dplyr::mutate(Year = lubridate::year(RecDate)) %>% 
  dplyr::group_by(Year) %>% 
  summarize(YearPPT = sum(ppt, na.rm = TRUE)) %>% 
  summarize(Mean1986_2021 = mean(YearPPT),
            StDev1986_2021 = sd(YearPPT))

# in case you want to know specific years
MeanYearPPT <- dt3 %>% 
  dplyr::mutate(Year = lubridate::year(RecDate)) %>% 
  dplyr::group_by(Year) %>% 
  summarize(YearPPT = sum(ppt, na.rm = TRUE))

# this calculates z score for the entire year. be wary though that this may not reflect
# 'dry' or 'wet' growing season (see Broderick et al. 2022 Global Change Biology)
dt3 %>% 
  dplyr::mutate(Year = lubridate::year(RecDate)) %>% 
  dplyr::group_by(Year) %>% 
  summarize(YearPPT = sum(ppt, na.rm = TRUE)) %>% 
  filter(Year > 2016) %>% 
  mutate(z = (YearPPT-849)/222)

# the bottom code calculates cumulative precipitation. it's VERY messy, fyi
DOYdataframe <- as.data.frame(seq(as.Date("1986-01-01"), as.Date("2021-12-31"),by="day"))
colnames(DOYdataframe)
names(DOYdataframe)[names(DOYdataframe) == "seq(as.Date(\"1986-01-01\"), as.Date(\"2021-12-31\"), by = \"day\")"] <- "RecDate"

DOYdataframe$DOY <- lubridate::yday(DOYdataframe$RecDate)

CumulativePPT <- left_join(DOYdataframe, dt3, by="RecDate")

CumulativePPT[is.na(CumulativePPT)] <- 0 # Add '0' to dates with no recorded precipitation

CumulativePPT[is.na(CumulativePPT)] <- 'HQ' # Add 'HQ' 

# Calculates average cumulative precipitation
df1 <- CumulativePPT %>%
  dplyr::mutate(Year = lubridate::year(RecDate)) %>%
  group_by(Year) %>% 
  dplyr::mutate(cum_ppt = cumsum(ppt)) %>% 
  dplyr::group_by(DOY) %>% 
  summarize(MeanCumPPT = mean(cum_ppt),
            STDCumPPT = sd(cum_ppt))

head(df1)

# i like to show 1 standard deviation of the cumulative precipitation on the figure
# this is how you can do that
df2 <- df1 %>% 
  mutate(high = MeanCumPPT+STDCumPPT,
         low = MeanCumPPT-STDCumPPT) %>% 
  filter(DOY != 366) # remove the 366th day of the year (leap years)

# replace negative 1 STD values with zero
df2[df2 < 0] <- 0 

tail(df2)

# plot for only growing season (April 1 - September 30)
class(df2$DOY)
ggplot()+
  geom_ribbon(data=df2,aes(ymin = low, ymax = high, x=DOY), alpha=0.2)+
  geom_line(data=df2,stat="identity",aes(x=DOY,y=MeanCumPPT),color="black")+
  scale_x_continuous(limit = c(91,273))+
  scale_y_continuous(expand=c(0,0),limit=c(0,1000))+
  labs(y="Cumulative precipitation (mm)",
       x="Day of Year")+
  theme_classic()

# facet so each year is highlighted
head(df2)
df3 <- df2 %>% 
  add_column(Year = "Average")


df4_Years <- CumulativePPT %>%
  dplyr::mutate(Year = lubridate::year(RecDate)) %>%
  group_by(Year) %>% 
  dplyr::mutate(MeanCumPPT = cumsum(ppt),
                STDCumPPT = 0,
                high = 0,
                low = 0) %>% 
  filter(Year > 2016) %>%
  filter(DOY != 366)

df4_Years$Year <- as.character(df4_Years$Year)
combined_dataframe <- full_join(df3, df4_Years)

# Check the data, this will make an "ugly" ggplot figure
ggplot()+
  geom_ribbon(data=combined_dataframe,aes(ymin = low, ymax = high, x=DOY), alpha=0.2)+
  geom_line(data=combined_dataframe,stat="identity",aes(x=DOY,y=MeanCumPPT, color = Year))+
  scale_x_continuous(limit = c(91,273))+
  scale_y_continuous(expand=c(0,0),limit=c(0,1000))+
  labs(y="Cumulative precipitation (mm)",
       x="Day of Year")+
  facet_grid(. ~ Year)+
  theme_classic()


SampleDates <- data.frame(Year = c(2017, 2017,
                                   2018, 2018,
                                   2019, 2019,
                                   2020, 2020,
                                   2021, 2021),
                          DOY = c(160, 230,
                                   157, 225,
                                   156, 232,
                                   163, 225,
                                   161, 223))

SampleDates$Year <- as.character(SampleDates$Year)


tiff(file="2017-2021CumulativePrecipitation.tiff",width = 9, height = 3.5, pointsize = 1/300, units = 'in', res = 300)
ggplot()+
  geom_ribbon(data=df2,aes(ymin = low, ymax = high, x=DOY), alpha=0.2)+
  geom_line(data=df2,stat="identity",aes(x=DOY,y=MeanCumPPT), color = "black", linewidth = 0.75)+
  geom_line(data=df4_Years,stat="identity",aes(x=DOY,y=MeanCumPPT), color = "#e7298a", linewidth = 1)+
  scale_x_continuous(limit = c(91,273))+
  scale_y_continuous(expand=c(0,0),limit=c(0,1000))+
  labs(y="Cumulative precipitation (mm)",
       x="Day of Year")+
  facet_grid(. ~ Year)+
  geom_vline(data = SampleDates, aes(xintercept = DOY), linetype = "dashed")+
  theme_bw()+
  theme(text = element_text(family = "sans", color = "black"),
        axis.text.y = element_text(colour="black", size=10),
        axis.text.x = element_text(colour="black", size=10),
        #axis.title.x=element_blank(),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        legend.title = element_blank(),
        legend.position = "none")
dev.off()

# END

