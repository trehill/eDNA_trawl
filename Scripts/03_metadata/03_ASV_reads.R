#Author: Tessa Rehill 
#determining number of initial reads, and reads after filtering for each sample
#creating an ASV by sample table 
#create an ASV by taxonomy table 

library(tidyverse)
library(here)
library(vegan)
library(usedist)
library(taxize)
library(janitor)
library(dplyr)

#Determining initial reads ####
#12se
reads_12se <- read.delim("Raw_data/eDNA/12s/12s_e/retained_reads.filterAndTrim_step.length_var.txt",
                    h=TRUE,
                    fill = TRUE) %>%
  `colnames<-`(c("sample", "reads_in", "reads_out", "percentage_retained")) %>%
  as.data.frame()

reads_12se$sample <- gsub("_L001_R1_001.fastq.gz", "", reads_12se$sample)

#numbers
mean(reads_12se$reads_in) #70844.27
mean(reads_12se$reads_out) #66632.72

#12su
reads_12su <- read.delim("Raw_data/eDNA/12s/12s_u/retained_reads.filterAndTrim_step.length_var.txt",
                         h=TRUE,
                         fill = TRUE) %>%
  `colnames<-`(c("sample", "reads_in", "reads_out", "percentage_retained")) %>%
  as.data.frame()

reads_12su$sample <- gsub("_L001_R1_001.fastq.gz", "", reads_12su$sample)

#numbers
mean(reads_12su$reads_in) #75406
mean(reads_12su$reads_out) #71488


#BOTH 
combined_df <- rbind(reads_12se, reads_12su)

mean(combined_df$reads_in)
mean(combined_df$reads_out)

write_csv(combined_df,
          here("Processed_data",
               "datasets",
               "reads_retained.csv"))


