#Comparing new and old analysis 

#Detections A ###
#read in files
eDNA_old <-  read.csv(here::here("Processed_data",
                                 "datasets",
                                 "2021",
                                 "detections_all_A.csv"), 
                      head=TRUE)

eDNA_new <-  read.csv(here::here("Processed_data",
                                 "datasets",
                                 "detections_all_A.csv"), 
                      head=TRUE)

eDNA_old <- eDNA_old %>%
  filter(gamma_detection_method %in% c('only eDNA', 'both eDNA/trawl'))

eDNA_new <- eDNA_new %>%
  filter(gamma_detection_method %in% c('only eDNA', 'both eDNA/trawl'))

## count number of unique names 
unique_names <- unique(eDNA_new$LCT)
num_unique_names <- length(unique_names)
print(num_unique_names)


unique_names <- unique(eDNA_old$LCT)
num_unique_names <- length(unique_names)
print(num_unique_names)

##
unique(eDNA_new$LCT)

unique(eDNA_old$LCT)

common_species <- intersect(unique(eDNA_new$LCT), unique(eDNA_old$LCT))
print(common_species)

species_only_in_new <- setdiff(unique(eDNA_new$LCT), unique(eDNA_old$LCT))
print(species_only_in_new)

species_only_in_old <- setdiff(unique(eDNA_old$LCT), unique(eDNA_new$LCT))
print(species_only_in_old)


#Detections B 

#read in files
eDNA_old <-  read.csv(here::here("Processed_data",
                                 "datasets",
                                 "2021",
                                 "detections_all_B.csv"), 
                      head=TRUE)

eDNA_new <-  read.csv(here::here("Processed_data",
                                 "datasets",
                                 "detections_all_B.csv"), 
                      head=TRUE)

eDNA_old <- eDNA_old %>%
  filter(gamma_detection_method %in% c('only eDNA', 'both eDNA/trawl'))

eDNA_new <- eDNA_new %>%
  filter(gamma_detection_method %in% c('only eDNA', 'both eDNA/trawl'))

##
unique(eDNA_new$LCT)

unique(eDNA_old$LCT)

common_species <- intersect(unique(eDNA_new$LCT), unique(eDNA_old$LCT))
print(common_species)

species_only_in_df <- setdiff(unique(eDNA_new$LCT), unique(eDNA_old$LCT))
print(species_only_in_df)

species_only_in_old <- setdiff(unique(eDNA_old$LCT), unique(eDNA_new$LCT))
print(species_only_in_old)


###
#whats the difference between detections and eDNA dataset 
detections <-  read.csv(here::here("Processed_data",
                                 "datasets",
                                 "detections_all_A.csv"), 
                      head=TRUE)

eDNA_df <-  read.csv(here::here("Processed_data",
                                "eDNA",
                                "datasets",
                                "eDNA_allsets_analysisA.csv"), #has eDNA index reads
                     head=TRUE)

detections <- detections %>%
  filter(gamma_detection_method %in% c('only eDNA', 'both eDNA/trawl'))


common_species <- intersect(unique(eDNA_df$LCT), unique(detections$LCT))
print(common_species)

species_only_in_df <- setdiff(unique(eDNA_df$LCT), unique(detections$LCT))
print(species_only_in_df)

species_only_in_old <- setdiff(unique(detections$LCT), unique(eDNA_df$LCT))
print(species_only_in_old)


####
detections <-  read.csv(here::here("Processed_data",
                                 "datasets",
                                 "detections_all_A.csv"), 
                      head=TRUE)

habitat <- df


species_only_in_hab <- setdiff(unique(habitat$LCT), unique(detections$LCT))
print(species_only_in_hab)

species_only_in_det <- setdiff(unique(detections$LCT), unique(habitat$LCT))
print(species_only_in_det)

##
#comparing gamma diversity (for figure 1)

old <-  read.csv(here::here("Processed_data",
                                   "datasets",
                                   "2021",
                                   "diversity",
                                   "gamma_spp_count_A.csv"), 
                        head=TRUE)

new <-  read.csv(here::here("Processed_data",
                            "datasets",
                            "diversity",
                            "gamma_spp_count_A.csv"), 
                 head=TRUE)


species_only_in_new <- setdiff(unique(new$LCT), unique(old$LCT))
print(species_only_in_new)

species_only_in_old <- setdiff(unique(old$LCT), unique(new$LCT))
print(species_only_in_old)
