
library(dplyr)
library(tidyverse)
library(readxl)
library(ggplot2)
library(ggsci)
library(lubridate)
library(patchwork)
library(viridis)

##############################
# 1. File Paths
##############################
path_data <- "./data/Ciitizen_SCN8A_UArizona_2024.02.09.xlsx"
path_classifier <- "./data/ciitizen_health_classifier.xlsx"
path_tc_index <- "./data/tonic-clonic_index.xlsx"
path_focal_index <- "./data/focal_index.xlsx"
path_myoclonic_index <- "./data/myoclonic_index.xlsx"
path_absence_index <- "./data/absence_index.xlsx"
path_tonic_index <- "./data/tonic_index.xlsx"
path_overlap_patients <- "./data/Overlap Patients corrected.xlsx"


##############################
# 3. Import & Clean Seizure Data
##############################
# Helper functions for cleaning symbols in seizure index values
convert_greater_than <- function(x) {
  if (grepl(">", x)) {
    x <- as.numeric(sub(">", "", x)) + 1
  }
  return(x)
}
convert_less_than <- function(x) {
  if (grepl("<", x)) {
    x <- as.numeric(sub("<", "", x)) - 1
  }
  return(x)
}
convert_greater_than_equal_to <- function(x) {
  if (grepl("≥", x)) {
    x <- as.numeric(sub("≥", "", x))
  }
  return(x)
}
convert_less_than_equal_to <- function(x) {
  if (grepl("≤", x)) {
    x <- as.numeric(sub("≤", "", x))
  }
  return(x)
}

# Import seizure history data and clean the seizure value column
df_sz <- read_excel(path_data, sheet = "seizure_history")
df_sz$seizure_history_value[is.na(df_sz$seizure_history_value)] <- 1
df_sz$seizure_history_value <- sapply(df_sz$seizure_history_value, convert_greater_than)
df_sz$seizure_history_value <- sapply(df_sz$seizure_history_value, convert_greater_than_equal_to)
df_sz$seizure_history_value <- sapply(df_sz$seizure_history_value, convert_less_than)
df_sz$seizure_history_value <- sapply(df_sz$seizure_history_value, convert_less_than_equal_to)
df_sz$seizure_history_value <- as.numeric(str_trim(df_sz$seizure_history_value))
names(df_sz) <- sub('^seizure_history_', '', names(df_sz))

# Import seizure classifier and restrict to the relevant seizure types
classifier <- read_excel(path_classifier)
classifier <- classifier[grepl("seizure_", names(classifier))]
names(classifier) <- sub('^seizure_', '', names(classifier))

df_sz <- df_sz %>% 
  filter(type %in% c(classifier$`tonic-clonic`, classifier$focal, 
                     classifier$absence, classifier$tonic, classifier$myoclonic)) %>%
  mutate(type = case_when(
    type %in% classifier$`tonic-clonic` ~ "Tonic-clonic",
    type %in% classifier$focal ~ "Focal",
    type %in% classifier$absence ~ "Absence",
    type %in% classifier$tonic ~ "Tonic",
    type %in% classifier$myoclonic ~ "Myoclonic",
    TRUE ~ type
  ))

# Import seizure index data from several sheets and combine them
tc_index <- read_excel(path_tc_index) %>% mutate(type = "Tonic-clonic")
focal_index <- read_excel(path_focal_index) %>% mutate(type = "Focal")
absence_index <- read_excel(path_absence_index) %>% mutate(type = "Absence")
tonic_index <- read_excel(path_tonic_index) %>% mutate(type = "Tonic")
myoclonic_index <- read_excel(path_myoclonic_index) %>% mutate(type = "Myoclonic")

index <- bind_rows(tc_index, focal_index, absence_index, tonic_index, myoclonic_index)

names(index) <- sub('^seizure_history_', '', names(index))
names(index) <- sub('^seizure_', '', names(index))
index <- subset(index, select = -c(3:5, 7))
index$value <- as.numeric(index$value)

# Merge index values into seizure history; set any missing index to 1
df_type <- left_join(df_type, unique(index), by = c("type", "value", "unit"), 
                     relationship = "many-to-many")
df_type$index[is.na(df_type$index)] <- 1

# Remove LOF patients
df_type <- df_type %>% filter(!patient_uuid %in% classifier$subgroup_lof)
df_duration <- df_duration %>% filter(!patient_uuid %in% classifier$subgroup_lof)

df_table_data <- df_type %>%
  select(patient_uuid,type,index,age_days)%>%
  mutate(age_in_months = age_days/30)

#############################
#Look at every patient and then count the number of seizure events
#############################

seizure_counts <- df_table_data %>%
  group_by(patient_uuid)%>%
  summarise(seizure_count = n(),
            mean_index = round(mean(index),2),
            seizure_types = n_distinct(type))

view(seizure_counts)

df_table_data %>%
  filter(patient_uuid == 'dc15a4e2-a7f1-42be-8418-d5da4f2abe40')%>%
  summarise(mean_index = mean(index))

##############################
#calculate the gaps
##############################

seizure_gaps <- df_table_data %>%
  group_by(patient_uuid) %>%  
  arrange(age_in_months) %>%  # Ensure data is sorted by time
  summarise(
    gap = ifelse(n() > 1, max(diff(age_in_months), na.rm = TRUE), NA_real_),
    Start_Age = ifelse(n() > 1, age_in_months[which.max(diff(age_in_months))], NA_real_),  
    End_Age = ifelse(n() > 1, age_in_months[which.max(diff(age_in_months)) + 1], NA_real_)
  ) %>%
  mutate(gap = replace_na(gap, 0),
         Start_Age = round(Start_Age),  # Round start age
         End_Age = round(End_Age),
         gap_period = ifelse(gap > 0, paste0(Start_Age, " - ", End_Age), "No gap"))

##############################
#join the dfs
##############################

combined_seizure_data <- seizure_counts %>%
  left_join(seizure_gaps, by = "patient_uuid")%>%
  mutate(gap = round(gap))

combined_seizure_data <- combined_seizure_data %>%
  select(patient_uuid,seizure_count,mean_index,seizure_types,gap,gap_period)

view(combined_seizure_data)


##############################
#Medication Data
##############################