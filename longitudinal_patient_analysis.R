
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

view(df_table_data)
#############################
#Look at every patient and then count the number of seizure events
#############################

seizure_counts <- df_table_data %>%
  group_by(patient_uuid)%>%
  summarise(seizure_count = n(),
            mean_index = round(mean(index),2),
            number_seizure_types = n_distinct(type),
            seizure_types = paste(unique(type), collapse = ", ") )

view(seizure_counts)


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
  select(patient_uuid,seizure_count,mean_index,number_seizure_types, seizure_types,gap,gap_period)

view(combined_seizure_data)


##############################
#Medication Data
##############################
 
appointment_data_all <- bind_rows(
  seizures_summary_combined %>% select(patient_uuid, appointment_age_days = start_med_age),
  seizures_summary_combined %>% select(patient_uuid, appointment_age_days = end_med_age)
) %>%
  distinct() %>%
  mutate(appointment_age_months = appointment_age_days / 30)

appointment_summary <- appointment_data_all %>%
  group_by(patient_uuid) %>%
  summarise(last_appointment = max(appointment_age_months, na.rm = TRUE), .groups = "drop")

view(seizures_summary_combined)

length(unique(seizures_summary_combined$patient_uuid))

medication_df <- seizures_summary_combined %>%
  group_by(patient_uuid) %>%
  summarise(number_med_types = n_distinct(medication),
            med_types = paste(unique(medication), collapse = ", "))

view(medication_df)

pt <- unique(seizures_summary_combined$patient_uuid)

combined_seizure_data <- combined_seizure_data %>%
  filter(patient_uuid %in% pt)

combined_seizure_data <- left_join(combined_seizure_data,medication_df, by='patient_uuid')

view(combined_seizure_data)

# Determine current vs. weened medications
current_weaned_df <- df_duration %>%
  mutate(
    end_age_months = end_med_age / 30,
    start_age_months = start_med_age / 30
  ) %>%
  left_join(appointment_summary, by = "patient_uuid") %>%
  group_by(patient_uuid, medication) %>%
  summarise(
    is_current = any(end_age_months >= last_appointment, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(patient_uuid) %>%
  summarise(
    current_medications = sum(is_current),
    weened_medications = sum(!is_current),
    .groups = "drop"
  )

#Meds in the gap
#look at the seizures_summary data and find the age in months
#round it then see how it corresponds with gap start and end
#return a list of number_meds_gap and meds_gap

view(df_combined)

# Assuming df_table_data contains seizures and df_medications contains medication details
meds_during_gap <- seizure_gaps %>%
  inner_join(seizures_summary_combined, by = "patient_uuid") %>%  # Join on patient ID
  filter(
    round(start_med_age/30) <= End_Age,  # Medication started before or during gap end
    round(end_med_age/30) >= Start_Age   # Medication ended after or during gap start
  ) %>%
  select(patient_uuid, medication, start_med_age, end_med_age, gap_period) %>%
  mutate(start_in_months = round(start_med_age/30),
         end_in_months = round(end_med_age/30)) %>%
  distinct()

view(meds_during_gap)


gap_meds_join <- meds_during_gap %>%
  group_by(patient_uuid) %>%
  summarise(number_med_types_gap = n_distinct(medication),
            med_types_gap = paste(unique(medication), collapse = ", "))

view(gap_meds_join)

anti <- anti_join(combined_seizure_data,gap_meds_join, by='patient_uuid')
view(anti)

combined_df <- combined_seizure_data %>%
  left_join(gap_meds_join, by = "patient_uuid")


combined_df <- combined_df %>%
  mutate(
    number_med_types_gap = as.character(number_med_types_gap),  # Convert to character
    number_med_types_gap = replace_na(number_med_types_gap, "None"),  # Replace NA with "None"
    med_types_gap = replace_na(med_types_gap, "None"),  # Replace NA in med_types_gap
    gap_period = gsub("-", " to ", gap_period)
  )

# Join current_weaned_df to combined_df
combined_df <- combined_df %>%
  left_join(current_weaned_df, by = "patient_uuid") %>%
  mutate(
    current_medications = replace_na(current_medications, 0),
    weened_medications = replace_na(weened_medications, 0)
  )

view(combined_df)

#write.csv(combined_df, "combined_longitudinal_table.csv", row.names = FALSE)
