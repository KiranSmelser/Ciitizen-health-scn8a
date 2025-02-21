##############################
# Seizure Index Analysis Script
##############################

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
# 2. Import & Clean Medication Data
##############################
df_med <- read_excel(path_data, sheet = "medication_aggregate")

meds_to_use <- paste(c("Adrenocorticotropin (ACTH 1-18),I-125 (TYR)", "Clonazepam", "Levetiracetam", 
                       "Phenytoin", "Oxcarbazepine", "Carbamazepine", "Phenobarbital", 
                       "Lamotrigine", "Briveracetam", "Cannabidiol", "Clobazam", "Epidiolex", 
                       "Eslicarbazepine", "Ethosuximide", "Felbamate", "Gabapentin", 
                       "Prednisolone", "Lacosamide", "Primidone", "Rufinamide", "Topiramate", 
                       "Valproate", "Vigabatrin", "Zonisamide", "Stiripentol", "Tiagabine", 
                       "Perampanel"), collapse = "|")

df_med <- df_med %>%
  mutate(medication = ifelse(grepl("ACTH", medication), "ACTH", medication),
         medication = recode(medication, Epidiolex = "Epidiolex/CBD", 
                             Cannabidiol = "Epidiolex/CBD")) %>%
  filter(grepl(meds_to_use, medication))

# Function for handling medication intervals
merge_intervals <- function(intervals_df) {
  intervals_df <- intervals_df %>%
    arrange(start_med_age)
  
  if (nrow(intervals_df) == 0) {
    return(tibble(
      start_med_age = numeric(0),
      end_med_age   = numeric(0)
    ))
  }
  
  # Initialize first interval
  current_start <- intervals_df$start_med_age[1]
  current_end   <- intervals_df$end_med_age[1]
  merged_list   <- list()
  
  if (nrow(intervals_df) > 1) {
    for (i in seq(2, nrow(intervals_df))) {
      s <- intervals_df$start_med_age[i]
      e <- intervals_df$end_med_age[i]
      
      if (s <= current_end) {
        # Overlapping intervals so merge them
        current_end <- max(current_end, e, na.rm = TRUE)
      } else {
        # No overlap so add the previous interval
        merged_list[[length(merged_list) + 1]] <- tibble(
          start_med_age = current_start,
          end_med_age   = current_end
        )
        # Resetinterval
        current_start <- s
        current_end   <- e
      }
    }
  }
  
  # Push last interval
  merged_list[[length(merged_list) + 1]] <- tibble(
    start_med_age = current_start,
    end_med_age   = current_end
  )
  
  bind_rows(merged_list)
}

# Compute medication durations per patient/medication
df_duration <- df_med %>%
  transmute(
    patient_uuid,
    medication,
    start_med_age = medication_age_days_firstDate,
    end_med_age   = medication_age_days_lastDate
  ) %>%
  distinct() %>%
  arrange(patient_uuid, medication, start_med_age, end_med_age) %>%
  group_by(patient_uuid, medication) %>%
  group_modify(~ merge_intervals(.x)) %>%
  ungroup() %>%
  
  # Remove any intervals that are 0 days or negative
  filter((end_med_age - start_med_age) > 0) %>%
  
  # Each patient/med group now has valid, non-overlapping intervals
  group_by(patient_uuid, medication) %>%
  mutate(
    interval_id = row_number(),
    intervals_for_this_med = n()
  ) %>%
  ungroup() %>%
  
  mutate(
    medication = if_else(
      intervals_for_this_med > 1,
      paste0(medication, " ", interval_id),
      medication
    )
  )

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

df_type <- df_sz %>% 
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

##############################
# 4. Compute Appointment Summary
##############################
# For weighting, each patient's first and last appointment (in months)
# These are determined from the medication, seizure, and diagnosis sheets.
df_med_apts <- read_excel(path_data, sheet = "medication_aggregate")
df_sz_apts <- read_excel(path_data, sheet = "seizure_history")
df_diagnosis_apts <- read_excel(path_data, sheet = "clinical_diagnosis")

appointment_data_all <- bind_rows(
  df_med_apts %>% select(patient_uuid, appointment_age_days = medication_age_days_firstDate),
  df_sz_apts %>% select(patient_uuid, appointment_age_days = seizure_history_age_days),
  df_diagnosis_apts %>% select(patient_uuid, appointment_age_days = clinical_diagnosis_age_days_firstDate)
) %>% distinct() %>% mutate(appointment_age_months = appointment_age_days / 30)

appointment_summary <- appointment_data_all %>%
  group_by(patient_uuid) %>%
  summarise(first_appointment = min(appointment_age_months, na.rm = TRUE),
            last_appointment  = max(appointment_age_months, na.rm = TRUE),
            .groups = "drop")

##############################
# 5. Calculate Seizure Index Comparisons
##############################
# Define a safe_mean function so that if no values are present, 0 is returned.
safe_mean <- function(x) {
  if(length(x) == 0) return(0)
  m <- mean(x, na.rm = TRUE)
  if(is.nan(m)) return(0) else return(m)
}

##############################
# Determine meds_to_keep across ALL seizure types
##############################
df_combined_all_types <- df_type %>%
  left_join(df_duration, by = "patient_uuid")

# Exclude meds with < 3-month duration
df_combined_all_types <- df_combined_all_types %>%
  filter((end_med_age - start_med_age) > 91)

# Adjust the medication start age by adding 91 days when possible
df_combined_all_types <- df_combined_all_types %>%
  mutate(
    adjusted_start_med_age = if_else(
      (start_med_age + 91) < end_med_age,
      start_med_age + 91,
      start_med_age
    )
  )

# Define logical flags
df_combined_all_types <- df_combined_all_types %>%
  mutate(
    on_med       = (age_days >= adjusted_start_med_age) & (age_days <= end_med_age),
    before_med   = age_days < start_med_age,
    after_med    = age_days > end_med_age,
    within_3months_before = age_days >= pmax(start_med_age - 91, 0) & age_days < start_med_age,
    within_3months_after  = age_days > end_med_age & age_days <= (end_med_age + 91)
  )

# Handle overlapping medication periods
df_assigned_all <- df_combined_all_types %>%
  arrange(patient_uuid, age_days, adjusted_start_med_age) %>%
  group_by(patient_uuid, age_days) %>%
  filter(on_med) %>%
  ungroup()

df_combined_all_types <- df_combined_all_types %>%
  left_join(
    df_assigned_all %>%
      select(patient_uuid, age_days, medication) %>%
      rename(assigned_med = medication),
    by = c("patient_uuid", "age_days")
  ) %>%
  group_by(patient_uuid, age_days) %>%
  mutate(on_med_final = medication %in% assigned_med) %>%
  ungroup()

# Summarize to figure out which meds have at least one seizure on-med or within 
# 3 months before/after, across ALL seizure types
meds_to_keep <- df_combined_all_types %>%
  group_by(patient_uuid, medication) %>%
  summarise(
    has_seizure_on_med      = any(on_med_final, na.rm = TRUE),
    has_seizure_before_after = any(within_3months_before | within_3months_after, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(has_seizure_on_med | has_seizure_before_after) %>%
  select(patient_uuid, medication)

seizure_types <- unique(df_type$type)
seizures_summary_list <- list()

for (s_type in seizure_types) {
  
  # Filter for current seizure type and merge with medication durations
  df_type_filtered <- df_type %>% filter(type == s_type)
  df_combined <- df_type_filtered %>%
    left_join(df_duration, by = "patient_uuid")
  
  # Exclude medications with duration less than 3 months (91 days)
  df_combined <- df_combined %>% filter((end_med_age - start_med_age) > 91)
  
  # Adjust the medication start age by adding 91 days (if possible)
  df_combined <- df_combined %>%
    mutate(
      adjusted_start_med_age = if_else((start_med_age + 91) < end_med_age,
                                       start_med_age + 91, start_med_age)
    )
  
  # Define seizure timing relative to the medication period
  df_combined <- df_combined %>%
    mutate(
      on_med = (age_days >= adjusted_start_med_age) & (age_days <= end_med_age),
      before_med = age_days < start_med_age,
      after_med = age_days > end_med_age
    )
  
  # For seizures that occur during overlapping medication periods, assign them to all active periods
  df_assigned <- df_combined %>%
    arrange(patient_uuid, age_days, adjusted_start_med_age) %>%
    group_by(patient_uuid, age_days) %>%
    filter(on_med) %>%
    ungroup()
  
  df_combined <- df_combined %>%
    left_join(df_assigned %>% select(patient_uuid, age_days, medication) %>% rename(assigned_med = medication),
              by = c("patient_uuid", "age_days"))
  
  df_combined <- df_combined %>%
    group_by(patient_uuid, age_days) %>%
    mutate(on_med_final = if_else(medication %in% assigned_med, TRUE, FALSE)) %>%
    ungroup()
  
  # Identify seizures within 3 months before or after medication periods
  df_combined <- df_combined %>%
    mutate(
      within_3months_before = age_days >= pmax(start_med_age - 91, 0) & age_days < start_med_age,
      within_3months_after  = age_days > end_med_age & age_days <= (end_med_age + 91)
    )
  
  # Retain only medication episodes with at least one seizure on med or within 3 months before/after
  df_combined <- df_combined %>%
    inner_join(meds_to_keep, by = c("patient_uuid", "medication"))
  
  # Join the first and last appointment (in months) for each patient
  df_combined <- df_combined %>% left_join(appointment_summary, by = "patient_uuid") %>%
    mutate(
      adjusted_start_med_age_months = adjusted_start_med_age / 30,
      end_med_age_months = end_med_age / 30,
      # Duration for "before" is from first appointment to the adjusted start of medication:
      duration_before = pmax(adjusted_start_med_age_months - first_appointment, 0),
      # "On" duration is between adjusted start and medication end:
      duration_on = pmax(end_med_age_months - adjusted_start_med_age_months, 0),
      # "After" duration is from medication end to the last appointment:
      duration_after = pmax(last_appointment - end_med_age_months, 0),
      duration_off = duration_before + duration_after
    )
  
  # Compute Weighted Averages
  seizures_summary <- df_combined %>%
    group_by(patient_uuid, medication, start_med_age, end_med_age,
             adjusted_start_med_age_months, end_med_age_months, first_appointment, last_appointment) %>%
    summarise(
      on_med_avg     = safe_mean(index[on_med_final]),
      before_med_avg = safe_mean(index[before_med]),
      after_med_avg  = safe_mean(index[after_med]),
      off_med_avg    = safe_mean(index[!on_med_final]),
      # Since durations are the same for all events in an episode, take the first:
      duration_before = first(duration_before),
      duration_on     = first(duration_on),
      duration_after  = first(duration_after),
      duration_off    = first(duration_off),
      .groups = "drop"
    ) %>%
    mutate(
      weighted_index_med    = if_else(duration_on > 0, on_med_avg / sqrt(duration_on), 0),
      weighted_index_before = if_else(duration_before > 0, before_med_avg / sqrt(duration_before), 0),
      weighted_index_after  = if_else(duration_after > 0, after_med_avg / sqrt(duration_after), 0),
      weighted_index_off    = if_else(duration_off > 0, off_med_avg / sqrt(duration_off), 0),
      diff_on_vs_off   = weighted_index_med - weighted_index_off,
      diff_on_vs_before = weighted_index_med - weighted_index_before,
      diff_on_vs_after  = weighted_index_med - weighted_index_after,
      type = s_type
    )
  
  seizures_summary_list[[s_type]] <- seizures_summary
}

# Combine results for all seizure types
seizures_summary_combined <- bind_rows(seizures_summary_list)

##############################
# 6. Additional Data for Patient Reports
##############################
genetics <- read_excel(path_data, sheet = "genetic_findings")
demographics <- read_excel(path_data, sheet = "demographics")
# (Re-read appointment data)
df_med_apts <- read_excel(path_data, sheet = "medication_aggregate")
df_sz_apts <- read_excel(path_data, sheet = "seizure_history")
df_diagnosis_apts <- read_excel(path_data, sheet = "clinical_diagnosis")
overlap_patients <- read_excel(path_overlap_patients)
hospitalizations <- read_excel(path_data, sheet = "hospital_admission")

# Filter adverse effects for Moderate and Severe events
adverse_effects <- read_excel(path_data, sheet = "adverse_effects")
adverse_effect_severity <- read_excel("./data/effects_severity.xlsx")
adverse_effects <- adverse_effects %>%
  inner_join(
    adverse_effect_severity %>% filter(severity_score %in% c("Moderate", "Severe")), 
    by = c("adverse_effect" = "adverse_effect")
  )

##############################
# 7. Heatmap Function
##############################
create_heatmap <- function(data, fill_var, title_suffix, limits = NULL) {
  ggplot(data, aes(x = type, y = medication, fill = !!sym(fill_var))) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(!!sym(fill_var), 2)), size = 3, color = "black", na.rm = TRUE) +
    scale_fill_gradient2(
      low = "#083681",
      mid = "#F7F7F7",
      high = "#C80813FF",
      midpoint = 0,
      na.value = "#F7F7F7",
      limits = limits
    ) +
    #scale_fill_distiller(palette = "RdBu", limits = limits) +
    # scale_fill_gradientn(
    #   colors = hcl.colors(3, palette = "Blue-Red"),
    #   limits = limits,
    #   na.value = "#F7F7F7"
    # ) +
    #scale_fill_viridis(discrete = FALSE, option = "E") +
    theme_classic() +
    labs(
      title = title_suffix,
      x = "Seizure Type",
      y = "Medication",
      fill = "Change in Seizure Index"
    )
}

##############################
# 8. Generate Patient Reports (Timelines, Heatmaps, & Line Plots)
##############################
# Save patient charts to single PDF
pdf_file <- "./figures/patients/combined_patients_report.pdf"
pdf(pdf_file, width = 16, height = 12)

#Import developmental data
df_dev <- read_excel(path_data, sheet='development')

for (pt in unique(seizures_summary_combined$patient_uuid)) {
  
  # Retrieve the patient’s protein mutation (from genetics or overlap file)
  protein_mutation <- genetics %>%
    filter(patient_uuid == pt) %>%
    pull(variant_protein) %>%
    first()
  if (is.na(protein_mutation)) {
    protein_mutation <- overlap_patients %>%
      filter(`Patient ID` == pt) %>%
      pull(`p.`) %>%
      first()
  }
  protein_mutation <- ifelse(is.na(protein_mutation), "NA", protein_mutation)
  
  # Check if the patient has a report of "Infantile spasms"
  has_infantile_spasms <- df_sz %>%
    filter(patient_uuid == pt, type == "Infantile spasms") %>%
    nrow() > 0
  title_suffix <- if (has_infantile_spasms) " - IF" else ""
  
  # Subset seizure summary data for the current patient
  pt_data <- seizures_summary_combined %>% filter(patient_uuid == pt)
  
  # Prepare medication duration data (convert days to months)
  pt_data_duration <- df_duration %>%
    filter(patient_uuid == pt) %>%
    mutate(
      start_med_age_months = start_med_age / 30,
      end_med_age_months   = end_med_age / 30,
      first_3_months_end   = pmin(start_med_age_months + 3, end_med_age_months),
      # Add a base name for plotting
      medication_base = sub(" \\d+$", "", medication)
    )
  
  # Compute total duration per medication for this patient (sum across intervals if needed)
  duration_order <- pt_data_duration %>%
    mutate(total_duration = end_med_age - start_med_age) %>%
    group_by(medication_base) %>%
    summarise(total_duration = sum(total_duration, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_duration))
  
  # Use 'medication_base' to track which medication it maps to
  pt_data <- pt_data %>%
    mutate(medication_base = sub(" \\d+$", "", medication))
  
  # Factor 'medication' in descending order
  pt_data <- pt_data %>%
    left_join(duration_order, by = "medication_base") %>%
    arrange(total_duration) %>%
    mutate(
      medication = factor(medication, levels = unique(medication))
    )
  
  # Create heatmaps for different seizure index comparisons
  heatmap_on_vs_off <- create_heatmap(pt_data, "diff_on_vs_off", "On vs. Off")
  heatmap_on_vs_before <- create_heatmap(pt_data, "diff_on_vs_before", "On vs. Before") +
    theme(legend.position = "none")
  
  # Compute legend limits for the on vs. after heatmap from the diff columns
  # legend_limits <- c(
  #   min(c(pt_data$diff_on_vs_after, pt_data$diff_on_vs_before, -0.5), na.rm = TRUE),
  #   max(c(pt_data$diff_on_vs_after, pt_data$diff_on_vs_before), na.rm = TRUE)
  # )
  legend_limits <- max(abs(pt_data$diff_on_vs_after), abs(pt_data$diff_on_vs_before), na.rm = TRUE) * c(-1, 1)
  
  
  heatmap_on_vs_after <- create_heatmap(pt_data, "diff_on_vs_after", "On vs. After", limits = legend_limits) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  # Prepare data for the smooth line plot (seizure index over time)
  pt_data_type <- df_type %>%
    filter(patient_uuid == pt) %>%
    mutate(age_months = age_days / 30)
  
  p_smooth <- ggplot(pt_data_type, aes(x = age_months, y = index, color = type)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_line(linetype = "dashed", alpha = 0.8) +
    facet_grid(type ~ .) +
    theme_light() +
    labs(x = "Age (months)", y = "Seizure Index", color = "Seizure Type") +
    guides(color="none")
  
  # Collect unique appointment dates for the timeline
  appointment_data <- bind_rows(
    df_med_apts %>% filter(patient_uuid == pt) %>% select(appointment_age_days = medication_age_days_firstDate),
    df_sz_apts %>% filter(patient_uuid == pt) %>% select(appointment_age_days = seizure_history_age_days),
    df_diagnosis_apts %>% filter(patient_uuid == pt) %>% select(appointment_age_days = clinical_diagnosis_age_days_firstDate)
  ) %>%
    distinct() %>%
    mutate(appointment_age_months = appointment_age_days / 30)
  
  # Determine most recent appointment
  pt_demographics <- demographics %>%
    filter(patient_uuid == pt) %>%
    mutate(most_recent_record_age_months = most_recent_records_age_days / 30)
  
  # Adverse effects (convert age to months)
  pt_data_adverse <- adverse_effects %>%
    filter(patient_uuid == pt) %>%
    mutate(age_months = adverse_effect_age_days_firstDate / 30)
  
  # Hospitalizations for "Status epilepticus"
  pt_data_status <- hospitalizations %>%
    filter(patient_uuid == pt, admission_diagnosis == "Status epilepticus") %>%
    mutate(age_months = admission_age_days_firstDate / 30)
 
  ####################################
  
  #development milestones
  
  # add developments to a line on the plot
  # filter by patient within for loop right below
  # group by milestone, determine aquired vs regressed milestones
  # where regressed is defined as able -> unable
  # aquired is able -> able    no change          
  # gained = green not digressed = red never gained = not on plot
  
  pt_dev_data <- df_dev %>%
    filter(patient_uuid == pt)%>%
    group_by(domain_milestone)%>%
    arrange(domain_age_days_firstDate)%>%
    mutate(
      age_months = domain_age_days_firstDate/30,
      status_change = case_when(
        domain_status == 'Unable' & lag(domain_status) == 'Able' ~ 'Loss',
        domain_status == 'Able' & lag(domain_status) == 'Unable' ~ 'Gained',
        TRUE ~ NA_character_
      )
    )%>%
    ungroup()
  
  pt_dev_data <- pt_dev_data %>%
    filter(!is.na(status_change)) %>%
    mutate(status_change = factor(status_change, levels = c("Gained", "Loss")))
  
  
 #####################################
  
  # Compute earliest start age per medication
  pt_data_duration <- pt_data_duration %>%
    group_by(medication_base) %>%
    mutate(earliest_start = min(start_med_age_months)) %>%
    ungroup()
  
  # Order medications by start age in ascending order
  med_order <- pt_data_duration %>%
    distinct(medication_base, earliest_start) %>%
    arrange(earliest_start) %>%
    pull(medication_base)
  
  # Timeline plot: medication periods, seizure events, adverse effects, status events, and appointments
  p_timeline <- ggplot() +
    # Medication timeline: first 3 months in gray, remainder in blue
    geom_segment(
      data = pt_data_duration,
      aes(x = start_med_age_months, xend = first_3_months_end,
          y = medication_base, yend = medication_base),
      size = 2, color = "#8A9197FF"
    ) +
    geom_segment(
      data = pt_data_duration,
      aes(x = first_3_months_end, xend = end_med_age_months,
          y = medication_base, yend = medication_base),
      size = 2, color = "#709AE1FF"
    ) +
    # Seizure events
    geom_point(
      data = pt_data_type,
      aes(x = age_months, y = type),
      color = "#C80813FF", size = pt_data_type$index + 1, alpha = 0.6
    ) +
    # Adverse effects
    geom_point(
      data = pt_data_adverse,
      aes(x = age_months, y = "Adverse Effects"),
      color = "#FD7446FF", size = 2, shape = 15, alpha = 0.8
    ) +
    # Status epilepticus events
    geom_point(
      data = pt_data_status,
      aes(x = age_months, y = "Status epilepticus"),
      color = "#FED439FF", size = 5, shape = 18, alpha = 0.9
    ) +
    # Development Loss and Gains
    geom_point(
      data = pt_dev_data,
      aes(x = age_months, y='Developmental Milestones',
      color = status_change, shape = domain), alpha = 0.7, size = 3
    ) +
    
    scale_color_manual(values = c("Loss" = "red", "Gained" = "green3"), guide = 'none') +
    
    scale_shape_discrete(name = "Developmental Domain") +
    
    # Appointment markers
    geom_point(
      data = appointment_data,
      aes(x = appointment_age_months, y = "Appointments"),
      color = "#1A9993FF", size = 3, shape = 17, alpha = 0.6
    ) +
    # Most recent appointment marker
    # geom_point(
    #   data = pt_demographics,
    #   aes(x = most_recent_record_age_months, y = "Appointments"),
    #   color = "#370335FF", size = 3, shape = 17, alpha = 0.6
    # ) +
    theme_linedraw() +
    labs(
      title = paste(pt, " (", protein_mutation, ")", title_suffix, sep = ""),
      x = "Age (months)",
      y = ""
    ) +
    scale_y_discrete(limits = c("Appointments", "Developmental Milestones",
                                "Status epilepticus",
                                rev(unique(pt_data_type$type)),
                                "Adverse Effects",
                                rev(med_order)))
  
  # Combine the plots
  combined_plot <- (p_smooth | heatmap_on_vs_before | heatmap_on_vs_after) /
    p_timeline +
    plot_layout(heights = c(1, 1))
  
  print(combined_plot)
}

dev.off()