library(dplyr)
library(tidyverse)
library(data.table)
library(readxl)
library(pROC)
library(caret)
library(ggplot2)

path_data="./data/Supplementary_Table - Study Data by Groups.csv"
path_classifier="./data/ciitizen_health_classifier.xlsx"
directory <- getwd()

data <- read_csv(path_data, show_col_types = FALSE)[, -1]

classifier <- read_excel(path_classifier)
classifier <- classifier[, grep("diagnosis_", colnames(classifier))]

# Combine one-hot encoded diagnosis columns
for (col_name in colnames(classifier)) {
  corresponding_cols <- paste0("clinical_diagnosis_", classifier[[col_name]])
  
  # Check if the columns exist in the data
  existing_cols <- corresponding_cols %in% colnames(data)
  if (any(existing_cols)) {
    data[[col_name]] <- ifelse(rowSums(data[corresponding_cols[existing_cols]]) > 0, 1, 0)
  }
}

# Group unique types
data$unique_types[data$unique_types %in% c(2, 3, 4, 5)] <- '2+'

# Ensure categorical variables are factors
data$Onset_group <- as.factor(data$Onset_group)
data$unique_types <- as.factor(data$unique_types)


# Predictors
predictors <- c("bilateral_tc", "Onset_group", "focal", "absence", "infantile", 
                "abnormal_eeg", "unique_types", "age_onset_m")

# Targets
diagnosis_types <- c("diagnosis_behavioral", "diagnosis_muscoloskeletal", 
                           "diagnosis_gastro", "diagnosis_immune", 
                           "diagnosis_neuro", "diagnosis_sensory", 
                           "diagnosis_respiratory", "diagnosis_cardio")

results_df <- data.frame()
significant_results_df <- data.frame()
for (diagnosis_type in diagnosis_types) {
  for (pred in predictors) {
    formula <- as.formula(paste0("`", diagnosis_type, "` ~", pred))
    
    # Fit the model
    logismod <- glm(formula, family = binomial(link = "logit"), data = data)
    summary_logismod <- summary(logismod)
    
    # Check if predictor categorical
    if (pred %in% c("Onset_group", "unique_types")) {
      levels <- c("(Intercept)", levels(data[[pred]])[-1])
      
      # Calculate confidence interval for each level
      for (level in levels) {
        level <- if(level == "(Intercept)") level else paste0(pred, level)
        
        beta_hat <- logismod$coefficients[level]
        se_beta_hat <- summary_logismod$coef[,"Std. Error"][level]
        zstar <- qnorm(0.975, 0, 1)
        lb <- beta_hat - zstar * se_beta_hat
        ub <- beta_hat + zstar * se_beta_hat
        odds <- exp(c(lb, ub))
        
        # Diagnostics
        roc_obj <- roc(logismod$model[, 1], fitted(logismod))
        auc_val <- auc(roc_obj)
        predicted <- ifelse(fitted(logismod) > 0.5, 1, 0)
        cm <- suppressWarnings(confusionMatrix(as.factor(predicted), as.factor(logismod$model[, 1])))
        
        # Extract and save key metrics
        coef <- beta_hat
        ci_lower <- exp(lb)
        ci_upper <- exp(ub)
        auc <- auc_val
        p_value <- summary_logismod$coef[,"Pr(>|z|)"][level]
        
        model_results <- data.frame(
          'diagnosis_type' = diagnosis_type,
          'predictor' = paste(pred, level, sep = "_"),
          'coefficient' = coef,
          'ci_lower' = ci_lower,
          'ci_upper' = ci_upper,
          'auc' = auc,
          'p_value' = p_value,
          row.names = NULL
        )
        results_df <- rbind(results_df, model_results)
        
        # If p-value < 0.05, save to significant_results_df
        if (p_value < 0.05) {
          significant_results_df <- rbind(significant_results_df, model_results)
        }
      }
    } else {
      # Confidence interval
      beta_hat <- logismod$coefficients[2]
      se_beta_hat <- summary_logismod$coef[,"Std. Error"][2]
      zstar <- qnorm(0.975, 0, 1)
      lb <- beta_hat - zstar * se_beta_hat
      ub <- beta_hat + zstar * se_beta_hat
      odds <- exp(c(lb, ub))
      
      # Diagnostics
      roc_obj <- roc(logismod$model[, 1], fitted(logismod))
      auc_val <- auc(roc_obj)
      predicted <- ifelse(fitted(logismod) > 0.5, 1, 0)
      cm <- suppressWarnings(confusionMatrix(as.factor(predicted), as.factor(logismod$model[, 1])))
      
      # Extract and save key metrics
      coef <- beta_hat
      ci_lower <- exp(lb)
      ci_upper <- exp(ub)
      auc <- auc_val
      p_value <- summary_logismod$coef[,"Pr(>|z|)"][2]
      
      model_results <- data.frame(
        'diagnosis_type' = diagnosis_type,
        'predictor' = pred,
        'coefficient' = coef,
        'ci_lower' = ci_lower,
        'ci_upper' = ci_upper,
        'auc' = auc,
        'p_value' = p_value
      )
      results_df <- rbind(results_df, model_results)
      
      # If p-value < 0.05, save to significant_results_df
      if (p_value < 0.05) {
        significant_results_df <- rbind(significant_results_df, model_results)
      }
    }
    # Save model summary to .txt file
    summary_file <- paste0("./results/diagnosis_type/model_summary_", diagnosis_type, "_", pred, ".txt")
    capture.output(summary_logismod, file = summary_file)
  }
}
write.csv(results_df, file = "./results/diagnosis_type/model_results.csv", row.names = FALSE)
write.csv(significant_results_df, file = "./results/diagnosis_type/significant_model_results.csv", row.names = FALSE)

# Create directory to save plots
dir.create("./results/diagnosis_type/plots")

# Create odds-ratio plots for predictors
for (predictor in unique(results_df$predictor)) {
  # Subset the data for the current predictor
  predictor_df <- results_df[results_df$predictor == predictor, ]
  
  # Determine n
  obs = c()
  for (type in predictor_df$diagnosis_type) {
    if (grepl("Onset_group", predictor)) {
      if (predictor == "Onset_group_(Intercept)") {
        subset_data <- data[data[["Onset_group"]] == "Onset_group_1-4 Months" & data[[type]] == 1, ]
        obs = append(obs, nrow(subset_data))
      } else {
        subset_data <- data[data[["Onset_group"]] == predictor & data[[type]] == 1, ]
        obs = append(obs, nrow(subset_data))
      }
    } else if (grepl("unique_types", predictor)) {
      if (predictor == "unique_types_(Intercept)") {
        subset_data <- data[data[["unique_types"]] == 0 & data[[type]] == 1, ]
        obs = append(obs, nrow(subset_data))
      } else {
        if (grepl("2+", predictor)) {
          subset_data <- data[data[["unique_types"]] == "2+" & data[[type]] == 1, ]
          obs = append(obs, nrow(subset_data))
        } else {
          subset_data <- data[data[["unique_types"]] == substring(predictor, nchar(predictor)) & data[[type]] == 1, ]
          obs = append(obs, nrow(subset_data))
        }
      }
    } else if (predictor == "age_onset_m") {
      subset_data <- data[data[[type]] == 1, ]
      obs = append(obs, nrow(subset_data))
    } else {
      subset_data <- data[data[[predictor]] == 1 & data[[type]] == 1, ]
      obs = append(obs, nrow(subset_data))
    }
  }
  predictor_df$n <- obs
  
  # Skip if predictor_df is empty or if none of the results are significant
  if (nrow(predictor_df) == 0 || all(predictor_df$p_value >= 0.05)) {
    next
  }
  
  # Calculate the midpoint of the error bars
  predictor_df$midpoint <- with(predictor_df, (ci_lower + ci_upper) / 2)
  
  # Color based on midpoint value
  predictor_df$color <- ifelse(predictor_df$midpoint < 1, "#a30234", "#0076c0")
  
  predictor_df$color <- ifelse(predictor_df$p_value > 0.05, "transparent", predictor_df$color)
  
  # Plot
  p <- ggplot(predictor_df, aes(x = diagnosis_type, y = midpoint, color = color)) +
    geom_point() +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_color_identity() +
    coord_flip(clip = "off") +
    labs(title = predictor, x = "Diagnosis Type", y = "Odds Ratio") +
    theme_bw()
  
  p <- p + geom_text(aes(y = -Inf, label = paste("n =", as.character(n))), vjust = 2.5, hjust = 1.75, size=3)
  
  if (predictor %in% c("Onset_group_Onset_group4-7 Months", "Onset_group_Onset_group7-12 Months", "Onset_group_Onset_group12+ Months")) {
    p <- p + ylim(0,1)
  } else if (predictor == "age_onset_m") {
    p <- p + ylim(0.925,1)
  } else if (predictor == "focal") {
    p <- p + ylim(0,13)
  } else if (predictor == "unique_types_unique_types2+") {
    p <- p + ylim(0,61)
  }
  
  ggsave(paste0("./results/diagnosis_type/plots/", predictor, ".pdf"), plot = p)
}