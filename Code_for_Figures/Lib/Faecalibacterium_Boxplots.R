### Code used to generate Figure 6d
# Data wrangling from strain assignment table
calc_n_obs <- function(fungal_df){
  phases <- list("Developmental",
                 "Transitional",
                 "Stable")
  for(i in phases){
    filtered <- fungal_df %>%
      filter(phase == i)
    n <- length(unique(filtered$Sample))
    print(n)
  }
}

# Statistics
library(lme4)
library(broom)

all_meta <- read.csv("../Fungal_plots/all_meta_for_Adrians_code.csv", stringsAsFactors = T)
colnames(all_meta)[2] <- "Age"

fig_6d <- merge(all_clades, all_meta[c(1,3)])

# Print the phase and summary of the model
print_model_summary <- function(phase, model) {
  cat("\n-----------------------------------------\n")
  cat("Phase:", phase, "\n")
  cat("-----------------------------------------\n")
  
  model_summary <- summary(model)
  print(model_summary)
  
  cat("\n-----------------------------------------\n\n\n")
}

# Your existing function to print the model summary remains the same
extract_pvalues <- function(phase, model, df_sub, reference_group = "Dartmouth") {
  library(dplyr)
  # Extract the summary for fixed effects
  model_summary <- summary(model)
  
  # Retrieve the fixed effect names and p-values
  fixed_effects_names <- rownames(model_summary$coefficients)
  p_values <- model_summary$coefficients[, "Pr(>|z|)"]
  
  # Prepare the group names by removing the 'Study' prefix
  group2_names <- sub("Cohort", "", fixed_effects_names)
  
  max_val <- max(df_sub$n_strains, na.rm = TRUE)
  
  # Create a dataframe for p-values
  pvalue_df <- data.frame(
    phase = rep(phase, length(group2_names)),
    group1 = reference_group,
    group2 = group2_names,
    p_value = p_values,
    signif = NA, # Placeholder, will calculate this next
    y.position = NA, # Placeholder, will set based on group2
    Cohort = reference_group
  )
  
  # Calculate significance codes based on p-values
  pvalue_df$signif <- ifelse(pvalue_df$p_value < .001, '***',
                             ifelse(pvalue_df$p_value < .01, '**',
                                    ifelse(pvalue_df$p_value < .05, '*',
                                           ifelse(pvalue_df$p_value < .1, '.', 'ns'))))
  
  # Set y.position based on group2 values
  pvalue_df$y.position <- ifelse(group2_names == "DIABIMMUNE", max_val + 0.5, 
                                 ifelse(group2_names == "TEDDY", max_val + 1, NA))  
  
  # Filter out the intercept (if present) and the reference group itself
  pvalue_df <- pvalue_df[!pvalue_df$group2 %in% c("(Intercept)", reference_group), ]
  
  return(pvalue_df)
}

# Initialize a list to store the p-value data frames for each phase
pvalue_data <- list()

for (phase in c("Developmental", "Stable","Transitional")) {
  df_sub <- filter(fig_6d, phase == phase)
  
  model <- glmer("n_strains ~ Cohort + (1|Individual)", data = df_sub, family = poisson, control = glmerControl(optimizer="bobyqa"))
  
  # Print model summary
  print_model_summary(phase, model)
  
  # Extract p-values and create a data frame for stat_pvalue_manual
  pvalue_df <- extract_pvalues(phase, model, df_sub)
  pvalue_data[[phase]] <- pvalue_df
}

# Combine all p-value data frames into one
pvalue_data <- do.call(rbind, pvalue_data)

## Plotting
faecalibacterium_boxplot <- function(df_list, pvalue_df, cols){
  library(ggplot2)
  library(ggpubr)
  plot <- ggplot(df_list, aes(Cohort, n_strains, fill = Cohort)) +
    geom_jitter(color = "gray20", alpha = 0.6, size = 0.4) +
    geom_boxplot(notch = FALSE, varwidth = TRUE, color = "gray20", outlier.shape = NA, alpha = 0.8) +
    stat_boxplot(geom ="errorbar", width = 0.6, linetype = 1, color = "gray20") +
    facet_wrap(~phase, ncol = 3) +
    labs(y = "Number of Faecalibacterium Clades") +
    theme_classic() +
    scale_y_continuous(limits = c(0,10)) +
    scale_color_manual(values = c("gray20", "gray20")) +
    scale_fill_manual(values = cols, name = "Study") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 8, color = "black"))
  
  pvalue_plot <- plot + stat_pvalue_manual(
    data = pvalue_df,
    y = "y.position", 
    group = "phase",
    label = "signif",
    tip.length = 0.01
  )  
  
  return(pvalue_plot)
}
