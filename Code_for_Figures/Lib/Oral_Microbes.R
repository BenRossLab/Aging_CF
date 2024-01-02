### Code used to generate Figure 5a
get_study_data <- function(metaphlan, study_metadata){
  library(tidyverse)
  sample_list <- study_metadata$Sample
  data <- metaphlan %>%
    select(c("clade", sample_list))
  data <- na.omit(data)
  return(data)
}

# Data wrangling
sum_oral_genera <- function(metaphlan_df, metadata_df, oral_bug_list, Study = "Dartmouth"){
  library(tidyverse)
  library(reshape2)
  data_melted <- reshape2::melt(metaphlan_df)
  data_melted$clade <- data_melted$clade %>%
    gsub(".*?g__", "", .) %>%
    gsub("\\|s__.*", "", .)
  filtered_genera <- data_melted %>%
    filter(clade %in% oral_bug_list)
  colnames(filtered_genera)[2] <- "Sample"
  merged <- merge(filtered_genera, metadata_df, by = "Sample")
  add_phase <- merged %>%
    mutate(Phase = case_when(
      Age <= 12 ~ "Developmental",
      Age > 12 & Age <= 18 ~ "Transitional",
      Age > 18 ~ "Stable",
      TRUE ~ as.character(Age)
    ))
  add_phase <- na.omit(add_phase)
  sum_genera <- add_phase %>%
    group_by(Sample, clade, Age, Phase) %>% 
    summarise_all(sum) %>%
    as.data.frame()
  sum_genera$value <- sum_genera$value/100
  sum_genera$Exp <- Study
  return(sum_genera)
}

calc_oral_fraction <- function(metaphlan_df, metadata_df, oral_bug_list, Study = "Dartmouth"){
  library(tidyverse)
  library(reshape2)
  data_melted <- reshape2::melt(metaphlan_df)
  data_melted$clade <- data_melted$clade %>%
    gsub(".*?g__", "", .) %>%
    gsub("\\|s__.*", "", .)
  filtered_genera <- data_melted %>%
    filter(clade %in% oral_bug_list)
  colnames(filtered_genera)[2] <- "Sample"
  merged <- merge(filtered_genera, metadata_df, by = "Sample")
  add_phase <- merged %>%
    mutate(Phase = case_when(
      Age <= 12 ~ "Developmental",
      Age > 12 & Age <= 18 ~ "Transitional",
      Age > 18 ~ "Stable",
      TRUE ~ as.character(Age)
    ))
  add_phase <- na.omit(add_phase)
  subset <- subset(add_phase, select = c("Sample",
                                         "Phase",
                                         "value"))
  sum_genera <- subset %>%
    group_by(Sample, Phase) %>% 
    summarise_all(sum) %>%
    as.data.frame()
  subset_sum <- subset(sum_genera, select = c("Phase",
                                              "value"))
  oral_fraction <- sum_genera %>%
    group_by(Phase) %>%
    summarise(Fraction_Oral_Bugs = (sum(value > 3.4)/n())*100) %>%
    as.data.frame()
  oral_fraction$Experiment <- Study
  return(oral_fraction)
}

# Statistics
library(lqmm)
library(dplyr)

print_model_summary <- function(phase, model) {
  cat("\n-----------------------------------------\n")
  cat("Phase:", phase, "\n")
  cat("-----------------------------------------\n")
  
  model_summary <- summary(model)
  print(model_summary)
  
  cat("\n-----------------------------------------\n\n\n")
}

extract_pvalues <- function(phase, model, df_sub, reference_group = "Dartmouth") {
  model_summary_text <- capture.output(summary(model))
  
  coefs_lines <- model_summary_text[grep("Study", model_summary_text)]
  
  pvalue_df <- data.frame(
    Phase = character(),
    group1 = character(),
    group2 = character(),
    p_value = numeric(),
    signif = character(),
    y.position = numeric(),
    Study = character(),
    stringsAsFactors = FALSE
  )
  
  for (line in coefs_lines) {
    matches <- regmatches(line, regexec("^(Study\\w+)\\s+([-+]?[0-9]*\\.?[0-9]+)\\s+([-+]?[0-9]*\\.?[0-9]+)\\s+([-+]?[0-9]*\\.?[0-9]+)\\s+([-+]?[0-9]*\\.?[0-9]+)\\s+([-+]?[0-9]*\\.?[0-9]+)", line))[[1]]
    if (length(matches) == 7) {  
      group_name <- matches[2]
      p_value <- as.numeric(matches[6])
      
      group2_name <- sub("Study", "", group_name)
      
      if (group2_name == reference_group) next
      
      max_val <- max(df_sub$value, na.rm = TRUE)  
      y_pos <- ifelse(group2_name == "DIABIMMUNE", max_val + 0.1, ifelse(group2_name == "TEDDY", max_val + 0.3, NA))
      
      signif_code <- ifelse(p_value < .001, '***', ifelse(p_value < .01, '**', ifelse(p_value < .05, '*', ifelse(p_value < .1, '.', 'ns'))))
      
      pvalue_df <- rbind(pvalue_df, data.frame(
        Phase = phase,
        group1 = reference_group,
        group2 = group2_name,
        p_value = p_value,
        signif = signif_code,
        y.position = y_pos,
        Study = reference_group
      ))
    }
  }
  
  return(pvalue_df)
}

pvalue_data <- list()

for (phase in c("Developmental", "Stable", "Transitional")) {
  df_sub <- filter(fig_5a, Phase == phase)
  
  if (nrow(df_sub) > 0) { 
    model <- lqmm(value ~ Study, random = ~ 1, group = Individual,
                  tau = c(0.5), data = df_sub, 
                  control = lqmmControl(LP_max_iter = 10000))
    
    print_model_summary(phase, model)
    
    pvalue_df <- extract_pvalues(phase, model, df_sub, reference_group = "Dartmouth")
    pvalue_data[[phase]] <- pvalue_df

    } else {
    cat("No data to model for phase:", phase, "\n")
  }
}

# Combine all p-value data frames into one data frame for plotting
pvalue_final <- do.call(rbind, pvalue_data)

# Plotting
oral_bug_boxplot <- function(df_list, pvalue_df, cols){
  library(ggplot2)
  library(ggpubr)
  plot <- ggplot(df_list, aes(Study, value, fill = Study)) +
    geom_jitter(color = "gray20", alpha = 0.6, size = 0.4) +
    geom_boxplot(notch = FALSE, varwidth = TRUE, color = "gray20", outlier.shape = NA, alpha = 0.8, size = 0.5) +
    stat_boxplot(geom ="errorbar", width = 0.6, linetype = 1, color = "gray20") +
    labs(y = "Oral Microbe Relative Abundance") +
    theme_classic() +
    facet_wrap(~Phase, ncol = 3, strip.position = "top") +
    scale_color_manual(values = c("gray20", "gray20")) +
    scale_fill_manual(values = cols) +
    scale_y_continuous(trans = "log10", labels = scales::percent) + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 8, color = "black"),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank())
  
  pvalue_plot <- plot + stat_pvalue_manual(
    data = pvalue_df,
    y = "y.position", 
    group = "Phase",
    label = "signif",
    tip.length = 0.01
  )  
  
  final_plot <- pvalue_plot + geom_hline(yintercept = 0.034, linetype = "dashed", color = "#E31A1C", size = 0.75)

  return(final_plot)
}
