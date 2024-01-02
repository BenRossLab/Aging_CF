## Functions to recreate Supplemental Figures 3B-E
tidy_metaphlan <- function(metaphlan_df, metadata_df, Study = "Dartmouth"){
  library(tidyverse)
  cohort_samples <- filter(metadata_df, Cohort == Study) %>% select(Sample)
  
  metaphlan_filtered <- metaphlan_df %>%
    select(clade, all_of(cohort_samples$Sample))

  df <- drop_na(metaphlan_filtered)
  
  dim(df)
  
  transpose <- data.frame(t(df))
  names(transpose) <- as.matrix(transpose[1,])
  transpose <- transpose[-1,]
  transpose$Sample <- rownames(transpose)

  merge <- merge(metadata_df, transpose, by = c("Sample"))
  merge[, 5:ncol(merge)] <- lapply(5:ncol(merge), function(x) as.numeric(as.character(merge[[x]])))
  return(merge)
} 

calc_shannon <- function(metaphlan_df){
  library(tidyverse)
  print(head(metaphlan_df))
  final_df <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("Sample", "Shannon_Div", "Type", "Age", "Cohort"))
  alpha <- (dplyr::select(metaphlan_df[c(1, 5:ncol(metaphlan_df))], -Sample) %>% 
              apply(1, function(x) -sum(log((x / sum(x))**(x / sum(x))))))
  new_rows <- data.frame(Sample = metaphlan_df[1],
                         Shannon_Div = alpha,
                         Type = metaphlan_df[4],
                         Age = metaphlan_df[2],
                         Cohort = metaphlan_df[3])
  final_df <- rbind(final_df, new_rows)
  return(final_df)
}

calc_nspecies <- function(metaphlan_df){
  library(tidyverse)
  final_df <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("Sample", "N_Species", "Type", "Age", "Cohort"))
  nspecies <- dplyr::select(metaphlan_df[c(1, 5:ncol(metaphlan_df))], -Sample) %>% 
    apply(1, function(x) sum(x > 0))
  new_rows <- data.frame(Sample = metaphlan_df[1],
                         N_Species = nspecies,
                         Type = metaphlan_df[4],
                         Age = metaphlan_df[2],
                         Cohort = metaphlan_df[3])
  final_df <- rbind(final_df, new_rows)
  return(final_df)
}

add_age_bins <- function(quant_df){
  library(tidyverse)
  age_bins <- quant_df %>%
    mutate(Age_Bin = case_when(
      Age <= 6 ~ "0-6",
      Age > 6 & Age <= 12 ~ "6-12",
      Age > 12 & Age <= 18 ~ "12-18",
      Age > 18 & Age <= 24 ~ "18-24",
      Age > 24 & Age <= 30 ~ "24-30",
      Age > 30 & Age <= 36 ~ "30-36"))
  
  final <- na.omit(age_bins)
  final$Age_Bin <- factor(final$Age_Bin, levels = c("0-6", 
                                                    "6-12", 
                                                    "12-18", 
                                                    "18-24", 
                                                    "24-30", 
                                                    "30-36"))
  return(final)
}

shannon_by_cohort <- function(metaphlan_df, metadata_df, cohorts) {
  all_results <- list()
  
  for (cohort in cohorts) {
    tidy_data <- tidy_metaphlan(metaphlan_df, metadata_df, Study = cohort)
    shannon_data <- calc_shannon(tidy_data)
    age_bins <- add_age_bins(shannon_data)
    all_results[[cohort]] <- age_bins
  }
  
  final_df <- bind_rows(all_results)
  return(final_df)
}


nspecies_by_cohort <- function(metaphlan_df, metadata_df, cohorts) {
  all_results <- list()
  
  for (cohort in cohorts) {
    tidy_data <- tidy_metaphlan(metaphlan_df, metadata_df, Study = cohort)
    nspecies_data <- calc_nspecies(tidy_data)
    age_bins <- add_age_bins(nspecies_data)
    all_results[[cohort]] <- age_bins
  }
  
  final_df <- bind_rows(all_results)
  return(final_df)
}


## Works for Shannon div and num of species
plot_boxplots <- function(quant_df, quant_variable = "Shannon_Div", 
                          y_lab = "Shannon Diversity", num_col = 6, 
                          color_list, stats_list){
  library(ggplot2)
  library(ggpubr)
  boxplot <- ggplot(quant_df, aes(Cohort, .data[[quant_variable]], fill = Cohort)) +
    geom_jitter(color = "gray20", alpha = 0.3, size = 0.5) +
    geom_boxplot(notch = FALSE, varwidth = TRUE, color = "gray20", outlier.shape = NA, alpha = 0.8, size = 0.5) +
    stat_boxplot(geom ="errorbar", width = 0.6, linetype = 1, color = "gray20") +
    labs(x = "Age (Months)", y = y_lab) +
    theme_classic() +
    facet_wrap(~ Age_Bin, ncol = num_col, strip.position = "bottom") +
    scale_fill_manual(name = "Study", values = color_list) +
    scale_color_manual(values = c("gray20", "gray20")) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.x = element_text(size = 8, color = "black"),
          axis.text=element_text(size=8, color = "black"),
          axis.title=element_text(size=10, color = "black"),
          legend.text=element_text(size=8, color = "black"))
  
  stats <- boxplot + stat_compare_means(comparisons = stats_list,
                                     method = "wilcox.test",
                                     paired = FALSE,
                                     label = "p.signif",
                                     hide.ns = FALSE,
                                     step.increase = 0.1,
                                     size = 3) +
    theme(strip.text = element_text(size = 10))
  return(stats)
}

plot_lineplots <- function(quant_df, quant_variable = "Shannon_Div", 
                           y_lab = "Shannon Diversity", color_list){
  library(ggplot2)
  library(scales)
  plot <- ggplot(quant_df) + 
    geom_jitter(aes(Age, .data[[quant_variable]], color = Cohort), alpha = 0.3, size = 0.2) +
    geom_smooth(aes(Age, .data[[quant_variable]], color = Cohort), level = 0.95, method = "loess", alpha = 0.2) +
    labs(x = "Age (Months)", y = y_lab) +
    scale_x_continuous(limits = c(-0.001, 36), breaks = pretty_breaks()) +
    theme_classic() +
    scale_color_manual(name = "Study", values = color_list) +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size = 8, color = "black"),
          axis.text=element_text(size=8, color = "black"),
          axis.title=element_text(size=10, color = "black"),
          legend.text=element_text(size=8, color = "black"))
  return(plot)
}