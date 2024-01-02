## Code to recreate Supplemental Figure 1
phyla_needed <- c("Actinobacteria",
                  "Bacteroidetes",
                  "Firmicutes",
                  "Proteobacteria")

get_phyla <- function(metaphlan_df, metadata_df, phyla_oi, age_max = 36){
  library(tidyverse)
  melted <- pivot_longer(metaphlan_df,
                         cols = -clade,
                         names_to = "Sample",
                         values_to = "Relab")
  melted$clade <- melted$clade %>%
    gsub(".*?p__", "", .) %>%
    gsub("\\|c__.*", "", .)
  filtered_phyla <- melted %>%
    filter(clade %in% phyla_oi)
  filtered_phyla <- na.omit(filtered_phyla)
  sum_phyla <- filtered_phyla %>%
    group_by(Sample, clade) %>%
    summarise_all(sum) %>%
    as.data.frame()
  merged <- merge(sum_phyla, metadata_df, by = "Sample")
  age_filter <- merged %>%
    filter(Age <= age_max)
  age_bins <- age_filter %>%
    mutate(Age_Bin = case_when(
      Age <= 6 ~ "0-6",
      Age > 6 & Age <= 12 ~ "6-12",
      Age > 12 & Age <= 18 ~ "12-18",
      Age > 18 & Age <= 24 ~ "18-24",
      Age > 24 & Age <= 30 ~ "24-30",
      Age > 30 & Age <= 36 ~ "30-36",
      TRUE ~ as.character(Age)
    ))
  age_bins <- na.omit(age_bins)
  final_df <- subset(age_bins, select = c("Sample", "clade", "Age_Bin", "Relab", "Cohort", "Type"))
  return(final_df)
}

phyla_boxplot <- function(phyla_df, test_list, color_list){
  library(ggplot2)
  library(ggpubr)
  plot <- ggplot(phyla_df, aes(Cohort, Relab, fill = Cohort)) +
    geom_jitter(color = "gray20", alpha = 0.6, size = 0.4) +
    geom_boxplot(notch = FALSE, varwidth = TRUE, color = "gray20", outlier.shape = NA, alpha = 0.8, size = 0.5) +
    stat_boxplot(geom ="errorbar", width = 0.6, linetype = 1, color = "gray20") +
    labs(y = "Relative Abundance") +
    theme_classic() +
    facet_wrap(clade~Age_Bin, ncol = 6, strip.position = "top") +
    scale_color_manual(values = c("gray20", "gray20")) +
    scale_fill_manual(values = color_list) +
    theme(axis.text.x = element_text(angle = 90, color = "black", size = 8), axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank())
  all_stat <- plot + stat_compare_means(comparisons = test_list, 
                                        method = "wilcox.test", 
                                        paired = FALSE, 
                                        label = "p.signif", 
                                        hide.ns = FALSE, 
                                        step.increase = 0.1, 
                                        size = 3) +
    theme(strip.text = element_text(size = 10, color = "black"))
  
  return(all_stat)
  
}

process_by_cohort <- function(metaphlan_df, metadata_df, phyla_oi, cohort_list) {
  library(tidyverse)
  process_cohort <- function(cohort) {
    cohort_metadata_df <- filter(metadata_df, Cohort == cohort)
    get_phyla(metaphlan_df, cohort_metadata_df, phyla_oi)
  }
  results <- lapply(cohort_list, process_cohort)
  final_df <- bind_rows(results)
  final_df$Age_Bin <- factor(final_df$Age_Bin, levels = c("0-6",
                                                          "6-12",
                                                          "12-18",
                                                          "18-24",
                                                          "24-30",
                                                          "30-36"))
  
  
  return(final_df)
}
