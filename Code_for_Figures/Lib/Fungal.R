# Functions used to create fungal plots in Figure 5B and 5C
## Data wrangling 
calc_n_obs <- function(fungal_df){
  library(tidyverse)
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

metaphlan_fungal_tidy <- function(metaphlan_df, metadata_df, study = "Dartmouth"){
  library(tidyverse)
  library(reshape2)
  fungal_filter <- metaphlan_df %>%
    filter(grepl("Eukaryota", clade))
  melt <- reshape2::melt(fungal_filter)
  colnames(melt)[2] <- "Sample"
  merge <- merge(melt, metadata_df, by = "Sample")
  bin <- merge %>%
    mutate(phase = case_when(
      Age <= 12 ~ "Developmental",
      Age > 12 & Age <= 18 ~ "Transitional",
      Age > 18 ~ "Stable",
      TRUE ~ as.character(Age)
    ))
  bin$Study <- study
  bin_sub <- subset(bin, select = c("Sample", "value", "phase", "Study", "clade"))
  bin_sub <- na.omit(bin_sub)
  bin_sub$clade <- bin_sub$clade %>%
    gsub(".*?g__", "", .) %>%
    gsub("\\|s__.*", "", .)
  return(bin_sub)
}

fungal_avg <- function(fungal_tidy_df){
  library(tidyverse)
  avg_fungal <- fungal_tidy_df[2:5] %>%
    group_by(clade, phase, Study) %>% 
    summarise_all(mean) %>%
    as.data.frame()
  return(avg_fungal)
}

fungal_sum_for_stats <- function(fungal_tidy_df){
  library(tidyverse)
  fungal_sum <- fungal_tidy_df[1:4] %>%
    group_by(Sample, phase, Study) %>% 
    summarise_all(sum) %>%
    as.data.frame()
  filtered <- fungal_sum %>%
    filter(value > 0)
  final <- filtered[2:3]
  return(final)
}

fungal_prevalence_genera <- function(fungal_tidy_df){
  library(tidyverse)
  genus_prevalence <-  fungal_tidy_df[2:5] %>% 
    group_by(clade, phase, Study) %>% 
    summarise(Prevalence = (sum(value > 0)/n())*100)
  return(genus_prevalence)
}

fungal_prevalence_all <- function(fungal_tidy_df){
  library(tidyverse)
  sums <- fungal_tidy_df[1:4] %>%
    group_by(Sample, phase, Study) %>% 
    summarise_all(sum) %>%
    as.data.frame()
  prev <-  sums[2:4] %>% 
    group_by(phase, Study) %>% 
    summarise(Prevalence = (sum(value > 0)/n())*100)
  return(prev)
}

fungal_frequency_genera <- function(fungal_tidy_df){
  library(tidyverse)
  genus_frequency <-  fungal_tidy_df[2:5] %>% 
    group_by(clade, phase, Study) %>% 
    summarise(frequency = (sum(value > 0)))
  return(genus_frequency)
}


fungal_frequency_all <- function(fungal_tidy_df){
  library(tidyverse)
  sums <- fungal_tidy_df[1:4] %>%
    group_by(Sample, phase, Study) %>% 
    summarise_all(sum) %>%
    as.data.frame()
  frequency <-  sums[2:4] %>% 
    group_by(phase, Study) %>% 
    summarise(Frequency = (sum(value > 0)))
  return(frequency)
}

fungal_presence_all <- function(fungal_tidy_df){
  library(tidyverse)
  sums <- fungal_tidy_df[1:4] %>%
    group_by(Sample, phase, Study) %>% 
    summarise_all(sum) %>%
    as.data.frame()
  presence <- sums[2:4] %>%
    mutate(Presence = case_when(
      value > 0 ~ "Present",
      value == 0 ~ "Not_Present",
      TRUE ~ as.character(value)
    ))
  presence$freq <- 1
  frequency_sum <- presence[c(1:2,4:5)] %>%
    group_by(phase, Study, Presence) %>%
    summarise_all(sum) %>%
    as.data.frame
  return(frequency_sum)
}

presence_all <- function(fungal_tidy_df){
  library(tidyverse)
  sums <- fungal_tidy_df[1:4] %>%
    group_by(Sample, phase, Study) %>% 
    summarise_all(sum) %>%
    as.data.frame()
  presence <- sums[1:4] %>%
    mutate(Presence = case_when(
      value > 0 ~ "Present",
      value == 0 ~ "Not_Present",
      TRUE ~ as.character(value)
    ))
  return(presence)
}

## Statistics
print_model_summary <- function(phase, model) {
  cat("\n-----------------------------------------\n")
  cat("Phase:", phase, "\n")
  cat("-----------------------------------------\n")
  
  model_summary <- summary(model)
  print(model_summary)
  
  cat("\n-----------------------------------------\n\n\n")
}

extract_pvalues <- function(phase, model, df_sub, reference_group = "Dartmouth") {
  library(dplyr)
  model_summary <- summary(model)
  
  fixed_effects_names <- rownames(model_summary$coefficients)
  
  p_values <- model_summary$coefficients[, "Pr(>|z|)"]
  
  group2_names <- sub("Study", "", fixed_effects_names)
  
  pvalue_df <- data.frame(
    phase = rep(phase, length(group2_names)),
    group1 = reference_group,
    group2 = group2_names,
    p_value = p_values,
    signif = NA, 
    y.position = NA, 
    Study = reference_group
  )
  
  pvalue_df$signif <- ifelse(pvalue_df$p_value < .001, '***',
                             ifelse(pvalue_df$p_value < .01, '**',
                                    ifelse(pvalue_df$p_value < .05, '*',
                                           ifelse(pvalue_df$p_value < .1, '.', 'ns'))))
  
  pvalue_df$y.position <- ifelse(pvalue_df$group2 == "DIABIMMUNE", 26.5,
                                 ifelse(pvalue_df$group2 == "TEDDY", 28.0, NA))
  
  pvalue_df <- pvalue_df[!pvalue_df$group2 %in% c("(Intercept)", reference_group), ]
  
  return(pvalue_df)
}

pvalue_data <- list()

for (phase in c("Developmental", "Stable", "Transitional")) {
  df_sub <- filter(fig_5b, phase == phase)
  model <- glmer(Presence_Bool ~ Study + (1|Individual), data = df_sub, family = binomial, control = glmerControl(optimizer="bobyqa"))
  print_model_summary(phase, model)
  
  pvalue_df <- extract_pvalues(phase, model, df_sub)
  pvalue_data[[phase]] <- pvalue_df
}

pvalue_data <- do.call(rbind, pvalue_data)

## Plotting 
### Figure 5B
fungal_bar <- function(df_list, pvalue_df, cols){
  library(ggplot2)
  library(ggpubr)
  plot <- ggplot(df_list, aes(Study, Prevalence, fill = Study)) +
    geom_col(color = "gray20", alpha = 0.8, size = 0.5) +
    labs(y = "Fungal Prevalence") +
    theme_classic() +
    facet_wrap(~phase, ncol = 3, strip.position = "top") +
    scale_y_continuous(limits = c(0, 30)) +
    scale_color_manual(values = c("gray20", "gray20")) +
    scale_fill_manual(values = cols) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 8, color = "black"),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank())
  
  pvalue_plot <- plot + stat_pvalue_manual(
    data = pvalue_df,
    y = "y.position", 
    group = "phase",
    label = "signif",
    tip.length = 0.01
  )  
  
  return(pvalue_plot)
}

### Figure 5C
fungal_genus_bar <- function(df_list){
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(RColorBrewer)
  final_df <- rbindlist(df_list)  
  final_df <- final_df %>%
    filter(clade != "Cryptococcus")
  final_df$phase <- factor(final_df$phase, levels = c("Developmental", "Transitional", "Stable"))
  final_df <- na.omit(final_df)
  plot <- ggplot(final_df, aes(Study, Prevalence, fill = clade)) +
    geom_col(color = "gray20", alpha = 0.8, size = 0.5) +
    labs(y = "Fungal Prevalence") +
    theme_classic() +
    facet_wrap(~phase, ncol = 3, strip.position = "top") +
    scale_color_manual(values = c("gray20", "gray20")) +
    scale_fill_manual(values = brewer.pal(2, "Dark2"), name = "Genus") +
    theme(axis.text.x = element_text(angle = 90, size = 8, color = "black"), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 8, color = "black"),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank())
  return(plot)  
}

