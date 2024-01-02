## Processing Metagenomics dataset to plot family line plot in Supplemental Figure 2A
make_family_df <- function(df, metadata, family_list, dataset = "Dartmouth", condition = "CF") {
  library(tidyverse)
  library(reshape2)
  ## Get Family data
  raw <- cbind(reshape::colsplit(df$clade, split = "\\|", names = c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')), df[2:ncol(df)])
  raw$family <- str_replace_all(raw$family, "f__", "")

  family_data <- raw[, c(5, 8:ncol(raw))]
  family_summaries <- family_data %>%
    group_by(family) %>%
    summarise_all(sum)
  family_t <- data.frame(t(family_summaries))
  names(family_t) <- as.matrix(family_t[1,])
  family_t <- family_t[-1,]
  family_t$Sample <- rownames(family_t)
  # if (dataset == "Dartmouth"){
  #   family_t$Sample <- sub("^", "D", family_t$Sample)
  # }
  family_subset <- subset(family_t, select = c(family_list, "Sample"))
  cols.num <- family_list
  family_subset[cols.num] <- lapply(cols.num, function(x) as.numeric(as.character(family_subset[[x]])))
  family_subset <- data.frame(family_subset[, family_list],
                              # sum_of_rows = apply(family_subset[, family_list], 1, sum) ,
                              Sample = family_subset$Sample)
  family_subset$Exp <- dataset
  family_subset$Disease_State <- condition
  ## Get required metadata
  meta <- metadata[,c("Sample", "Age")]
  meta_merge <- merge(family_subset, meta, by = "Sample")
  meta_merge$Age <- as.character(meta_merge$Age)
  final_df <- melt(meta_merge)
  final_df$Age <- as.numeric(as.character(final_df$Age))
  final_df <- na.omit(final_df)
  print(length(unique(final_df$Sample)))
  return(final_df)
}

family_plot <- function(combined_data, color_list, max_age = 36, ncol = 5){
  library(ggplot2)
  library(scales)
  plot <- ggplot(combined_data) + 
    geom_smooth(aes(Age, value, group = variable, color = variable), level = 0.95, method = "loess", alpha = 0.2) +
    facet_wrap(~Exp, ncol=ncol, scales = "free_x") +
    labs(x = "Age (Months)", y = "Relative Abundance") +
    theme_classic() +
    scale_x_continuous(limits = c(-0.001, max_age), breaks = pretty_breaks()) +
    scale_color_manual(name = "Family", values = color_list) +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size = 8, color = "black"),
          axis.text=element_text(size=8, color = "black"),
          axis.title=element_text(size=10, color = "black"),
          legend.text=element_text(size=8, color = "black"))
  
  return(plot)    
}

