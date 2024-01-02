## Processing Metagenomics dataset to plot Phylum plot in Figure 1A
make_phylum_df <- function(df, phyla_list, species_oi = "Escherichia_coli", dataset = "Dartmouth", condition = "CF") {
  library(tidyverse)
  library(reshape2)
  ## Get Phylum data
  raw <- cbind(reshape::colsplit(df$clade, split = "\\|", names = c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')), df[2:ncol(df)])
  raw$phylum <- str_replace_all(raw$phylum, "p__", "")
  raw$species <- str_replace_all(raw$species, "s__", "")
  
  phylum_data <- raw[, c(2, 8:ncol(raw))]
  phylum_summaries <- phylum_data %>%
    group_by(phylum) %>%
    summarise_all(sum)
  phylum_t <- data.frame(t(phylum_summaries))
  names(phylum_t) <- as.matrix(phylum_t[1,])
  phylum_t <- phylum_t[-1,]
  phylum_t$ID <- rownames(phylum_t)
  if (dataset == "Dartmouth"){
    phylum_t$ID <- sub("^", "D", phylum_t$ID)
  }
  phylum_subset <- subset(phylum_t, select = c(phyla_list, "ID"))
  cols.num <- phyla_list
  phylum_subset[cols.num] <- lapply(cols.num, function(x) as.numeric(as.character(phylum_subset[[x]])))
  phylum_subset <- data.frame(phylum_subset[, phyla_list],
                              sum_of_rows = apply(phylum_subset[, phyla_list], 1, sum) ,
                              ID = phylum_subset$ID )  # Divide each percentage by the sum of each row, to normalize the value to only the phyla we care about
  normalized_data <- sweep(phylum_subset[, cols.num], 1, phylum_subset$sum_of_rows, "/")
  normalized_data <- as.data.frame(lapply(normalized_data, function(x) x * 100))
  normalized_data$ID <- phylum_subset$ID
  ordered_data <- normalized_data[order(normalized_data$Proteobacteria, decreasing = TRUE),]
  ordered_data$Exp <- dataset
  ordered_data$Disease_State <- condition
  
  ## Get Species data
  species_data <- raw[, 7:ncol(raw)]
  species_data <- species_data[species_data$species == species_oi, ]
  species_t <- data.frame(t(species_data))
  names(species_t) <- as.matrix(species_t[1,])
  species_t$ID <- rownames(species_t)
  species_t$Exp <- dataset
  species_t$Disease_State <- condition
  species_t <- data.frame(species_t[-1,])
  if (dataset == "Dartmouth"){
    species_t$ID <- sub("^", "D", species_t$ID)
  }
  merged_data <- merge(species_t, ordered_data, by = c("Exp", "ID", "Disease_State"))
  final_order <- merged_data[order(merged_data$Proteobacteria, decreasing = TRUE),]
  return(final_order)
}

phyla_df_w_other <- function(df, phyla_list, species_oi = "Escherichia_coli", dataset = "Dartmouth", condition = "CF") {
  library(tidyverse)
  library(reshape2)
  ## Get Phylum data
  raw <- cbind(reshape::colsplit(df$clade, split = "\\|", names = c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')), df[2:ncol(df)])
  raw$phylum <- str_replace_all(raw$phylum, "p__", "")
  raw$species <- str_replace_all(raw$species, "s__", "")
  
  raw$phylum <- ifelse(raw$phylum %in% phyla_list, raw$phylum, 'Other')
  
  phylum_data <- raw[, c(2, 8:ncol(raw))]
  phylum_summaries <- phylum_data %>%
    group_by(phylum) %>%
    summarise_all(sum)
  phylum_t <- data.frame(t(phylum_summaries))
  names(phylum_t) <- as.matrix(phylum_t[1,])
  phylum_t <- phylum_t[-1,]
  phylum_t$ID <- rownames(phylum_t)
  if (dataset == "Dartmouth"){
    phylum_t$ID <- sub("^", "D", phylum_t$ID)
  }
  
  if (!"Other" %in% names(phylum_t)) {
    phylum_t$Other <- 0
    warning("No 'Other' phylum present. An 'Other' column has been added with zeros.")
  }
  
  phyla_list <- c("Other", phyla_list)
  
  phylum_subset <- subset(phylum_t, select = c(phyla_list, "ID"))
  cols.num <- phyla_list
  phylum_subset[cols.num] <- lapply(cols.num, function(x) as.numeric(as.character(phylum_subset[[x]])))
  phylum_subset <- data.frame(phylum_subset[, phyla_list],
                              sum_of_rows = apply(phylum_subset[, phyla_list], 1, sum) ,
                              ID = phylum_subset$ID )  # Divide each percentage by the sum of each row, to normalize the value to only the phyla we care about
  normalized_data <- sweep(phylum_subset[, cols.num], 1, phylum_subset$sum_of_rows, "/")
  normalized_data <- as.data.frame(lapply(normalized_data, function(x) x * 100))
  normalized_data$ID <- phylum_subset$ID
  ordered_data <- normalized_data[order(normalized_data$Proteobacteria, decreasing = TRUE),]
  ordered_data$Exp <- dataset
  ordered_data$Disease_State <- condition
  
  ## Get Species data
  species_data <- raw[, 7:ncol(raw)]
  species_data <- species_data[species_data$species == species_oi, ]
  species_t <- data.frame(t(species_data))
  names(species_t) <- as.matrix(species_t[1,])
  species_t$ID <- rownames(species_t)
  species_t$Exp <- dataset
  species_t$Disease_State <- condition
  species_t <- data.frame(species_t[-1,])
  if (dataset == "Dartmouth"){
    species_t$ID <- sub("^", "D", species_t$ID)
  }
  merged_data <- merge(species_t, ordered_data, by = c("Exp", "ID", "Disease_State"))
  final_order <- merged_data[order(merged_data$Proteobacteria, decreasing = TRUE),]
  return(final_order)
}


## Plotting function
# Takes combined df of all datasets of interest after running `make_phylum_df`
phylum_plot <- function(dfs, phyla_list, species_oi = "Escherichia_coli", datasets, color_list, ncol = 3) {
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  data <- dfs %>%
    select(all_of(c(phyla_list,
                    species_oi,
                    "Exp",
                    "ID",
                    "Disease_State")))
  data$Exp_Disease <- paste(data$Exp, data$Disease_State)
  data$Exp <- factor(data$Exp, levels = datasets)
  sorted <- data %>% arrange(Proteobacteria) %>%
    mutate(ID = factor(ID, levels = unique(ID[order(Proteobacteria, decreasing = TRUE)]))) %>% 
    gather(variable, value, -species_oi, -ID, -Exp, -Disease_State, -Exp_Disease, factor_key = TRUE) %>%
    na.omit()
  sorted[[species_oi]] <- as.numeric(as.character(sorted[[species_oi]]))
  # Make plot
  plot <- ggplot() + 
    geom_bar(data = sorted, position = "stack", stat = "identity", na.rm = TRUE, width = 1, aes(x = ID, y = value, fill = factor(variable, levels = phyla_list))) +
    geom_jitter(data = sorted, na.rm = TRUE, aes(x = ID, y = Escherichia_coli), color = "black", size = 0.25) +
    facet_wrap(~Exp, scales = "free_x", ncol = ncol) +
    scale_x_discrete(breaks = NULL, labels = NULL) +
    labs(y = "Relative Abundance", x = "Sample", size = 10, color = "black") +
    theme_minimal() +
    theme(strip.text.x = element_text(size = 8, color = "black"),
          axis.text=element_text(size=8, color = "black"),
          axis.title=element_text(size=10, color = "black"),
          legend.text=element_text(size=8, color = "black")) +
    scale_fill_discrete(name = "Phylum")
  
  plot <- plot + scale_fill_manual(values = c(color_list), name = "Phylum") 
  return(plot)
}
