## Format metadata
humann_meta <- function(metadata, control_clusters, control_study = "TEDDY", control_cond = "Non-CF",
                        cf_clusters, cf_study = "Dartmouth", cf_cond = "CF"){
  library(tidyverse)
  cluster_merge <- rbind(control_clusters, cf_clusters)
  colnames(cluster_merge)[1] <- "Sample"
  meta_merged <- metadata %>%
    filter(Cohort == control_study | Cohort == cf_study)
  meta_merged <- meta_merged %>%
    filter(Age >= 0) %>%
    filter(!is.na(Age))
  all_merged <- merge(cluster_merge, meta_merged, by = "Sample")
  binned <- all_merged %>%
    mutate(phase = case_when(
      Age <= 12 ~ "Developmental",
      Age > 12 & Age <= 18 ~ "Transitional",
      Age > 18 ~ "Stable",
      TRUE ~ as.character(Age)
    ))
  condition <- binned %>%
    mutate(Disease_State = case_when(
      Cohort == control_study ~ control_cond,
      Cohort == cf_study ~ cf_cond
    ))
  condition <- data.frame(condition)
  return(condition)
}

## Filter quant data for MaAsLin
humann_tidy <- function(control_df, cf_df, humann_meta_df, control_cond = "Non-CF"){
  library(tidyverse)
  rownames(humann_meta_df) <- humann_meta_df$Sample
  humann_meta_df$Disease_State_mod <- humann_meta_df$Disease_State
  control_cond <- as.character(control_cond)
  humann_meta_df$Disease_State_mod[humann_meta_df$Disease_State_mod == control_cond] = paste0("a_", control_cond)
  sample_list <- as.character(humann_meta_df$Sample)
  sample_list <- c("Gene", sample_list)
  quant <- merge(control_df, cf_df, by = c("# Gene Family"), all = TRUE)
  colnames(quant)[1] <- "Gene"
  quant_filter <- quant[, which((names(quant) %in% sample_list)==TRUE)]
  quant_filter[is.na(quant_filter)] = 0
  gene_list <- quant_filter$Gene[!grepl("s__", quant_filter$Gene) & !grepl("unclassified", quant_filter$Gene)]
  gene_filter <- quant_filter %>%
    filter(Gene %in% gene_list)
  KO_split <- separate(gene_filter, Gene, into = c("KO", "Gene_Family"), sep = ":", remove = FALSE, extra = "merge")
  output <- list(quant_filter = quant_filter, ko_split = KO_split)
  return(output)
}

## Get KO info from Humann Data
KOs_needed <- function(humann_tidy_df){
  library(tidyverse)
  KOs <- humann_tidy_df[2]
  return(KOs)
}

## Format for MaAsLin
split_kegg <- function(kegg_db){
  library(tidyverse)
  colnames(kegg_db) <- c("Module_Name",
                         "KO")
  kegg_db$Module_Name <- gsub("md:", "", kegg_db$Module_Name)
  kegg_db$KO <- gsub("ko:", "", kegg_db$KO)
  return(kegg_db)
}

humann_quant <- function(humann_tidy_df, split_kegg_df, humann_meta_df){
  quant_sub <- humann_tidy_df[,c(2,4:ncol(humann_tidy_df))]
  kegg_merge <- merge(split_kegg_df, quant_sub, by = "KO")
  melted <- reshape2::melt(kegg_merge, id.vars = c("KO",
                                                 "Module_Name"))
  colnames(melted)[3] <- "Sample"
  quant_meta <- merge(melted, humann_meta_df, by = "Sample")
  quant_meta$phase <- factor(quant_meta$phase, levels = c("Developmental",
                                                          "Transitional",
                                                          "Stable"))
  quant_meta <- na.omit(quant_meta)
  quant_meta_sub <- subset(quant_meta, select = c("Sample",
                                                  "Module_Name",
                                                  "phase",
                                                  "Disease_State",
                                                  "value"))
  return(quant_meta_sub)
}

maaslin_quant <- function(humann_quant_df){
  library(tidyverse)
  library(reshape2)
  quant_sums <- humann_quant_df %>%
    group_by(Module_Name, Sample, Disease_State, phase) %>%
    summarise_all(sum) %>%
    as.data.frame()
  maaslin_subset <- subset(quant_sums, select = c("Sample",
                                                  "Module_Name",
                                                  "value"))
  cast <- reshape2::dcast(maaslin_subset, 
                            Sample ~ Module_Name, 
                            value.var = "value")
  rownames(cast) <- cast$Sample
  maaslin_final <- cast[-1]
  maaslin_final <- data.frame(apply(maaslin_final, 2, 
                                    function(x) as.numeric(as.character(x))))
  rownames(maaslin_final) <- cast$Sample
  return(maaslin_final)
}

maaslin_meta <- function(humann_meta_df){
  library(tidyverse)
  humann_meta_df <- data.frame(humann_meta_df)
  rownames(humann_meta_df) <- humann_meta_df$Sample
  return(humann_meta_df)
}

run_maaslin <- function(maaslin_quant_df, maaslin_meta_df, 
                        outdir = "Dart_teddy_phases_maaslin_output_220518", 
                        fixed = c("Disease_State_mod", "phase"), 
                        abundance = 0.0001){
  library(Maaslin2)
  fit_data <- Maaslin2(
    maaslin_quant_df, maaslin_meta_df, outdir,
    fixed_effects = fixed,
    normalization = "NONE",
    plot_scatter = F,
    plot_heatmap = F,
    min_abundance = abundance)
  return(fit_data)
}

maaslin_means <- function(maaslin_results, humann_quant_df){
  library(tidyverse)
  filtered <- maaslin_results %>%
    filter(qval <= 0.05)
  unq <- data.frame(unique(maaslin_results$feature))
  colnames(unq)[1] <- "Module_Name"
  merge_quant <- merge(unq, humann_quant_df, by = "Module_Name")
  sub <- subset(merge_quant, select = c("Module_Name",
                                        "phase",
                                        "Disease_State",
                                        "value"))
  means <- sub %>%
    group_by(Module_Name, Disease_State, phase) %>%
    summarise_all(mean)
  means <- na.omit(means)
  means$col_name <- paste(means$Disease_State, 
                          means$phase,
                          sep = ": ")
  means$col_name <- factor(means$col_name,
                           levels = c("CF: Developmental",
                                      "CF: Transitional",
                                      "CF: Stable",
                                      "Non-CF: Developmental",
                                      "Non-CF: Transitional",
                                      "Non-CF: Stable"))
  
  return(means)
}

heatmap_tidy <- function(maaslin_means_df){
  library(reshape2)
  cast <- reshape2::dcast(maaslin_means_df, 
                          Module_Name ~ col_name, 
                          value.var = "value")
  cast[is.na(cast)] = 0
  return(cast)
}

split_paths <- function(kegg_paths_df){
  library(tidyverse)
  path_sp <- separate(kegg_paths_df, 
                      Module_Name, 
                      into = c("Module_Name",
                               "Module_Description"),
                      sep = 6, 
                      remove = FALSE)
  
  path_sp$Module_Description <- str_replace(path_sp$Module_Description, 
                                            "\\s*\\(\\d+\\)\\s*\\(complete \\d+/\\d+\\)$", 
                                            "")
  path_sp$Module_Description <- str_replace(path_sp$Module_Description, 
                                            "\\s*\\(\\d+\\)\\s*\\(\\d+\\sblocks?\\smissing\\s\\d+\\/\\d+\\)$", 
                                            "")
  return(path_sp)
}
  
heatmap_data <- function(split_paths_df, heatmap_data_df){
  library(tidyverse)
  library(scales)
  merged <- merge(split_paths_df, heatmap_data_df, by = "Module_Name")
  final <- merged[,3:9]
  rownames(final) <- final$Module_Description
  mat <- as.matrix(final[,c("CF: Developmental",
                            "CF: Transitional",
                            "CF: Stable",
                            "Non-CF: Developmental",
                            "Non-CF: Transitional",
                            "Non-CF: Stable")])
  mat <- t(apply(mat, 1L, scales::rescale))
  return(mat)
}

pathway_annot <- function(split_paths_df){
  library(tidyverse)
  library(data.table)
  path_sub <- subset(split_paths_df, select = c("Module_Description",
                                                "Overall_Pathway"))
  unq <- unique(path_sub[, c("Module_Description",
                             "Overall_Pathway")])
  rownames(unq) <- unq$Module_Description
  annot <- data.table(unq[2])
  rownames(annot) <- unq$Module_Description
  return(annot)
}

group_annot <- function(maaslin_means_df){
  library(tidyverse)
  group_sub <- subset(maaslin_means_df, select = c("col_name",
                                                   "Disease_State"))
  unq <- unique(group_sub[, c("col_name",
                           "Disease_State")])
  groups <- data.frame(unq[2])
  rownames(groups) <- unq$col_name
  colnames(groups)[1] <- "Disease_State"
  return(groups)
}

plot_heatmap <- function(heatmap_data_df, 
                         pathway_annot_df, 
                         group_annot_df,
                         num_rows = 3,
                         num_cols = 3,
                         font_size = 5){
  library(pheatmap)
  library(viridis)
  plot <- pheatmap(heatmap_data_df,
                   annotation_row = pathway_annot_df,
                   annotation_col = group_annot_df,
                   cluster_cols = T,
                   color = viridis(25),
                   border_color = "gray50",
                   clustering_distance_rows = "correlation",
                   cellwidth = 5,
                   cellheight = 5,
                   gaps_col = c(3),
                   cutree_rows = c(num_rows),
                   cutree_cols = c(num_cols),
                   labels_col = c("CF: Developmental (0-12)",
                                  "CF: Trantitional (12-18)",
                                  "CF: Stable (18-36)",
                                  "Non-CF: Developmental (0-12)",
                                  "Non-CF: Transitional (12-18)",
                                  "Non-CF: Stable (18-36)"),
                   angle_col = 90,
                   fontsize = font_size)
  return(plot)
} 

save_heatmap_svg <- function(heatmap, 
                             file = "~/Documents/Ross_Lab_Thesis_Work/MiGS_Metagenomics_210308/Manuscript/CF_Metagenomics_Manuscript/ps_figure_svgs/Dart_TEDDY_phases_heatmap_230325.svg",
                             height = 11,
                             width = 7,
                             dpi = 1500){
  library(ggplot2)
  ggsave(file, heatmap, height = height, width = width, dpi = dpi)
}

## Species based figure functions
sig_species <- function(humann_tidy, 
                        split_kegg_df, 
                        split_paths, 
                        maaslin_means){
  library(tidyverse)
  species_list <- humann_tidy$Gene[grepl("s__", humann_tidy$Gene)]
  genes <- humann_tidy %>%
    filter(Gene %in% species_list)
  split_genes <- separate(genes, 
                    Gene, 
                    into = c("KO", "Gene_Family"), 
                    sep = ": ",
                    remove = F, 
                    extra = "merge")
  split_genes2 <- separate(split_genes, 
                           Gene_Family, 
                           into = c("Gene_Family", "Species"),
                           sep = ".s__",
                           remove = T,
                           extra = "merge")
  mod_merge <- merge(split_kegg_df, split_genes2, by = "KO")
  path_merge <- merge(split_paths, mod_merge, by = "Module_Name")
  sig <- maaslin_means[1] 
  sig_mod_merge <- merge(path_merge, sig, by = "Module_Name")
  unq <- unique(sig_mod_merge)
  path_sub <- subset(path_merge, select = c(2, 7:ncol(path_merge)))
  return(path_sub)
}

melt_sig_species <- function(sig_species_df, humann_meta_df){
  library(tidyverse)
  library(reshape2)
  melted <- melt(sig_species_df, id.vars = c("Overall_Pathway",
                                          "Species"))
  colnames(melted)[3] <- "Sample"
  final <- merge(melted, humann_meta_df, by = "Sample")
  return(final)
}

sig_species_relab <- function(all_metaphlan, sig_species_df, humann_meta_df){
  library(tidyverse)
  library(reshape2)
  species <- sig_species_df$Species
  all_metaphlan$clade <- gsub(".*s__", "", all_metaphlan$clade)
  sp_metaphlan <- all_metaphlan %>%
    filter(clade %in% species)
  sp_metaphlan[is.na(sp_metaphlan)] <- 0
  melted <- melt(sp_metaphlan)
  colnames(melted) <- c("Species",
                        "Sample",
                        "relab")
  final <- merge(melted, humann_meta_df, by = "Sample")
  return(final)
}

# Stats for top species 
wilcoxon <- function(melted_df1, melted_df2){
  library(tidyverse)
  final_df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), 
                       c("Species", "phase", "p_value"))
  for (species in unique(melted_df1$Species)){
    for (Phase in unique(melted_df1$phase)){
      v1 <- filter(melted_df1, (Species==species) & (phase==Phase)) %>%
        pull(relab)
      v2 <- filter(melted_df2, (Species==species) & (phase==Phase)) %>% 
        pull(relab)
      wilcox <- print(wilcox.test(v1,v2))
      res <- wilcox$p.value
      new_row <- data.frame(Species = c(species), 
                            phase = c(Phase), 
                            p_value = c(res))
      final_df <- rbind(final_df, new_row)
      
    }
  }
  final_df$p.adjust <- p.adjust(final_df$p_value, method = "bonferroni", n = length(final_df$p_value))  
  return(final_df)
}

calc_fc <- function(melted_df1, 
                    type_df1 = "Non-CF", 
                    melted_df2, 
                    type_df2 = "CF"){
  library(tidyverse)
  library(gtools)
  melted_df1$Disease_State <- type_df1
  melted_df2$Disease_State <- type_df2
  tmp <- rbind(melted_df1, melted_df2)
  log_sum <- tmp %>%
    group_by(Disease_State, phase, Species) %>%
    summarize(mean = mean(relab))
  
  log_calc <- log_sum %>%
    group_by(phase, Species) %>%
    summarize(FC = foldchange(mean[Disease_State == type_df2],
                              mean[Disease_State == type_df1]),
              Log2FC = foldchange2logratio(FC))
  return(log_calc)
} 

all_stats <- function(wilcoxon_res, fc_res){
  library(tidyverse)
  stats <- merge(wilcoxon_res, fc_res, by = c("phase",
                                              "Species"))
  return(stats)
}

# For volcano plots

volcano_data <- function(all_stats_df, alpha = 0.05, fc_cutoff = 1, 
                         CF_col = "#1F78B4", 
                         control_col = "#FF7F00"){
  tmp <- as.data.frame(all_stats_df)
  tmp$cols <- c()
  tmp$cols <- NA
  for(i in 1:nrow(tmp)){
    if(is.na(tmp$p_value[i])){
      tmp$cols[i] <- NA
    }
    else if(tmp$p_value[i] <= alpha & tmp$Log2FC[i] > fc_cutoff){
      tmp$cols[i] <- CF_col
    } 
    else if(tmp$p_value[i] <= alpha & tmp$Log2FC[i] < -fc_cutoff){
      tmp$cols[i] <- control_col
    } 
    else if(tmp$p_value[i]<=alpha & tmp$Log2FC[i]>-fc_cutoff & tmp$Log2FC[i]<fc_cutoff){
      tmp$cols[i] <- "gray47"
    } 
    else if(tmp$p_value[i]>alpha & tmp$Log2FC[i] > fc_cutoff){
      tmp$cols[i] <- "gray47" 
    }
    else if(tmp$p_value[i]>alpha & tmp$Log2FC[i] < -fc_cutoff){
      tmp$cols[i] <- "gray47" 
    }
    else if(tmp$p_value[i]>alpha & tmp$Log2FC[i] < fc_cutoff){
      tmp$cols[i] <- "gray10" 
    }
  }
  tmp_filter <- na.omit(tmp)
  tmp_filter <- tmp_filter[!(tmp_filter$Log2FC == "-Inf" | tmp_filter$Log2FC == "Inf"),]
  tmp_filter$phase <- factor(tmp_filter$phase, levels = c("Developmental",
                                                          "Transitional",
                                                          "Stable"))
  tmp_filter$Species_edit <- paste(substr(tmp_filter$Species, 1, 1), 
                                   gsub(".*_", "", tmp_filter$Species), 
                                   sep = ". ")
  return(tmp_filter)
}

top_species <- function(volcano_data_df, 
                         n_species = 5, 
                         CF_col = "#1F78B4", 
                         control_col = "#FF7F00",
                         treat_lab = "CF",
                         control_lab = "Non-CF",
                         fc_cutoff = 1){
  library(tidyverse)
  library(ggplot2)
  high_sig <- volcano_data_df %>%
    filter(cols == CF_col) %>%
    filter(Log2FC > 0)
  pval_greater <- high_sig[order(high_sig$p_value),]
  treat_pval <- Reduce(rbind,
                    by(pval_greater,
                       pval_greater["phase"],
                       head,
                       n = n_species))
  treat_pval$Disease_State <- treat_lab
  low_sig <- volcano_data_df %>%
    filter(cols == control_col) %>%
    filter(Log2FC < 0)
  pval_lower <- low_sig[order(low_sig$p_value),]
  control_pval <- Reduce(rbind,
                         by(pval_lower,
                            pval_lower["phase"],
                            head,
                            n = n_species))
  control_pval$Disease_State = control_lab
  top_species <- rbind(treat_pval, control_pval)
  return(top_species)
}

plot_volcano <-  function(volcano_data_df, top_species_df,
                          alpha = 0.05, fc_cutoff = 1){
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  phases <- list("Developmental", "Transitional", "Stable")
  plot <- ggplot(volcano_data_df, 
                 aes(Log2FC, 
                     -log10(p_value))) +
    geom_point(aes(col=col), 
               alpha = 0.5, 
               size = 1.5, 
               colour = volcano_data_df$cols, 
               fill = volcano_data_df$cols) +
    xlab("Log2 Fold Change") + ylab("-Log10 P-value") +
    ylim(0, 50) +
    xlim(-21, 21) +
    geom_hline(yintercept = -log10(alpha),
               color = "black",
               linetype = "dashed",
               size = 0.4) +
    theme(legend.key = element_blank()) +
    theme_classic() +
    facet_wrap(~phase, ncol = 1) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 10,
                                    color = "black"),
          axis.title = element_text(size = 10,
                                    color = "black"),
          axis.text.x = element_text(size = 8,
                                     color = "black"))
  annot_plot <- plot +
    geom_label_repel(data = subset(data.frame(top_species_df),
                                   phase == "Developmental"),
                                   aes(label = Species_edit),
                                   box.padding   = 0.25,
                                   nudge_x = -0.25,
                                   nudge_y = 0.25,
                                   label.size = 0.1,
                                   segment.size = 0.3,
                                   segment.color = 'grey50',
                                   size = 1.5,
                                   max.overlaps = Inf,
                                   min.segment.length = 0) +
    geom_label_repel(data = subset(data.frame(top_species_df),
                                   phase == "Transitional"),
                     aes(label = Species_edit),
                     box.padding   = 0.25,
                     nudge_x = -0.25,
                     nudge_y = 0.25,
                     label.size = 0.1,
                     segment.size = 0.3,
                     segment.color = 'grey50',
                     size = 1.5,
                     max.overlaps = Inf,
                     min.segment.length = 0) +
    geom_label_repel(data = subset(top_species_df,
                                   phase == "Stable"),
                     aes(label = Species_edit),
                     box.padding   = 0.25,
                     nudge_x = -0.25,
                     nudge_y = 0.25,
                     label.size = 0.1,
                     segment.size = 0.3,
                     segment.color = 'grey50',
                     size = 1.5,
                     max.overlaps = Inf,
                     min.segment.length = 0) +
    geom_vline(xintercept = fc_cutoff,
               color = "black",
               linetype = "dotted") +
    geom_vline(xintercept = -fc_cutoff,
               color = "black",
               linetype = "dotted")
  return(annot_plot)
}

top_species_unique <- function(top_species_df){
  library(tidyverse)
  species_sub <- subset(top_species_df, select = c("Species",
                                                   "phase",
                                                   "Disease_State"))
  unq <- unique(species_sub[c("Species",
                              "phase",
                              "Disease_State")])
  return(unq)
}

# For the top species boxplots
top_species_metaphlan <- function(metaphlan_df, metadata_df, pval_df){
  library(tidyverse)
  library(reshape2)
  metaphlan_df$clade <- gsub(".*s__", "", metaphlan_df$clade)
  melted <- melt(metaphlan_df)
  colnames(melted) <- c("Species", "Sample", "relab")
  merged <- merge(melted, metadata_df[c("Sample", 
                                        "phase", 
                                        "Disease_State")], 
                  by = "Sample")
  filtered <- merge(merged, pval_df, by = c("Species",
                                            "phase",
                                            "Disease_State"))
  species_sub <- subset(filtered, select = c("Species",
                                             "phase",
                                             "Disease_State",
                                             "relab"))
  species_sub <- na.omit(species_sub)
  prevalence <- species_sub %>%
    group_by(Species, phase, Disease_State) %>%
    summarise(Prevalence = (sum(relab > 0)/n())*100)
  final <- merge(filtered, prevalence, by = c("Species",
                                              "phase",
                                              "Disease_State"))
  return(final)
}

pval_boxplot_data <- function(top_metaphlan_df){
  library(tidyverse)
  filtered <- top_metaphlan_df %>%
    filter(relab > 0)
  filtered$phase <- factor(filtered$phase, levels = c("Developmental",
                                                      "Transitional",
                                                      "Stable"))
  filtered$Disease_State <- factor(filtered$Disease_State, levels = c("Non-CF",
                                                                      "CF"))
  filtered$Species_edit <- paste(substr(filtered$Species, 1, 1),
                                 gsub(".*_", "", filtered$Species), sep = ". ")
  final <- filtered %>%
    mutate(annot = paste0(Species_edit, " (", round(Prevalence, digits = 2), ")"))
  return(final)
}

boxplot_plot <- function(boxplot_data, color_list = c("#FF7F00", "#1F78B4")){
  library(ggplot2)
  plot <- ggplot(boxplot_data, aes(x = relab, 
                                   y = reorder(annot, relab),
                                   fill = Disease_State)) +
    geom_jitter(color = "gray20", 
                alpha = 0.5, 
                size = 0.1) +
    geom_boxplot(notch = F, 
                 varwidth = T, 
                 color = "gray20", 
                 outlier.shape = NA, 
                 alpha = 0.8, 
                 size = 0.5) +
    stat_boxplot(geom = "errorbar",
                 width = 0.6,
                 linetype = 1,
                 color = "gray20") +
    labs(x = "Relative Abundance",
         size = 10,
         color = "black") +
    theme_classic() +
    facet_wrap(phase~Disease_State,
               ncol = 2,
               strip.position = "top",
               scales = "free_y") +
    scale_color_manual(values = c("gray20", "gray20")) +
    scale_fill_manual(values = color_list) +
    theme(axis.text = element_text(size = 8,
                                   color = "black"),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          legend.position = "none",
          strip.placement = "outside")
  return(plot)
}

## For the top species bar charts
top_species_humann <- function(top_species_df, sig_species_df){
  library(tidyverse)
  all <- merge(top_species_df, sig_species_df)
  unique(all[c("Species",
               "phase",
               "Disease_State")])
  all_sub <- subset(all, select = c("Species",
                                    "phase",
                                    "Disease_State",
                                    "Overall_Pathway",
                                    "value"))
  all_sum <- all_sub %>%
    filter(value > 0) %>%
    group_by(Species, phase, Disease_State, Overall_Pathway) %>%
    summarise_all(mean) %>%
    as.data.frame()
  sums <- all_sum[,c(1:3,5)] %>%
    group_by(Species, phase, Disease_State) %>%
    summarise_all(sum) %>%
    as.data.frame()
  colnames(sums)[4] <- "totals"
  final <- merge(all_sum, sums, by = c("Species",
                                       "phase",
                                       "Disease_State"))
  final <- final[order(final$totals),]
  final$Species_edit <- paste(substr(final$Species, 1, 1), gsub(".*_", "", final$Species), sep = ". ")
  final$Disease_State <- factor(final$Disease_State, levels = c("Non-CF",
                                                                "CF"))
  final$phase <- factor(final$phase, levels = c("Developmental",
                                                "Transitional",
                                                "Stable"))
  return(final)
}

humann_bar_plot <- function(top_species_humann_df){
  library(ggplot2)
  library(RColorBrewer)
  plot <- ggplot(top_species_humann_df, aes(value,
                                            reorder(Species_edit,
                                                    value,
                                                    sum),
                                            fill = Overall_Pathway)) +
    geom_bar(position = "stack", 
             stat = "identity", 
             color = "black",
             fun.y = sum,
             size = 0.1) +
    theme_classic() +
    labs(x = "Pathway Abundance") +
    facet_wrap(phase ~ Disease_State, 
               ncol = 2, 
               scales = "free_y") +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 10,
                                    color = "black"),
          axis.text.y = element_text(size = 8,
                                     color = "black"),
          axis.text.x = element_text(size = 8,
                                     angle = 45,
                                     color = "black",
                                     hjust = 1),
          axis.title.y = element_blank()) +
    scale_x_continuous(labels = function(x) format(x, scientific = T)) +
    scale_fill_manual(values = brewer.pal(10, "RdYlBu"), 
                      name = "Pathway")
  return(plot)
}

## Mean abundance plots
relab_meta_merge <- function(metaphlan_df, metadata_df){
  library(tidyverse)
  library(reshape2)
  melted <- melt(metaphlan_df, id.vars = "clade")
  colnames(melted)[2] <- "Sample"
  merged <- merge(melted, metadata_df)
  merged$clade <- gsub(".*?s__", "", merged$clade)
  return(merged)
}

calc_means <- function(relab_meta_df){
  library(tidyverse)
  sub <- subset(relab_meta_df, select = c("clade",
                                          "phase",
                                          "Disease_State",
                                          "value"))
  means <- sub %>%
    filter(value > 0) %>%
    group_by(clade, phase, Disease_State) %>%
    summarise_all(mean) %>%
    as.data.frame()
  ordered <- means[order(-means$value),]
  return(ordered)
}

get_top_species <- function(calc_means_df, n_species = 5, control = "Non-CF", treated = "CF"){
  library(tidyverse)
  top <- Reduce(rbind,
                by(calc_means_df,
                   calc_means_df[c("phase", "Disease_State")],
                   head,
                   n = n_species))
  colnames(top)[1] <- "Species"
  return(top)
}
