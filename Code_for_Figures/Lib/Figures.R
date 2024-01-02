

improve_plot <- function() {
  #####
  # Add to any ggplot2 plot to make it look better
  #####

  library('ggplot2')
  t <- theme(axis.text = element_text(color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(linetype="dashed"))
  return(t)
}


variable_importance_overlap <- function(df_model_hayden, df_model_diabimmune, ntree=1000) {
  library('RRF')
  library('tidyverse')
  source("Lib/InfantAgeModels.R")
  df_list <- add_species_zeros(list(df_model_hayden, df_model_diabimmune))
  df_list <- add_species_zeros(list(df_model_healthy, df_model_diabimmune))
  
  df_model_hayden_full <- df_list[[1]]
  df_model_diabimmune_full <- df_list[[2]]
  
  model_hayden <- RRF(age ~ ., df_model_hayden_full, ntree = ntree)
  model_diabimmune <- RRF(age ~ ., df_model_diabimmune_full, ntree = ntree)
  df_importance_hayden <- get_importance_df(model_hayden)
  df_importance_diabimmune <- get_importance_df(model_diabimmune)
  df_importance_hayden <- rename(df_importance_hayden, IncNodePurity_ds1 = IncNodePurity)
  df_importance_diabimmune <- rename(df_importance_diabimmune, IncNodePurity_ds2 = IncNodePurity)
  
  df_importance <- full_join(df_importance_diabimmune, df_importance_hayden, on = "species") %>% mutate(IncNodePurity_ds1 = replace_na(IncNodePurity_ds1, 0), IncNodePurity_ds2 = replace_na(IncNodePurity_ds2, 0))
  return(df_importance)
}

simplify_name_standard <- function(sp) {
    # Simplifies to F. prausnitzii
    str_replace(sp, "s__([A-Z]{1})[a-z]+","\\1") %>% str_replace("_" ,". ")
}


importance_overlap_plot <- function(df_importance, ds1, ds2, limits.use = c(1e-3, 4e5), species_oi=NULL, outfile=NULL) {
  library('ggrepel')
  library('scales')
  library('tidyverse')
  library('extrafont')

  t <- improve_plot()
  xlab.use <- paste0(ds1, " species importance")
  ylab.use <- paste0(ds2, " species importance")

  cor.val <- cor(df_importance$IncNodePurity_ds1, df_importance$IncNodePurity_ds2, method = "spearman")
  df_importance$min <- apply(df_importance, 1, function(x) min(as.numeric(c(x["IncNodePurity_ds2"],x["IncNodePurity_ds1"]))))
  df_importance <- arrange(df_importance, desc(min))
  print("Spec oi")
  print(species_oi)
  # select a subset of important species
  if (!is.null(species_oi)) {
    df_importance_sub <- dplyr::filter(df_importance, species%in%species_oi)
  } else {
    df_importance_sub <- head(df_importance,5)
  }

  # Add the plotmath to make the species names italic
  df_importance_sub <- mutate(df_importance_sub, species=simplify_name_standard(species))
  df_importance_sub$species <- paste0("italic('", df_importance_sub$species, "')")
  

  p <- ggplot(df_importance, aes(x = IncNodePurity_ds1+1e-3, y = IncNodePurity_ds2+1e-3))+ geom_point(fill="steelblue", pch=21, size=3.5) +
    theme_bw(10)+ annotate("text", x=0.003, y=limits.use[2] * 0.8, size=8/.pt, label=sprintf("ρ=%.2f", cor.val)) +
    scale_x_continuous(labels = label_comma(), trans='log10', limits=limits.use) +
    scale_y_continuous(labels = label_comma(), trans='log10', limits=limits.use) +
    geom_label_repel(data = df_importance_sub, aes(x = IncNodePurity_ds1+1e-3, y=IncNodePurity_ds2+1e-3, label=species), force=100, nudge_y=0, nudge_x=-1,parse=TRUE, size=8/.pt) +
    t + ylab(ylab.use)+ xlab(xlab.use)+ theme(panel.grid.major.x = element_line(linetype="dashed"))
  
  if (!is.null(outfile)) {
    grDevices::cairo_pdf(outfile, width=3.5, height=3.5, family="Arial")
    print(p)
    dev.off()

    png(str_replace(outfile,".pdf",".png"), width=3.5, height=3.5, units="in", res=300, family="Arial")
    print(p)
    dev.off()

  }else {
    print(p)
  }
}


cross_dataset_predictions <- function(df_model_hayden, df_model_diabimmune, ds1, ds2, outfile, max.use=12, ntree=1000) {
  source("Lib/InfantAgeModels.R")
  library('gridExtra')
  library('extrafont')
  
  df_list <- add_species_zeros(list(df_model_hayden, df_model_diabimmune))
  df_model_hayden_full <- df_list[[1]]
  df_model_diabimmune_full <- df_list[[2]]
  
  r <- cross_dataset_model(df_model_hayden_full, df_model_diabimmune_full, ds1, ds2, method = "quantitative", ntree=ntree)
  
  t <- improve_plot()
  
  x.loc.denom <- 14

  p1 <- ggplot(r$df_plot.within.1, aes(x = actual, y = prediction)) + geom_point(pch=21, fill="steelblue") + ylab("Within-dataset predicted age") + xlab(paste0(ds1," age (months)")) +
    theme_bw(10) + t + ylim(c(0, max.use)) + xlim(c(1, max.use)) + annotate("text", x=max.use/x.loc.denom, y=max.use, size=8/.pt, label=sprintf("ρ=%.2f", r$df_results[1, 2]))
  p2 <- ggplot(r$df_plot.cross.1, aes(x = actual, y = prediction)) + geom_point(pch=21, fill="steelblue") + ylab("Cross-dataset predicted age") + xlab(paste0(ds1," age (months)")) +
    theme_bw(10) + t + ylim(c(0, max.use)) + xlim(c(1, max.use)) + annotate("text", x=max.use/x.loc.denom, y=max.use, size=8/.pt, label=sprintf("ρ=%.2f", r$df_results[3, 2]))
  p3 <- ggplot(r$df_plot.within.2, aes(x = actual, y = prediction)) + geom_point(pch=21, fill="steelblue") + ylab("Within-dataset predicted age") + xlab(paste0(ds2," age (months)")) +
    theme_bw(10) + t + ylim(c(0, max.use)) + xlim(c(1, max.use)) + annotate("text", x=max.use/x.loc.denom, y=max.use, size=8/.pt, label=sprintf("ρ=%.2f", r$df_results[2, 2]))
  p4 <- ggplot(r$df_plot.cross.2, aes(x = actual, y = prediction)) + geom_point(pch=21, fill="steelblue") + ylab("Cross-dataset predicted age") + xlab(paste0(ds2," age (months)")) +
    theme_bw(10) + t + ylim(c(0, max.use)) + xlim(c(1, max.use)) + annotate("text", x=max.use/x.loc.denom, y=max.use, size=8/.pt, label=sprintf("ρ=%.2f", r$df_results[4, 2]))
  
  grDevices::cairo_pdf(outfile, width=7.5, height=7.5, family="Arial")
  print(grid.arrange(p1, p2, p3, p4, ncol=2))
  dev.off()

  png(str_replace_all(outfile,".pdf",".png"), width=7.5, height=7.5, units="in", res=300, family="Arial")
  print(grid.arrange(p1, p2, p3, p4, ncol=2))
  dev.off()

}


venn_diagram_species_overlap <- function(df_importance, outfile, ds1='Hayden',ds2='Diabimmune') {
  library("VennDiagram")
  library('extrafont')
  n_include <- 50
  cex.use <- 1.0
  species.ds1 <- df_importance %>% arrange(desc(IncNodePurity_ds1)) %>% head(n_include) %>% pull(species)
  species.ds2 <- df_importance %>% arrange(desc(IncNodePurity_ds2)) %>% head(n_include) %>% pull(species)
  venn_data <- list(ds1=species.ds1, ds2=species.ds2)
  
  grid.newpage()
  p <- draw.pairwise.venn(area1 = n_include,
                          area2 = n_include,
                          cross.area = length(venn_data[[1]][venn_data[[1]]%in%venn_data[[2]]]), fill = c("blue", "red"), lwd=c(3, 3), category=c(ds1, ds2), cat.cex=c(cex.use, cex.use), cex=c(cex.use, cex.use, cex.use))
  
  for (i in 5:9) {
    p[[i]]$gp$fontfamily <- "Arial"
  }

  p[[8]]$vjust <- -3.25
  p[[8]]$hjust <- -0.1
  
  p[[9]]$vjust <- -3.25
  p[[9]]$hjust <- 1.25
  
  pdf(outfile, width=3.5, height=3.5, family="Arial")
  grid.draw(p)
  dev.off()

  png(str_replace(outfile,".pdf",".png"), width=3.5, height=3.5, units="in", res=300, family="Arial")
  grid.draw(p)
  dev.off()

}


binomial.se <- function(x) {
  p <- sum(x) / length(x)
  q <- 1- p
  return(sqrt(p * q / length(x)))
}
se <- function(x) {
  if (length(x) == 1)
    return(NA)
  return(sd(x)/sqrt(length(x)))
}


age_limits_hayden <- function(df_model_hayden_Healthy) {
  age.limits <- list()
  ages.test <- unique(df_model_hayden_Healthy$age)
  
  for (i in 1:length(ages.test)) {
    age.limits[[i]] <- c(ages.test[i]-1, ages.test[i]+1)
  }
  
  return(age.limits)
  
}

calc_prevalence_bins <- function(df_model1, species.oi, df_model2=NULL, age.limits=NULL, max_age=36, bin_size=3, multiple_corection="fdr") {
  library('tidyverse')
  if (is.null(df_model2))
    colnames.use <- c("species", "age_range", "n", "prevalence", "prevalence_se","abundance", "abundance_se", "nonzero_abundance", "nonzero_abundance_se")
  else
    colnames.use <- c("species", "age_range", "n_df1", "n_df2", "prevalence_df1","prevalence_df2", "prevalence_se_df1", "prevalence_se_df2", "y.text", "pval.wilcox","pval.fisher")
  
  if (is.null(age.limits)) {
    age.bins <- seq(from=0, to=max_age, by=bin_size)
    age.limits <- list()
    for (i in 1:(length(age.bins)-1)) {
      age.limits[[i]] <- c(age.bins[i], age.bins[i+1])
    }
  }
  df_out <- NULL
  for (i in 1:length(age.limits)) {
    df_model_sub1 <- dplyr::filter(df_model1, (age >= (age.limits[[i]][1])) & (age < (age.limits[[i]][2])))
    if (!is.null(df_model2))
      df_model_sub2 <- dplyr::filter(df_model2, (age >= (age.limits[[i]][1])) & (age < (age.limits[[i]][2])))
    for (species in species.oi) {
      vals1 <- df_model_sub1[, species]
      n1 <- nrow(df_model_sub1)
      d1 <- sum(vals1 > 0)
      p1 <- d1/n1
      p1.se <- binomial.se(vals1>0)
      
      if (!is.null(df_model2)) {
        vals2 <- df_model_sub2[, species]
        n2 <- nrow(df_model_sub2)
        d2 <- sum(vals2 > 0)
        p2 <- d2/n2
        p2.se <- binomial.se(vals2>0)
        
        if ((n1 < 3) | (n2 < 3)) {
          pval.wilcox <- NA
          pval.fisher <- NA
        }else {
          pval.wilcox <- wilcox.test(vals1, vals2)$p.value
          pval.fisher <- fisher.test(matrix(c(d1,(n1-d1), d2,(n2-d2)), ncol=2))$p.value
        }
        df_out <- bind_rows(df_out, data.frame(species, sprintf("%i-%i", age.limits[[i]][1],age.limits[[i]][2]), n1, n2, p1, p1.se, p2, p2.se, max(p1+p1.se, p2+p2.se)+0.02, pval.wilcox, pval.fisher))
      }else {
          df_out <- bind_rows(df_out, data.frame(species, sprintf("%i-%i", age.limits[[i]][1],age.limits[[i]][2]), n1, p1, p1.se, mean(vals1), se(vals1), mean(vals1[vals1>0]), se(vals1[vals1>0])))
      }
    }
  }
  colnames(df_out) <- colnames.use
  df_out$age <- df_out %>% pull(age_range) %>% sapply(function(x) mean(as.numeric(str_split(x, "-")[[1]])))
  if (!is.null(df_model2)) {
    df_out$pval.wilcox.mc <- p.adjust(df_out$pval.wilcox, multiple_corection)
    df_out$pval.fisher.mc <- p.adjust(df_out$pval.fisher, multiple_corection)
    df_out$label.text <- "*"
  }
  df_out$age_range <- factor(df_out$age_range, ordered=TRUE, levels=unique(df_out$age_range))
  
  return(df_out)
}


age_model_species_properties <- function(df_model_Healthy, species.oi, age.min=22.5, age.max=25.5) {
  library('tidyverse')
  # Looking properties of these 10 species
  df_plot <- data.frame()
  
  df_model_Healthy_sub <- df_model_Healthy %>% na.omit %>% filter((age > age.min) & (age < age.max)) %>% select(-Individual, -age)
  mean_abundance <- apply(df_model_Healthy_sub, 2, function(x) log10(mean(x[x > 0])))
  prevalence <- apply(df_model_Healthy_sub, 2, function(x) sum(x > 0) / length(x))
  df_plot <- bind_rows(df_plot, data.frame(log10abundance=mean_abundance, prevalence=prevalence, type="all species"))
  
  df_model_Healthy_sub <- df_model_Healthy %>% na.omit %>% filter((age > age.min) & (age < age.max))%>% select(-Individual, -age)
  mean_abundance <- apply(df_model_Healthy_sub[, species.oi], 2, function(x) log10(mean(x[x > 0])))
  prevalence <- apply(df_model_Healthy_sub[, species.oi], 2, function(x) sum(x > 0) / length(x))
  df_plot <- bind_rows(df_plot, data.frame(log10abundance=mean_abundance, prevalence=prevalence, type="10 species"))
  
  df_plot <- filter(df_plot, prevalence > 0)
  
  return(df_plot)
  
}

age_model_simplified_species <- function(df_model_Healthy, species_oi, species_oi_alt=NULL, ntree=1000, outfile=NULL) {
  library('gridExtra')
  library('extrafont')
  library('tidyverse')
  source("Lib/InfantAgeModels.R")

  predictions.1 <- predict_in_cv(df_model_Healthy, response_var="age", ntree=ntree)
  predictions.2 <- predict_in_cv(select(df_model_Healthy, age, one_of(species.oi)), response_var="age", ntree=ntree)

  df_prediction_plot <- data.frame(age=na.omit(df_model_Healthy)$age, predictions_full=predictions.1, predicitons_sub=predictions.2)
  
  cor.1 <- cor(df_prediction_plot$age, df_prediction_plot$predictions_full, method="spearman")
  cor.2 <- cor(df_prediction_plot$age, df_prediction_plot$predicitons_sub, method="spearman")
  
  p1 <- ggplot(df_prediction_plot, aes(x = age, y = predictions_full)) + geom_point(pch=21, fill="steelblue") + ylab("Predicted age") + xlab("Age (months)") + ggtitle("All species") + theme_bw(10) + t + annotate("text", x=4, y=38, size=8/.pt, label=sprintf("ρ=%.2f", cor.1))
  p2 <- ggplot(df_prediction_plot, aes(x = age, y = predicitons_sub)) + geom_point(pch=21, fill="steelblue") + ylab("Predicted age") + xlab("Age (months)") + ggtitle(paste0(length(species_oi)," best species")) + theme_bw(10) + t + annotate("text", x=4, y=38, size=8/.pt, label=sprintf("ρ=%.2f", cor.2))


  if (!is.null(species_oi_alt)) {
    predictions.3 <- predict_in_cv(select(df_model_Healthy, age, one_of(species.species_oi_alt)), response_var="age", ntree=ntree)
    df_prediction_plot <- data.frame(age=na.omit(df_model_Healthy)$age, predictions_full=predictions.1, predicitons_sub=predictions.2, predicitons_alt=predictions.3)
    cor.3 <- cor(df_prediction_plot$age, df_prediction_plot$predicitons_sub, method="spearman")
    p3 <- ggplot(df_prediction_plot, aes(x = age, y = predicitons_alt)) + geom_point(pch=21, fill="steelblue") + ylab("Predicted age") + xlab("Age (months)") + ggtitle("Alt best species") + theme_bw(12) + t + annotate("text", x=9, y=38, size=6, label=sprintf("ρ=%.2f", cor.3))
    grDevices::cairo_pdf(outfile, width=9, height=5.5, family="Arial")
    grid.arrange(p1, p2, p3, ncol=3)
    dev.off()

  } else {

    if (!is.null(outfile)) {
      grDevices::cairo_pdf(outfile, width=8, height=4.25, family="Arial")
      grid.arrange(p1, p2, ncol=2)
      dev.off()

      png(str_replace(outfile,".pdf",".png"), width=8, height=4.25, units="in", res=300, family="Arial")
      grid.arrange(p1, p2, ncol=2)
      dev.off()

    } else {
      grid.arrange(p1, p2, ncol=2)
    }
  }
}

choose_age_and_abundance <- function(df_prevalence, sp.oi, colonization_healthy, vals_healthy, use_static=TRUE, renorm_prevalence=FALSE, min.n=20) {
  if (use_static)
    return(list(age=colonization_healthy[[sp.oi]], abundance=vals_healthy[[sp.oi]]))
  
  df_prevalence_sub <- filter(df_prevalence, n>=min.n) %>% filter(species == sp.oi)
  
  if (renorm_prevalence) {
    df_prevalence_sub <- filter(df_prevalence_sub, prevalence <= 0.5)
  }
  # You need to normalize the prevalence for species that do not appear at 100% prevalence
  df_prevalence_sub$prevalence <- df_prevalence_sub$prevalence / max(df_prevalence_sub$prevalence)
  df_random <- df_prevalence_sub %>% filter(prevalence > runif(1))
  age.edit <- pull(df_random, age)[1]
  abundance.edit <- pull(df_random, nonzero_abundance)[1]
  return(list(age=age.edit, abundance=abundance.edit))
}

all_species_rescue <- function(species.oi, n_comb, df_model_CF, df_model_Healthy, indv, age.model.full=NULL, prevalence_bin_size=3, col_frac_use=0.25, ntree=10000, use_static=TRUE, renorm_prevalence=FALSE, summary.use="mean", max_do=500) {
  if (use_static)
    total_combs <- combn(species.oi, n_comb, simplify=FALSE)
  else {
    total_combs <- c()
    total.reps <- ceiling(max_do / length(combn(species.oi, n_comb, simplify=FALSE)))
    for (i in 1:total.reps) {
      total_combs <- c(total_combs, combn(species.oi, n_comb, simplify=FALSE))
    }
    total_combs <- sample(total_combs)
    stopifnot(length(total_combs) > max_do)
  }
  df_prevalence <- calc_prevalence_bins(df_model_Healthy, species.oi, bin_size=prevalence_bin_size)
  colonization_healthy <- identify_colonization_points(df_prevalence, species.oi, col_frac_use)
  ab_healthy <- calc_average_abundances(df_model_Healthy, species.oi)
  
  if (is.null(age.model.full))
    age.model.full <- relative_microbial_age(df_model_Healthy[, c(species.oi, "age")], df_model_CF[, c(species.oi, "age")], ntree=ntree)
  
  results <- list()
  i <- 1
  for (j in 1:length(total_combs)) {
    if ((j %% 100) == 0)
      print(sprintf("Currently on rep %i", j))
    species.test <- total_combs[[j]]
    df_model_CF_sub <- df_model_CF %>% filter(Individual == indv) %>% arrange(age)
    
    for (sp.oi in species.test) {
      vals <- choose_age_and_abundance(df_prevalence, sp.oi, colonization_healthy, vals_healthy, use_static, renorm_prevalence)
      df_model_CF_sub[which(df_model_CF_sub$age >= vals$age), sp.oi]<- vals$abundance
    }
    
    age.model.sub <- relative_microbial_age(df_model_Healthy[, c(species.oi, "age")], df_model_CF_sub[, c(species.oi, "age")], age.model.full$model_components, ntree=ntree)
    
    if (summary.use == "mean")
      summary.value <- age.model.sub$df_plot$relative_age %>% mean
    else
      summary.value <- age.model.sub$df_plot$relative_age %>% median
    if ((summary.value > 0) & (check_consecutive(age.model.sub$df_plot$relative_age, 2))) {
      p <- ggplot(age.model.sub$df_plot, aes(x=age, y=relative_age)) + geom_line() + geom_point(pch=21, size=3)+ theme_bw(16) +
        t + ylab("Relative age (months)") + xlab("Age (months)") + ggtitle(paste(sapply(species.test, simplify_name_standard), collapse="\n")) + theme(plot.title = element_text(size=12))
      results[[paste0("graph", i)]] <- p
      results[[paste0("species", i)]] <- species.test
      i <- i + 1
    }
    if (j >= max_do)
      break
  }
  return(results)
}


determine_color <- function(label1,label2, legend1, legend2) {
    color.use <- NULL
    if (is.na(label1))
        label1 <- "NA"
    if (is.na(label2))
        label2 <- "NA"

    if ((label1 == "*") & (label2 == "*")) {
        color.use <- "Both"
    } else if (label1 == "*") {
        color.use <- legend1
    } else if (label2 == "*") {
        color.use <- legend2
    }
    return(color.use)
}


prevalence_plot <- function(df_model_Healthy, df_model_CF, species.oi, hayden_bins=FALSE, outfile=NULL, df_model_Healthy_alt=NULL, healthy.alt.label="Healthy (Alt)", mode = "graph", min_n=10, title.use=NULL, healthy.label="Healthy",cf_label="CF", max_age_use=42, bin_size_use=3) {
  source("Lib/InfantAgeModels.R")
  library('extrafont')
  t <- improve_plot()
  
  df_model_CF <- add_missing_species(df_model_CF, species.oi)
  df_model_Healthy <- na.omit(df_model_Healthy)

  if (hayden_bins) {
    age.limits <- age_limits_hayden(df_model_CF)
    stopifnot(length(age.limits) < 25)
    df_prevalence_wt <- calc_prevalence_bins(df_model_Healthy, species.oi, age.limits=age.limits)
    df_prevalence_cf <- calc_prevalence_bins(df_model_CF, species.oi, age.limits=age.limits)
    df_prevalence_wt$type <- "Diabimmune"
    df_prevalence_cf$type <- "Hayden"
    # For the p-values
    df_prevalence_compare <- calc_prevalence_bins(df_model_Healthy, species.oi, df_model_CF, age.limits=age.limits)
  } else {
    df_prevalence_wt <- calc_prevalence_bins(df_model_Healthy, species.oi, max_age=max_age_use, bin_size=bin_size_use)
    df_prevalence_cf <- calc_prevalence_bins(df_model_CF, species.oi, max_age=max_age_use, bin_size=bin_size_use)
    df_prevalence_wt$type <- healthy.label
    df_prevalence_cf$type <- cf_label

    # For the p-values
    df_prevalence_compare <- calc_prevalence_bins(df_model_Healthy, species.oi, df_model_CF, max_age=max_age_use, bin_size=bin_size_use)

    if (!is.null(df_model_Healthy_Alt)) {
      df_prevalence_wt_alt <- calc_prevalence_bins(df_model_Healthy_alt, species.oi, max_age=max_age_use, bin_size=bin_size_use)
      df_prevalence_wt_alt$type <- healthy.alt.label
      df_prevalence_compare_alt <- calc_prevalence_bins(df_model_Healthy_alt, species.oi, df_model_CF, max_age=max_age_use, bin_size=bin_size_use)
      df_prevalence_wt <- bind_rows(df_prevalence_wt, df_prevalence_wt_alt)

    }
  }
  df_prevalence_wt$species <- str_replace_all(df_prevalence_wt$species, "s__","") %>% str_replace_all("_"," ")
  df_prevalence_cf$species <- str_replace_all(df_prevalence_cf$species, "s__","") %>% str_replace_all("_"," ")
  
  species.order.use <- str_replace_all(species.oi, "s__","") %>% str_replace_all("_"," ")
  df_prevalence <- bind_rows(df_prevalence_wt,df_prevalence_cf)
  df_prevalence$species <- factor(df_prevalence$species, ordered=TRUE, levels=species.order.use)
  df_prevalence <- filter(df_prevalence, n >= min_n)

  df_prevalence_compare$species <- str_replace_all(df_prevalence_compare$species, "s__","") %>% str_replace_all("_"," ")
  df_prevalence_compare$species <- factor(df_prevalence_compare$species, ordered=TRUE, levels=species.order.use)
  df_text <- dplyr::filter(df_prevalence_compare, pval.fisher.mc < 0.05)
  if (!is.null(df_model_Healthy_Alt)) {
    df_prevalence_compare_alt$species <- str_replace_all(df_prevalence_compare_alt$species, "s__","") %>% str_replace_all("_"," ")
    df_prevalence_compare_alt$species <- factor(df_prevalence_compare_alt$species, ordered=TRUE, levels=species.order.use)

    df_prevalence_compare_sub <- dplyr::filter(df_prevalence_compare, pval.fisher.mc < 0.05) %>% dplyr::select(species, age, y.text, label.text)
    df_prevalence_compare_alt_sub <- dplyr::filter(df_prevalence_compare_alt, pval.fisher.mc < 0.05) %>% dplyr::select(species, age, y.text, label.text)    
    df_text <- full_join(df_prevalence_compare_sub, df_prevalence_compare_alt_sub, by=c("species","age"))
    df_text$color <- apply(df_text,1, function(x) determine_color(x["label.text.x"],x["label.text.y"], healthy.label, healthy.alt.label))
    df_text$y.text <- apply(df_text,1, function(x) as.numeric(max(x["y.text.x"],x["y.text.y"],na.rm=TRUE)))

    df_text$label <- "*"
    df_text$species <- str_replace_all(df_text$species,"s__","") %>% str_replace_all("_"," ")                 
    df_text$species <- factor(df_text$species, ordered=TRUE, levels=species.order.use)
    
      save(df_prevalence_compare,df_prevalence_compare_alt,df_text, file = "prevalence_data.RData")

    p <- ggplot(df_prevalence, aes(x = age, y = prevalence)) + geom_line(lwd=0.75,alpha=0.75, aes(group=type, color=type)) + scale_color_manual(values = c("black","#1f78b4","#ff7f00", "#c33b3b")) + scale_y_continuous(labels = scales::percent) + xlab("Age (months)") + ylab("Prevalence") +
        geom_point(aes(color=type), size=0.75) + geom_errorbar(aes(ymin=prevalence-prevalence_se, ymax=prevalence+prevalence_se,group=type,color=type), width=1.25) + facet_wrap(species ~ ., ncol=4) + theme_bw(10) + t + theme(strip.text = element_text(size = 8, face = "italic"), legend.title=element_blank())
    p <- p + geom_text(data=df_text, aes(x=age, y=y.text,color=color,label=label), size=5.5)
    print(df_text$color)
  }
  else if (mode != "large") {
    p <- ggplot(df_prevalence, aes(x = age, y = prevalence)) + geom_line(lwd=0.75,alpha=0.75, aes(group=type, color=type)) + scale_color_manual(values = c("steelblue", "#c33b3b")) + scale_y_continuous(labels = scales::percent) + xlab("Age (months)") + ylab("Prevalence") +
      geom_point(aes(color=type), size=0.75) + geom_errorbar(aes(ymin=prevalence-prevalence_se, ymax=prevalence+prevalence_se,group=type,color=type), width=1.25) + facet_wrap(species ~ ., ncol=4) + theme_bw(10) + t + theme(strip.text = element_text(size = 8, face = "italic"), legend.title=element_blank())
    p <- p + geom_text(data=df_text, aes(x=age, y=y.text, label=label.text),size=5.5)
  } else {
    p <- ggplot(df_prevalence, aes(x = age, y = prevalence)) + geom_line(lwd=1.25,aes(group=type, color=type)) + scale_color_manual(values = c("steelblue", "#c33b3b")) + scale_y_continuous(labels = scales::percent) + xlab("Age (months)") + ylab("Prevalence") +
      geom_errorbar(aes(ymin=prevalence-prevalence_se, ymax=prevalence+prevalence_se,group=type,color=type), width=1.25) + facet_wrap(species ~ .) + theme_bw(16) + t + theme(strip.text = element_text(size = 10), legend.title=element_blank())
    p <- p + geom_text(data=df_text, aes(x=age, y=y.text, label=label.text),size=5.5)
    
  }
  
  if (!is.null(title.use))
    p <- p + ggtitle(title.use)
  
  if (!is.null(outfile)) {
    #pdf(outfile, width=8, height=3.5, family="Arial")
    pdf(outfile, width=9, height=4.75, family="Arial")
    print(p)
    dev.off()
  }
  else {
    print(p)
  }
}
