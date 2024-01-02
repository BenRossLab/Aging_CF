subset_to_species <- function(df_model, species_oi, min_ab=0.1) {
    library('tidyverse')
    for (sp in species_oi) {
        if (!sp%in%colnames(df_model))
            df_model[,sp] <- 0
    }
    df_species <- select(df_model, one_of(species_oi))

    df_species_relabl <- apply(df_species, 1, function(x) x/sum(x)) %>% t
    df_species <- df_species[apply(df_species_relabl,1,function(x) sum(x) > min_ab),]
    return(na.omit(df_species))
}


convert_relab_to_counts <- function(df_species) {
    # Very hacky
    return(apply(df_species, 2, function(x) round(x * 10000)))
}


try_n_clusters <- function(df_species, clusters.test=1:15, verbose=TRUE) {
    library('DirichletMultinomial')
    models <- lapply(clusters.test, dmn, count=data.matrix(df_species), verbose=verbose)
    names(models) <- as.character(clusters.test)
    return(models)
}


identify_best_n_clusters <- function(models, method="laplace", clusters.test=1:15) {
    library('DirichletMultinomial')
    library('ggplot2')
    
    t <- theme(axis.text = element_text(color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(linetype="dashed"))

    if (method == "AIC")
        scores <- sapply(models, AIC)
    else if (method == "BIC")
        scores <- sapply(models, BIC)
    else if (method == "laplace")
        scores <- sapply(models, laplace)
    names(scores) <- clusters.test
    df_plot <- data.frame(n_test=clusters.test, score=scores)
    p <- ggplot(df_plot, aes(x=n_test, y=score)) + geom_line() + geom_point(pch=21, size=3, fill="steelblue") + theme_bw(16) + t
    
    n.clusters <- names(scores)[which.min(unlist(scores))]
    best <- models[[as.character(n.clusters)]]

    return(list(scores=scores, df_plot=df_plot, plot=p, model.best=best, n.clusters=which.min(unlist(scores))))
}


order_clusters <- function(df_plot) {
    return(df_plot %>% apply(2, function(x) which(x == max(x))[1]) %>% order)
}

assign_samples <- function(df_species, model, predict=FALSE) {

    # Ensure your columns are not factors
    for (col in colnames(df_species)) {
        stopifnot(!is.factor(df_species[,col]))
    }
    if (predict)
        predictions <- predict(model, newdata=data.matrix(df_species))
    else
        predictions <- mixture(model)
    assignments <- apply(predictions, 1, which.max)
    
    names(assignments) <- rownames(df_species)
    return(assignments)
}


calc_samples_per_cluster <- function(df_model, ass, n.clusters, min_datapts=25, age.range = seq(from=0,to=36,by=3)) {
    library('tidyverse')
    df_out <- data.frame()
    df_model$sample <- rownames(df_model)
    df_model <- filter(df_model, sample%in%names(ass))
    for (i in 1:(length(age.range)-1)) {
        df_model_sub <- dplyr::filter(df_model, (age >= age.range[i]) & (age < age.range[i+1]))
        if (nrow(df_model_sub) < min_datapts) {
            print(sprintf("skipping %i-%i",age.range[i],age.range[i+1]))
            next
        }
        fraction.clusters <- sapply(1:n.clusters, function(j) sum(df_model_sub$sample%in%names(ass)[ass==j])/nrow(df_model_sub))
        df_out <- bind_rows(df_out, data.frame(bin=sprintf("%i-%i",age.range[i],age.range[i+1]),cluster=1:n.clusters,fraction.clusters))
    }

    # Add on all pts > max
    i <- length(age.range)
    df_model_sub <- dplyr::filter(df_model, (age >= age.range[i]))
    fraction.clusters <- sapply(1:n.clusters, function(j) sum(df_model_sub$sample%in%names(ass)[ass==j])/nrow(df_model_sub))
    df_out <- bind_rows(df_out, data.frame(bin=sprintf(">%i",age.range[i]),cluster=1:n.clusters,fraction.clusters))


    df_out_plot <- df_out %>% pivot_wider(names_from=cluster, values_from=fraction.clusters) 
    df_out_plot <- data.frame(df_out_plot)
    rownames(df_out_plot) <- df_out_plot$bin
    df_out_plot <- df_out_plot[,2:ncol(df_out_plot)]
    colnames(df_out_plot) <- 1:n.clusters
    
    return(df_out_plot)
}



# Make a point for every cell greater than the min_plot
make_graph <- function(df_plot, min_plot=0.04, outfile=NULL, tick.names=NULL, asp.use=NULL, vertex.scale.use=NULL, df_plot_nrow_override=NULL) {
    library('igraph')
    library('extrafont')
    n.clusters <- ncol(df_plot)
    cols <- inlmisc::GetColors(n = n.clusters)

    layout <- data.frame()
    cols.full <- c()
    k <- 1
    graph.vector <- c()
    vertex_weights <- c()
    for (i in 1:nrow(df_plot)) {
        for (j in 1:ncol(df_plot)) {
            if (df_plot[i,j] < min_plot)
                next
            layout <- bind_rows(layout, data.frame(i=i,j=n.clusters-j+1))
            cols.full <- c(cols.full,cols[j])
            graph.vector <- c(graph.vector, c(k,k))
            vertex_weights <- c(vertex_weights,df_plot[i,j])
            k <- k + 1
        }
    }
    
    # This puts a tiny circle at the end.
    if (!1%in%unique(layout$j)) {
        layout <- bind_rows(layout, data.frame(i=1,j=1))
        graph.vector <- c(graph.vector, c(k,k))
        vertex_weights <- c(vertex_weights, 0)
        cols.full <- c(cols.full,"#00FFFFFF")
        k <- k + 1
    }

    # This puts a tiny circle at the beginning
    
    if (!ncol(df_plot)%in%unique(layout$j)) {
      layout <- bind_rows(layout, data.frame(i=1,j=ncol(df_plot)))
      graph.vector <- c(graph.vector, c(k,k))
      vertex_weights <- c(vertex_weights, 0)
      cols.full <- c(cols.full,"#00FFFFFF")
      
    }
    
    if (is.null(tick.names)) {
        tick.names <- rownames(df_plot)

    }

    if (is.null(df_plot_nrow_override))
        df_plot_nrow_override <- nrow(df_plot)

    if (is.null(asp.use))
        asp.use <- ( 0.65 - 1.75) / (10 - 3) * df_plot_nrow_override + 2.35
    if (is.null(vertex.scale.use))
        vertex.scale.use <- (25 - 50) / (10 - 3) * df_plot_nrow_override + 65

    g <- igraph::make_graph(graph.vector, directed = FALSE)
    if (!is.null(outfile)) {
        pdf(outfile, width=5, height=5, family="Arial")
        plot.igraph(
          x                = g,
          xlim             = c(-1, 1.0),
          ylim             = c(-1, 1),
        asp=asp.use,

            axes=FALSE,
          layout           = as.matrix(layout),
          vertex.label     = NA,
          vertex.size      = scales::rescale(sqrt(vertex_weights / 3.14)) * vertex.scale.use,
          vertex.color     = cols.full,
          edge.width       = 1e-10,
          edge.color       = "#00000000")
        axis(3, at = seq(-1, 1, 2/(length(tick.names)-1)), labels=tick.names, cex.axis=8/12, las=2)
        axis(4, at = seq(-1, 1, 2/(dim(df_plot)[2]-1)), labels=dim(df_plot)[2]:1, cex.axis=8/12, las=2)
        mtext("Age range (months)", side=3, line=3)
        mtext("Cluster number", side=4, line=2)
        dev.off()
    } else {
        plot.igraph(
          x                = g,
          xlim             = c(-1, 1),
          ylim             = c(-1, 1),
        asp=asp.use,
            axes=FALSE,    
          layout           = as.matrix(layout),
          vertex.label     = NA,
          vertex.size      = scales::rescale(sqrt(vertex_weights / 3.14)) * vertex.scale.use,
          vertex.color     = cols.full,
          edge.width       = 1e-10,
          edge.color       = "#00000000",
            )
            axis(3, at = seq(-1, 1, 2/(length(tick.names)-1)), labels=tick.names, cex.axis=0.85)
        axis(4, at = seq(-1, 1, 2/(dim(df_plot)[2]-1)), labels=dim(df_plot)[2]:1, cex.axis=0.85)
        mtext("Age range (months)", side=3, line=3)
        mtext("Cluster number", side=4, line=2)


    }
}



cluster_species_heatmap <- function(df_species, assignment, order.use, n.clusters, outfile=NULL) {
    library('ComplexHeatmap')
    library('extrafont')
    df_species <- apply(df_species,1, function(x) x / sum(x)) #We need relative abundance and species X samples
    cols.order <- c()
    cols.cluster <- c()
    for (i in 1:n.clusters) {
        loc.cluster <- which(assignment==order.use[i])
        cols.order <- c(cols.order, loc.cluster)
        cols.cluster <- c(cols.cluster,rep(i, times=length(loc.cluster)))
    }
    stopifnot(length(cols.order) == ncol(df_species))
    df_plot <- df_species[,cols.order]
    rownames(df_plot) <- str_replace_all(rownames(df_plot),"s__","") %>% str_replace_all("_", " ")

    # From the teddy code
    colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")), bias=3)(100)

    # Cluster columns to colors
    cols <- inlmisc::GetTolColors(n = n.clusters)

    color.map <- cols
    names(color.map) <- 1:n.clusters
    cols.cluster <- factor(cols.cluster, levels=unique(cols.cluster), ordered=TRUE)

    ha = HeatmapAnnotation(cluster=cols.cluster, col=list(cluster=color.map), border=TRUE)
    rownames(df_species) <- str_replace_all(rownames(df_species),"s__","") %>% str_replace_all("_", " ")
    p <- Heatmap(df_plot, name="abundance", show_column_names=FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
                                            col=colors, top_annotation = ha, row_names_gp = gpar(fontsize = 8, fontface="italic"))
    
    if (!is.null(outfile)) {
        pdf(outfile, width=9.1, height=5, family="Arial")
        print(p)
        dev.off()
    } else
        print(p)
    
}



do_dmm_clustering <- function(df_model, species_use, clusters.test=1:15, method_use="laplace") {

    df_species <- subset_to_species(df_model, species_use)
    models <- try_n_clusters(df_species, clusters.test=clusters.test, verbose=FALSE)

    results <- identify_best_n_clusters(models, method_use, clusters.test)
    assignments <- assign_samples(df_species, results$model)
    n.clusters <- dim(results$model@group)[2]
    df_plot <- calc_samples_per_cluster(df_model, assignments, n.clusters)
    order.use <- df_plot %>% apply(2, function(x) which(x == max(x))[1]) %>% order
    df_plot <- df_plot[,order.use]
    return(list(df_plot=df_plot, order.use=order.use, plot=results$plot))

}


cluster_second_dataset <- function(df_model_second, model, species_oi, order.use) {
    df_model_dart <- add_missing_species(df_model_second, species_oi)
    df_model_dart_sub <- select(df_model_dart,age,Individual,one_of(species_oi))
    assignments <- assign_samples(df_model_dart_sub, model,TRUE)
    n.clusters <- dim(model@group)[2]
    df_plot_dart <- calc_samples_per_cluster(df_model_dart_sub, assignments, n.clusters)
    # Use the same order as before
    df_plot_dart <- df_plot_dart[,order.use]
    make_graph(df_plot_dart)
}


clusters_prevalence <- function(df_species, n.clusters, species_oi, assignments=NULL, model=NULL, order.use=NULL, reduce_names=FALSE) {
    source("/home/averster/Documents/RossLab/CF_Paper/CF_Metagenomics_Manuscript/av_code_for_figures/Lib/Figures.R")

    df_species <- subset_to_species(df_species, species_use)

    if (is.null(assignments))
        assignments <- assign_samples(df_species, model)

    stopifnot(length(assignments) == nrow(df_species))
    df_out <- NULL
    for (i in 1:n.clusters) {
        if (is.null(order.use)) {
            df_model_sub1 <- df_species[assignments==i,]
        } else {
            df_model_sub1 <- df_species[assignments==order.use[i],]
        }
        for (species in species_oi) {
              df_out <- bind_rows(df_out, prevalence_stats(df_model_sub1,species,i))
        }
    }
    colnames(df_out) <- c("species","cluster","n","prevalence","prevalence_se","abundance","abundance_se","nonzero_abundance","nonzero_abundance_se")

    df_out$species <- clean_up_metaphlan_names(df_out$species, reduce_names)
    species.order.use <- clean_up_metaphlan_names(species_oi, reduce_names)

    df_out <- replace(df_out, is.na(df_out), 0)

    df_out$species <- factor(df_out$species, ordered=TRUE, levels=species.order.use)
    df_out$cluster <- as.factor(df_out$cluster)

    return(df_out)

}



group_prevalence <- function(df_species, species_oi, df_groups, reduce_names=FALSE, min.n.data=10) {
    source("/home/averster/Documents/RossLab/CF_Paper/CF_Metagenomics_Manuscript/av_code_for_figures/Lib/Figures.R")

    df_species <- subset_to_species(df_species, species_oi)

    # Verify that the assignments correspond to the species
    stopifnot(nrow(df_species) == nrow(df_groups))
    stopifnot(all(rownames(df_species) == df_groups$samples))

    df_out <- NULL
    for (group in unique(df_groups$assignment)) {
        samples_oi <- filter(df_groups, assignment == group) %>% pull(samples)
        if (length(samples_oi) <= min.n.data)
            next
        
        df_species_sub <- df_species[samples_oi,]
        for (species in species_oi) {
              df_out <- bind_rows(df_out, prevalence_stats(df_species_sub,species,group))
        }
    }
    colnames(df_out) <- c("species","assignment","n","prevalence","prevalence_se","abundance","abundance_se","nonzero_abundance","nonzero_abundance_se")

    df_out$species <- clean_up_metaphlan_names(df_out$species, reduce_names)
    species.order.use <- clean_up_metaphlan_names(species_oi, reduce_names)

    df_out <- replace(df_out, is.na(df_out), 0)

    df_out$species <- factor(df_out$species, ordered=TRUE, levels=species.order.use)
    return(df_out)

}


clean_up_metaphlan_names <- function(sp_vector, reduce=TRUE) {
        sp_vector <- sp_vector %>% str_replace_all("s__","") %>% str_replace_all("_"," ")
        if (reduce)
            sp_vector <- str_replace(sp_vector, "([A-Z]{1})[a-z]+ ", "\\1. ")
        return(sp_vector)
}

prevalence_stats <- function(df_model, sp,i) {
    vals1 <- df_model[,sp]
    n1 <- nrow(df_model)
    d1 <- sum(vals1 > 0)
    p1 <- d1/n1
    p1.se <- binomial.se(vals1>0)
    return(data.frame(sp,i,n1,p1,p1.se,mean(vals1), se(vals1),mean(vals1[vals1>0]), se(vals1[vals1>0])))
}

order_species_by_cluster_differentiation <- function(df_model_CF, df_assignments_CF, species_use, groups_use=NULL, columns_group="cluster_ordered") {
    df_model_CF$samples <- rownames(df_model_CF)
    df_test <- left_join(df_model_CF, df_assignments_CF)
    if(!is.null(groups_use))
        df_test <- dplyr::filter(df_test, !!as.name(columns_group)%in%groups_use)

    df_results <- data.frame()
    for (sp in species_use[species_use%in%colnames(df_test)]) {
        pval <- kruskal.test(as.formula(paste0(sp, " ~ ", columns_group)), df_test)$p.value
        df_results <- bind_rows(df_results, data.frame(sp, pval))
    }

    df_results$padj <- p.adjust(df_results$pval, method="bonferroni")
    species_oi <- filter(df_results, padj < 0.05) %>% arrange(padj) %>% pull(sp)
    return(species_oi)
}

determine_species_for_group_differentiation <- function(df_cluster_assignments_cf, df_model_cf, species_use, method="clusters", min_n=25, max_return=10) {
    #TODO: CF vs non-CF?
    if (method == "clusters") {
        clusters_sub <- count(df_cluster_assignments_cf, cluster_ordered) %>% filter(n > min_n) %>% pull(cluster_ordered)
        species_oi <- order_species_by_cluster_differentiation(df_model_cf, df_cluster_assignments_cf, species_use, clusters_sub)
    } else {
        species_oi <- order_species_by_cluster_differentiation(df_model_cf, df_cluster_assignments_cf, species_use, columns_group="assignment")
    }
    if (length(species_oi) > max_return)
        species_oi <- species_oi[1:max_return]
    return(species_oi)
}

most_important_species_curve <- function(df_cluster_assignments_cf, df_cluster_assignments_healthy, df_model_cf, df_model_healthy, species_use, outfile, method="clusters", min_n=25) {
        library('extrafont')
	    
        t <- improve_plot()
        
        species_use <- determine_species_for_group_differentiation(df_cluster_assignments_cf, df_model_cf, species_use, method)

        n_clusters <- max(df_cluster_assignments_healthy$cluster_ordered)


        if (method == "clusters") {
            df_cluster_prevalence <- clusters_prevalence(df_model_cf, n_clusters, species_use, df_cluster_assignments_cf$cluster_ordered, reduce_names=TRUE) 
            df_cluster_prevalence$type <- "CF"
            df_cluster_prevalence_2 <- clusters_prevalence(df_model_healthy, n_clusters, species_use, df_cluster_assignments_healthy$cluster_ordered, reduce_names=TRUE)
            df_cluster_prevalence_2$type <- "Healthy"
        } else {
            df_cluster_prevalence <- group_prevalence(df_model_cf, species_use, df_cluster_assignments_cf, reduce_names=TRUE) 
            df_cluster_prevalence$type <- "CF"
            df_cluster_prevalence_2 <- group_prevalence(df_model_healthy, species_use, df_cluster_assignments_healthy, reduce_names=TRUE)
            df_cluster_prevalence_2$type <- "Healthy"
        }
        df_plot <- bind_rows(df_cluster_prevalence, df_cluster_prevalence_2)

        if (method == "clusters")  {
            clusters_sub <- count(df_cluster_assignments_cf, cluster_ordered) %>% filter(n > min_n) %>% pull(cluster_ordered)
            df_plot <- filter(df_plot, cluster%in%clusters_sub)
            p <- ggplot(df_plot, aes(x = cluster, y = prevalence))
        } else {
            p <- ggplot(df_plot, aes(x = assignment, y = prevalence))
        }
        p <- p + geom_line(lwd=0.75,aes(group=type, color=type)) + geom_point(aes(color=type), size=0.75) + scale_y_continuous(labels = scales::percent) + xlab("Cluster") + ylab("Prevalence") + scale_color_manual(values = c("steelblue", "#c33b3b")) +
                    geom_errorbar(aes(ymin=prevalence-prevalence_se, ymax=prevalence+prevalence_se,group=type, color=type), width=0.22) + facet_wrap(. ~ species, ncol=5) + theme_bw(14) + t + theme(strip.text = element_text(size = 12, face="italic"), legend.title=element_blank())

        if (method != "clusters") {
            p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) )
        }

        pdf(outfile, width = 10, height=4, family="Arial")
        print(p)
        dev.off()

}


model_selection_curve <- function(models, outfile) {

    source("/home/averster/Documents/RossLab/CF_Paper/CF_Metagenomics_Manuscript/av_code_for_figures/Lib/Figures.R")
    library('gridExtra')
    library('extrafont')

    t <- improve_plot()
    r <- identify_best_n_clusters(models,"laplace")
    best.score <- filter(r$df_plot, n_test == r$n.clusters) %>% pull(score)

    p1 <- ggplot(r$df_plot, aes(x=n_test, y=score)) + geom_line() + geom_point(pch=21, size=3, fill="steelblue") + theme_bw(10) + t + xlab("Number of clusters")
    p1 <- p1 + ylab("Laplace") + annotate("text", x=r$n.clusters, y=best.score*1.005, size=16,label="*") + scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
            theme(axis.text=element_text(size=8))

    r <- identify_best_n_clusters(models,"AIC")
    best.score <- filter(r$df_plot, n_test == r$n.clusters) %>% pull(score)

    p2 <- ggplot(r$df_plot, aes(x=n_test, y=score)) + geom_line() + geom_point(pch=21, size=3, fill="steelblue") + theme_bw(10) + t + xlab("Number of clusters")
    p2 <- p2 + ylab("AIC") + annotate("text", x=r$n.clusters, y=best.score*1.005, size=16,label="*") + scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
            theme(axis.text=element_text(size=8))

    r <- identify_best_n_clusters(models,"BIC")
    best.score <- filter(r$df_plot, n_test == r$n.clusters) %>% pull(score)

    p3 <- ggplot(r$df_plot, aes(x=n_test, y=score)) + geom_line() + geom_point(pch=21, size=3, fill="steelblue") + theme_bw(10) + t + xlab("Number of clusters")
    p3 <- p3 + ylab("BIC") + annotate("text", x=r$n.clusters, y=best.score*1.005, size=16,label="*") + scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
            theme(axis.text=element_text(size=8))

    pdf(outfile, width=6, height=5.1, family="Arial")
    grid.arrange(p1,p2,p3,ncol=2)
    dev.off()
}



multiple_dmm_runs <- function(df_species, outfile, df_plot=NULL) {
  library('DirichletMultinomial')
  library('gridExtra')
  library('extrafont')
  t <- improve_plot()

  df_species <- subset_to_species(na.omit(df_model_healthy), species_use)

  set.seed(123)

  clusters.test <- 1:15
  
  if (is.null(df_plot)) {
    df_plot <- data.frame()
    for (i in 1:5) {

        fit <- lapply(clusters.test, dmn, count=data.matrix(df_species), verbose=TRUE)
        lplc <- sapply(fit, BIC)
        names(lplc) <- clusters.test
        names(fit) <- clusters.test
        print(sprintf("Best according to BIC is %i",which(lplc == min(lplc))))

        df_plot <- bind_rows(df_plot, data.frame(n_test=clusters.test, laplace=lplc, rep=i))
    }
  }

  t <- improve_plot()

  df_annotate <- group_by(df_plot, rep) %>% summarize(y=min(laplace),x=which(laplace == min(laplace)))
  df_annotate <- df_annotate %>% group_by(x) %>% summarize(y=y,adjustment=cumsum(x == x)-1)

  p <- ggplot(df_plot,aes(x=n_test,y=laplace)) + geom_line(aes(group=rep), alpha=0.3) + theme_bw(12) + ylab("BIC") +
      geom_text(data=df_annotate, aes(x=x, y=y+500+adjustment*1000, label="*"), colour="red", size=16) + xlab("Number of clusters") + t

    pdf(outfile, width=3.01, height=2.65, family="Arial")
    print(p)
    dev.off()


}
