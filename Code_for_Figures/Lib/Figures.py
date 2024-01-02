
import argparse
from pathlib import Path
import pickle
import Figure_Data

import numpy as np
import pandas as pd
import sys
import metagenomic_strains

from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()

NTREE=10000
ro.globalenv['NTREE'] = NTREE


def identify_agemodel_species(df_model_Healthy, metric_cutoff=0.95):
    ro.globalenv['df_model'] = df_model_Healthy
    ro.globalenv['metric_cutoff'] = metric_cutoff
    ro.globalenv['ntree'] = NTREE

    ro.r('''

        library('tidyverse')
        library('RRF')
        source("/home/averster/Documents/RossLab/Lib/InfantAgeModels.R")

        set.seed(123)
        species.oi <- identify_subspecies(na.omit(df_model), metric.cutoff=metric_cutoff, ntree=NTREE)
    ''')
    return(ro.globalenv['species.oi'])


def calc_relative_age(df_model_Healthy, df_model_CF, species_oi=None):
    ro.globalenv['df_model_Healthy'] = df_model_Healthy
    ro.globalenv['df_model_CF'] = df_model_CF
    if species_oi is None:
        species_oi = ro.r('NULL')
    ro.globalenv['species.oi'] = species_oi

    ro.r('''
        source("Lib/InfantAgeModels.R")
        library('tidyverse')

        df_model_Healthy <- na.omit(df_model_Healthy)
        if (!is.null(species.oi)) {
            df_model_CF <- add_missing_species(df_model_CF, species.oi)
            df_model_Healthy <- select(df_model_Healthy,one_of(species.oi),age,Individual)
            df_model_CF <- select(df_model_CF,one_of(species.oi),age, Individual)
        }

        df_tax <- add_species_zeros(list(df_model_Healthy, df_model_CF))
        df_model_Healthy_sub <- df_tax[[1]]
        df_model_CF_sub <- df_tax[[2]]
        save(list=ls(), file="debug.RData")

        age.model.full <- relative_microbial_age(df_model_Healthy_sub, df_model_CF_sub, ntree=NTREE)
        df_relative_age_CF <- select(df_model_CF, Individual) %>% rownames_to_column %>% left_join(rownames_to_column(age.model.full$df_plot))

        age.model.full <- relative_microbial_age(df_model_Healthy_sub, df_model_Healthy_sub, ntree=NTREE)
        df_relative_age <- select(df_model_Healthy, Individual) %>% rownames_to_column %>% left_join(rownames_to_column(age.model.full$df_plot))

    ''')
    return ro.globalenv['df_relative_age'], ro.globalenv['df_relative_age_CF']


def age_model_simplified_species(df_model_Healthy, species_oi, outfile, species_oi_alt=None):
    ro.globalenv['df_model_Healthy'] = df_model_Healthy
    ro.globalenv['species.oi'] = species_oi
    ro.globalenv['outfile'] = outfile
    if species_oi_alt is not None:
        ro.globalenv['species_oi_alt'] = species_oi_alt
    else:
        ro.globalenv['species_oi_alt'] = ro.r('NULL')

    ro.r('''
        library('gridExtra')
        source("Lib/Figures.R")
        t <- improve_plot()

        df_model_Healthy <- na.omit(df_model_Healthy)
        age_model_simplified_species(df_model_Healthy, species_oi, species_oi_alt=species_oi_alt, ntree=NTREE, outfile=outfile)
    ''')

def plot_relative_age(df_relative_age, outfile, indv=None):

    ro.globalenv['df_relative_age'] = df_relative_age
    ro.globalenv['outfile'] = outfile
    if indv is not None:
        ro.globalenv['indv'] = indv

    # TODO: Justify age < 30 only

    ro.r('''
        library('extrafont')

        # Get the pvalues
        pval.1 <- summary(lm(relative_age ~ age, filter(df_relative_age, age<18)))[[4]][2,4]
        pval.2 <- summary(lm(relative_age ~ age, filter(df_relative_age, (age>=18) & (age < 30))))[[4]][2,4]

        t <- improve_plot()

        if (exists("indv")) {
            df_plot <- df_relative_age %>% filter(Individual == indv)
            p <- ggplot(df_plot, aes(x = age, y = relative_age)) + geom_line(alpha=0.5) + geom_point(size=3, pch=21, fill="steelblue")  + theme_bw(10) + t + ylab("Relative Microbiome Age") + xlab("Age (months)") + ggtitle(sprintf("Indiv: %i",indv.oi))
        } else {
            df_plot <- df_relative_age
            p <- ggplot(filter(df_relative_age, age<30), aes(x = age, y = relative_age)) + geom_point(pch=21, fill = "steelblue") + theme_bw(10) + t + geom_smooth(lwd=2, color="purple",se=FALSE, method="loess") + ylab("Relative microbiota age") + xlab("Age (months)")
            p <- p + geom_vline(xintercept = 18, linetype="dashed", color="red") + annotate("text", x=4.25, y=7.5, size=8/.pt, label=sprintf("p = %s",formatC(pval.1, format="E", digits=2))) + annotate("text", x=27, y=7.5, size=8/.pt, label=sprintf("p = %.2f",pval.2))
        }

    pdf(outfile, width=3.3,height=3.3,family="Arial")
    print(p)
    dev.off()

    png(str_replace_all(outfile,".pdf",".png"), width=3.3, height=3.3, units="in", res=300, family="Arial")
    print(p)
    dev.off()

    ''')


def species_agemodel_abundanceplot(df_model_Healthy, species_oi, outfile):
    ro.globalenv['df_model_Healthy'] = df_model_Healthy
    ro.globalenv['species.oi'] = species_oi
    ro.globalenv['outfile'] = outfile

    ro.r('''
        library('gridExtra')
        library('extrafont)
        source('Lib/Figures.R')
        library('tidyverse')
        t <- improve_plot()

        df_plot <- age_model_species_properties(df_model_Healthy, species.oi)
        p1 <- ggplot(df_plot, aes(x=type,y=log10abundance)) + geom_boxplot(fill="steelblue", outlier.shape=NA) + geom_jitter(size=2.5, fill="grey", pch=21) + 
                            theme_bw(16) + t + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")
        p2 <- ggplot(df_plot, aes(x=type,y=prevalence)) + geom_boxplot(fill="steelblue", outlier.shape=NA) + geom_jitter(size=2.5, fill="grey", pch=21) + 
                            theme_bw(16) + t + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + scale_y_continuous(labels = scales::percent)
        pdf(outfile, width=5, height=4, family="Arial")
        grid.arrange(p1,p2, ncol=2)
        dev.off()
    ''')

def prevalence_plot(df_model_healthy, species_oi, outfile):


    ro.globalenv['df_model_healthy'] = df_model_healthy
    ro.globalenv['species.oi'] = species_oi
    ro.globalenv['outfile'] = outfile
    ro.r('''
        source("Lib/Figures.R")
        library('extrafont')
        df_prevalence_wt <- calc_prevalence_bins(na.omit(df_model_healthy), species.oi, max_age=36)

        t <- improve_plot()
        df_prevalence_wt$species <- str_replace_all(df_prevalence_wt$species, "s__","") %>% str_replace_all("_"," ")

        species.order.use <- df_prevalence_wt %>% group_by(species) %>% summarize(prevalence=mean(prevalence)) %>% arrange(desc(prevalence)) %>% pull(species)

        df_prevalence_wt$age_range <- factor(df_prevalence_wt$age_range, ordered=TRUE, levels=unique(df_prevalence_wt$age_range))
        df_prevalence_wt$species <- factor(df_prevalence_wt$species, ordered=TRUE, levels=species.order.use)

        p_errorbars <- ggplot(df_prevalence_wt, aes(x = age_range, y = prevalence)) + geom_line(lwd=1.25,aes(group=species, color=species)) + scale_y_continuous(labels = scales::percent) +
                geom_errorbar(aes(ymin=prevalence-prevalence_se, ymax=prevalence+prevalence_se,group=species,color=species), width=1.25) + theme_bw(16) + t + scale_colour_brewer(palette = "Set3")

        p <- ggplot(df_prevalence_wt, aes(x = age, y = prevalence)) + geom_line(lwd=1.25,aes(group=species, color=species)) + scale_y_continuous(labels = scales::percent) + ylab("Prevalence") +
                 theme_bw(14) + t + scale_colour_brewer(palette = "Set3") + theme(legend.text=element_text(face="italic", size=8))
        pdf(outfile, width=6, height=3.5, family="Arial")
        print(p)
        dev.off()
    ''')


def individual_relative_age(df_relative_age, df_model_dart, indv_oi, species_oi, outfile):

    ro.globalenv['df_relative_age'] = df_relative_age
    ro.globalenv['df_model_dart'] = df_model_dart
    ro.globalenv['indv.oi'] = indv_oi
    ro.globalenv['species_oi'] = species_oi
    ro.globalenv['outfile'] = outfile
    ro.r('''
        library('tidyverse')
        library('extrafont')

        t <- improve_plot()
        df_model_dart <- add_missing_species(df_model_dart, species.oi)
        df_plot <- df_relative_age %>% filter(Individual == indv.oi)

        df_plot_species <- df_model_dart %>% filter(Individual==indv.oi) %>% select(one_of(species_oi), age) %>% arrange(age) %>% pivot_longer(!age)
        df_plot_species$name <- str_replace_all(df_plot_species$name, "s__","") %>% str_replace_all("_", " ")

        df_plot_species$facet <- "species"
        df_plot$facet <- "relative_age"

        df_plot_sub <- select(df_plot, age, Individual, relative_age, facet) %>% rename(name=Individual, value=relative_age)

        df_plot_full <- bind_rows(df_plot_sub, df_plot_species)

        facet_names <- c('relative_age' = "Relative Age (Months)",
                            'species' = "Species Abundance")

        sp.order.use <- df_plot_species %>% group_by(name) %>% summarize(ab=max(value)) %>% arrange(desc(ab)) %>% pull(name)
        sp.order.use <- c(sp.order.use,"116")
        df_plot_full$name <- factor(df_plot_full$name, ordered=TRUE, levels=sp.order.use)

        p <- ggplot(df_plot_full, aes(x = age,y=value)) + geom_line(data=filter(df_plot_full,facet=="relative_age"),alpha=0.5) + geom_point(data=filter(df_plot_full,facet=="relative_age"), pch=21, fill="steelblue",size=3) +
                        geom_line(data=filter(df_plot_full,facet=="species"), aes(group=name, color=name), lwd=1.5) + facet_grid(facet~.,scales="free_y",labeller=as_labeller(facet_names)) + theme_bw(12) + t + ylab("") + xlab("Age (months)") + theme(legend.title=element_blank()) + ggtitle(sprintf("Individual %i",indv.oi)) + scale_colour_brewer(palette = "Set3")


        pdf(outfile, width=9, height=4.5, family="Arial")
        print(p)
        dev.off()
    ''')


def prevalence_plot_compare(df_model_Healthy, df_model_CF, species_oi, outfile=None, df_model_Healthy_Alt=None, hayden_bins=False):

    ro.globalenv['df_model_Healthy'] = df_model_Healthy
    ro.globalenv['df_model_CF'] = df_model_CF
    ro.globalenv['species.oi'] = species_oi
    ro.globalenv['hayden_bins'] = hayden_bins
    if outfile is None:
        outfile = ro.r('NULL')
    ro.globalenv['outfile'] = outfile

    if df_model_Healthy_Alt is None:
        df_model_Healthy_Alt = ro.r('NULL')
    
    ro.globalenv['df_model_Healthy_Alt'] = df_model_Healthy_Alt

    ro.r('''
        source("Lib/Figures.R")
        prevalence_plot(df_model_Healthy, df_model_CF, species.oi, hayden_bins, outfile, df_model_Healthy_Alt)

    ''')


def cf_indv_rescue(outfile_barplot, outfile_gridplot, indv=116, n_comb=3):

    ro.globalenv['df_model_Healthy'] = df_model_Healthy
    ro.globalenv['df_model_CF'] = df_model_CF
    ro.globalenv['species.oi'] = species_oi
    ro.globalenv['indv'] = indv
    ro.globalenv['n_comb'] = n_comb
    ro.globalenv['outfile_barplot'] = outfile_barplot
    ro.globalenv['outfile_gridplot'] = outfile_gridplot

    ro.r('''

        library('gridExtra')
        library('extrafont')

        n_plot <- 10

        results <- all_species_rescue(species.oi, n_comb, df_model_CF, na.omit(df_model_Healthy), indv, age.model.full)
        graphs_list <- results[str_detect(names(results),"graph")]
        species_detected <- results[str_detect(names(results),"species")]
        df_plot <- summarize_species_detected(species_detected, species.oi, n_comb)
        p <- ggplot(df_plot, aes(x = sp, y = frac)) + geom_bar(stat="identity", color="black",fill="steelblue") + theme_bw(20) + t + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab(sprintf("Percentage of %i species combinations",n_comb)) + xlab("") + scale_y_continuous(labels = scales::percent)
        pdf(outfile_barplot, width=3.75, height=3.75)
        print(p)
        dev.off()


        if (length(graphs_list) > n_plot)
            graphs_list <- graphs_list[sample(length(graphs_list), n_plot, replace=FALSE)]
        else
            grid.arrange(grobs=graphs_list, ncol=3)

        pdf(outfile_gridplot, width=7, height=7, family="Arial")
        print(grid.arrange(grobs=graphs_list, ncol=3))
        dev.off()

    ''')


def variable_importance_overlap(df_model_healthy, df_model_secondary, species_oi, outfile_overlap, outfile_venn, ds1="Hayden",ds2="Diabimmune", filter_age=True):
    ro.globalenv['df_model_healthy'] = df_model_healthy
    ro.globalenv['df_model_secondary'] = df_model_secondary
    ro.globalenv['outfile_overlap'] = outfile_overlap
    ro.globalenv['outfile_venn'] = outfile_venn
    ro.globalenv['filter_age'] = filter_age
    ro.globalenv['species_oi'] = species_oi
    ro.globalenv['ds1'] = ds1
    ro.globalenv['ds2'] = ds2

    ro.r('''
        library('tidyverse')
        source("Lib/Figures.R")
        limits.use <- c(1e-3, 4e5)
        if (filter_age) {
            df_model_secondary <- dplyr::filter(df_model_secondary, age <= 12)
            df_model_healthy <- dplyr::filter(df_model_healthy, age <= 12)
            limits.use <- c(1e-3, 4e3)
        }
        df_importance <- variable_importance_overlap(na.omit(df_model_healthy), na.omit(df_model_secondary), NTREE)
        print(df_importance)
        importance_overlap_plot(df_importance, ds1, ds2, limits.use, species_oi=species_oi, outfile=outfile_overlap)
        venn_diagram_species_overlap(df_importance, outfile_venn, ds1, ds2)
    ''')


def predictions_crossdataset(df_model_healthy, df_model_secondary, ds1, ds2, outfile, filter_age=True):
    ro.globalenv['df_model_healthy'] = df_model_healthy
    ro.globalenv['df_model_secondary'] = df_model_secondary
    ro.globalenv['ds1'] = ds1
    ro.globalenv['ds2'] = ds2
    ro.globalenv['outfile'] = outfile
    ro.globalenv['filter_age'] = filter_age

    ro.r('''
        source("Lib/Figures.R")
        df_model_healthy <- na.omit(df_model_healthy)
        df_model_secondary <- na.omit(df_model_secondary)
        max.use <- 45
        if (filter_age) {
         df_model_secondary <- df_model_secondary %>% filter(age <= 12)
         df_model_healthy <- df_model_healthy %>% filter(age <= 12)
         max.use <- 12
        }
        
        cross_dataset_predictions(df_model_healthy, df_model_secondary, ds1,ds2, outfile, max.use,NTREE)
    ''')


def dmm_clustering(df_model_healthy, df_model_cf, species_use, outfile_graph, outfile_assignments, outfile_graph_cf, outfile_assignments_cf, outfile_heatmap, outfile_model_selection, outfile_model_selection_repeat, outfile_linegraph, outfile_model, method_clusters_use="BIC"):
    ro.globalenv['df_model_healthy'] = df_model_healthy
    ro.globalenv['df_model_cf'] = df_model_cf
    ro.globalenv['species_use'] = species_use
    ro.globalenv['method_clusters_use'] = method_clusters_use
    ro.globalenv['outfile_graph'] = outfile_graph
    ro.globalenv['outfile_assignments'] = outfile_assignments
    ro.globalenv['outfile_heatmap'] = outfile_heatmap
    ro.globalenv['outfile_model_selection'] = outfile_model_selection
    ro.globalenv['outfile_model_selection_repeat'] = outfile_model_selection_repeat    
    ro.globalenv['outfile_linegraph'] = outfile_linegraph
    ro.globalenv['outfile_graph_cf'] = outfile_graph_cf
    ro.globalenv['outfile_assignments_cf'] = outfile_assignments_cf
    ro.globalenv['outfile_model'] = outfile_model
    ro.r('''
    
        source("Lib/DMM_Clustering.R")
        source("Lib/Figures.R")

        # The only way to get the graphs to come out the same size is to hardcode this. You might need to mess with it a little to improve your plots.
        nrow_use <- 11

        set.seed(121)
        df_species <- subset_to_species(na.omit(df_model_healthy), species_use)

        models <- try_n_clusters(df_species, clusters.test=1:15, verbose=FALSE)
        r <- identify_best_n_clusters(models,method_clusters_use, clusters.test=1:15)
        print(sprintf("Number of clusters is %i", r$n.clusters))
        
        #model_selection_curve(models, outfile_model_selection)
        #multiple_dmm_runs(df_species, outfile_model_selection_repeat, df_plot=NULL)

        model_best <- r$model.best
        n_clusters <- r$n.clusters
        assignments <- assign_samples(df_species, model_best)
        df_plot <- calc_samples_per_cluster(df_model_healthy, assignments, n_clusters)

        order.use <- order_clusters(df_plot)
        df_plot <- df_plot[,order.use]

        save.image("DMM_RObj.RData")
        make_graph(df_plot, outfile = outfile_graph, df_plot_nrow_override=nrow_use)

        cluster_species_heatmap(df_species, assignments, order.use, n_clusters, outfile_heatmap)

        df_cluster_assignments <- data.frame(samples=rownames(df_species), cluster=assignments)
        df_cluster_orders <- data.frame(cluster_ordered=1:length(order.use), cluster=order.use)
        df_cluster_assignments <- left_join(df_cluster_assignments, df_cluster_orders)

        write_csv(df_cluster_assignments, outfile_assignments)
        
        # Show the CF version
        df_cf_species <- subset_to_species(df_model_cf, species_use)

        assignments.cf <- assign_samples(df_cf_species, model_best, predict=TRUE)

        age.range.use <- seq(from=0,to=24,by=3)
        df_plot <- calc_samples_per_cluster(df_model_cf, assignments.cf, n_clusters, min_datapts=15, age.range=age.range.use)
        df_plot <- df_plot[,order.use]
        make_graph(df_plot, outfile = outfile_graph_cf, df_plot_nrow_override=nrow_use)

        save(df_species, df_plot, n_clusters, model_best, order.use, df_model_cf, assignments.cf, file=outfile_model)


        df_cluster_assignments_cf <- data.frame(samples=rownames(df_cf_species), cluster=assignments.cf)
        df_cluster_assignments_cf <- left_join(df_cluster_assignments_cf, df_cluster_orders)

        write_csv(df_cluster_assignments_cf, outfile_assignments_cf)


        most_important_species_curve(df_cluster_assignments_cf, df_cluster_assignments, df_cf_species, df_species, species_use, outfile_linegraph)
        print("On the igraph, shrink the axis labels from 12 to 10 font size. The code wont do it for some reason.")

    ''')


def dmm_important_species(df_model_healthy, df_model_cf, species_use, infile_assignments, infile_assignments_cf, infile_map, outfile_linegraph):

    ro.globalenv['df_model_healthy'] = df_model_healthy
    ro.globalenv['df_model_cf'] = df_model_cf
    ro.globalenv['species_use'] = species_use


    ro.globalenv['infile_assignments'] = infile_assignments
    ro.globalenv['infile_assignments_cf'] = infile_assignments_cf
    ro.globalenv['infile_map'] = infile_map

    ro.globalenv['outfile_linegraph'] = outfile_linegraph

    ro.r('''
        library('tidyverse')
        source("Lib/DMM_Clustering.R")

        df_cluster_assignments <- read_csv(infile_assignments)
        df_cluster_assignments_cf <- read_csv(infile_assignments_cf)
        df_cluster_assignments_map <- read_csv(infile_map)

        df_cluster_assignments <- left_join(df_cluster_assignments, df_cluster_assignments_map)
        df_cluster_assignments_cf <- left_join(df_cluster_assignments_cf, df_cluster_assignments_map)

        df_species <- subset_to_species(na.omit(df_model_healthy), species_use)
        df_cf_species <- subset_to_species(df_model_cf, species_use)

        most_important_species_curve(df_cluster_assignments_cf, df_cluster_assignments, df_cf_species, df_species, species_use, outfile_linegraph, method="groups")
    ''')


def determine_age_group(age, groups = [(0,12),(12,18),(18,np.inf)]):
    names = {0: "Developmental",1:"Transitional",2:"Stable"}
    for i, g in enumerate(groups):
        if (age > g[0]) & (age <= g[1]):
            return names[i]
    assert False #We should not get here


def fp_strains_barplot(df_healthy, df_cf, infile_strains_healthy, infile_strains_cf, plot_do="healthy", outfile="temp.pdf"):

    df_assignments = pd.DataFrame({"samples": df_healthy.index, "assignment":[determine_age_group(x) for x in df_healthy["age"]]})
    df_Fp = metagenomic_strains.load_fpmarker_counts(infile_strains_healthy)
    df_plot = metagenomic_strains.combine_data(df_Fp, None, df_assignments)
    df_plot["type"] = "Healthy"

    df_assignments_cf = pd.DataFrame({"samples": df_cf.index, "assignment":[determine_age_group(x) for x in df_cf["age"]]})
    df_Fp_cf = metagenomic_strains.load_fpmarker_counts(infile_strains_cf)
    df_plot_cf = metagenomic_strains.combine_data(df_Fp_cf, None, df_assignments_cf)
    df_plot_cf["type"] = "CF"

    df_plot = pd.concat([df_plot, df_plot_cf])

    if plot_do == "healthy":
        df_plot = df_plot.query('type == "Healthy"')
        df_plot_summary = metagenomic_strains.summarize_clades(df_plot)
        df_plot_summary["plot_var"] = df_plot_summary["category"]
        df_pvals = metagenomic_strains.calculate_pvals(df_plot)
    elif plot_do == "middle":
        df_plot = df_plot.query('assignment == "middle"')
        df_plot_summary = metagenomic_strains.summarize_clades(df_plot)
        df_plot_summary["plot_var"] = df_plot_summary["type"]
        df_pvals = metagenomic_strains.calculate_pvals(df_plot, "type")
    elif plot_do == "late":
        df_plot = df_plot.query('assignment == "late"')
        df_plot_summary = metagenomic_strains.summarize_clades(df_plot)
        df_plot_summary["plot_var"] = df_plot_summary["type"]
        df_pvals = metagenomic_strains.calculate_pvals(df_plot, "type")
    elif plot_do == "all":
        df_plot_summary = metagenomic_strains.summarize_clades(df_plot, group_vars = ["type"])
        df_plot_summary["plot_var"] = df_plot_summary["type"]
        df_pvals = metagenomic_strains.calculate_pvals(df_plot, "type")

    ro.globalenv['df_plot_summary'] = df_plot_summary
    ro.globalenv['df_pvals'] = df_pvals
    ro.globalenv['outfile'] = outfile

    ro.r('''
    
        library('tidyverse')
        library('ggsignif')
        library('extrafont')

        source("Lib/Figures.R")
        t <- improve_plot()

        df_pvals$pval_adjust <- p.adjust(df_pvals$pval, "fdr")
        print(df_pvals)

        df_plot_summary_long <- pivot_longer(df_plot_summary, cols=A:L, names_to="strain", values_to="detected")
        if ("CF"%in%df_plot_summary_long$plot_var) {
            df_plot_summary_long$plot_var <- factor(df_plot_summary_long$plot_var, levels=c("Healthy","CF"))
        } else {
            stopifnot(all(c("Developmental","Transitional","Stable")%in%df_plot_summary_long$plot_var))
            df_plot_summary_long$plot_var <- factor(df_plot_summary_long$plot_var, levels=c("Developmental","Transitional","Stable"))

        }

        print(head(df_plot_summary_long))
        p <- ggplot(df_plot_summary_long, aes(x=strain,y=detected)) + geom_bar(aes(fill=plot_var,group=plot_var),stat="identity",color="black", position="dodge") +
            theme_bw(10) + t + scale_y_continuous(labels = scales::percent, limits=c(0,0.9)) + ylab("Prevalence") +  
            theme(panel.border=element_rect(color="black",size=1.15), legend.title=element_blank(), legend.text=element_text(size=8),legend.position=c(0.775,0.9), legend.background = element_rect(fill = "white", color = "black"), legend.key.size = unit(0.3, "cm"))
        if ("CF"%in%df_plot_summary_long$plot_var) {
            p <- p + scale_fill_manual(values = c("#ff7f00","#1f78b4"))
        } else {
            p <- p + scale_fill_brewer(palette = "Set2")
        }
        df_pvals_sub <- filter(df_pvals, pval_adjust < 0.05)
        if (nrow(df_pvals_sub) > 0) {
            for (i in 1:nrow(df_pvals_sub)) {
                
                clade_oi <- df_pvals_sub[i,"clade"]
                pval.use <- df_pvals_sub[i,"pval_adjust"]
                x.use <- which(df_pvals$clade == clade_oi)
                
                y.use <- filter(df_plot_summary_long, strain == clade_oi) %>% pull(detected) %>% max
                y.use <- y.use + 0.03
                
                if (pval.use < 1e-4) {
                    annotation.use <- "***"
                } else if (pval.use < 1e-2) {
                    annotation.use <- "**"
                } else if (pval.use < 0.05) {
                    annotation.use <- "*"
                }
                
                p <- p + geom_signif(annotation = annotation.use,
                    y_position = y.use, xmin = x.use-0.3, xmax = x.use+0.3,
                    tip_length = c(0.00, 0.00), textsize=6, vjust=0.5) 

            }
        }

        pdf(outfile, width=4, height=3, family="Arial")
        print(p)
        dev.off()
            
    ''')


def write_fp_tree():

    infile_fp_markers = "/home/averster/Documents/RossLab/Sequences/faecalibacterium_12clades_markers_FINAL_DeFilippis2020.fna"
    indir = Path("/home/averster/Documents/RossLab/Genomes/Fp_Genbank")
    clade_marker_counts = metagenomic_strains.load_marker_count(infile_fp_markers)

    infile_blast_save = "df_blast_save.p"
    if Path(infile_blast_save).exists():
        df_blast = pickle.load(open(infile_blast_save,"rb"))
    else:
        df_blast = metagenomic_strains.blast_for_markers(clade_marker_counts, infile_fp_markers, indir)
        pickle.dump(df_blast, open(infile_blast_save,"wb"))

    clade_dict, color_dict = metagenomic_strains.define_dicts(df_blast)

    assert df_blast.query('cladeL > 0.95').shape[0] > 0, "We should be finding at least one cladeL genome"
    metagenomic_strains.write_tree(df_blast, color_dict, clade_dict, "Figure_Fp_Tree_full.pdf", "Figure_Fp_Tree_zoom.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s","--infile_species_saved")
    parser.add_argument('--infile_species_alt')
    parser.add_argument("-f","--frac_cutoff",default=0.95)
    parser.add_argument("-g","--graphs_do",action="append",default=[])
    parser.add_argument("-c","--cf_dataset")
    parser.add_argument("--healthy_dataset")
    parser.add_argument("--healthy_secondary", default="")
    parser.add_argument("-a","--restrict_age", default=False, action="store_true")
    parser.add_argument("--infile_strains_healthy", default="TEDDY_counts_FpMarkersDeFilippis_combined.csv")
    parser.add_argument("--infile_strains_cf", default="dart_CF_counts_FpMarkersDeFilippis_combined.csv")
    args = parser.parse_args()

    if args.healthy_dataset == "diabimmune":
        df_model_healthy = Figure_Data.load_diabimmune_data()
    elif args.healthy_dataset == "hayden":
        df_model_healthy = Figure_Data.load_hayden("Healthy")
    elif (args.healthy_dataset == "teddy") | (args.healthy_dataset == "TEDDY"):
        df_model_healthy = Figure_Data.load_TEDDY()
    else:
        raise Exception("Invalid Healthy Dataset")

    if args.healthy_secondary == "diabimmune":
        df_model_healthy_secondary = Figure_Data.load_diabimmune_data()
    elif args.healthy_secondary == "hayden":
        df_model_healthy_secondary = Figure_Data.load_hayden("Healthy")
    elif (args.healthy_secondary == "teddy") | (args.healthy_secondary == "TEDDY"):
        df_model_healthy_secondary = Figure_Data.load_TEDDY()

    if args.cf_dataset == "dart":
        df_model_cf = Figure_Data.load_dartmouth()
        df_model_cf["Individual"] = df_model_cf["Individual"].astype(str)
    elif args.cf_dataset == "hayden":
        df_model_cf = Figure_Data.load_hayden("CF")

    if args.infile_species_saved is not None:
        if not Path(args.infile_species_saved).exists():
            species_oi = identify_agemodel_species(df_model_healthy, args.frac_cutoff)
            pickle.dump( species_oi, open( args.infile_species_saved, "wb" ) )
        else:
            species_oi = pickle.load( open( args.infile_species_saved, "rb" ) )
    else:
        species_oi = identify_agemodel_species(df_model_healthy)

    species_oi_alt = None
    if args.infile_species_alt is not None:
        species_oi_alt = pickle.load( open( args.infile_species_alt, "rb" ) )

    if ("hayden" in args.healthy_dataset) | ("hayden" in args.healthy_secondary):
        assert args.restrict_age


    if ("variable_importance" in args.graphs_do) | ("all" in args.graphs_do):
        outfile_overlap="Figure_VariableImportanceOverlap_{}_{}.pdf".format(args.healthy_dataset,args.healthy_secondary)
        outfile_venn="Figure_VariableImportanceOverlap_{}_{}_Venn.pdf".format(args.healthy_dataset,args.healthy_secondary)
        if args.restrict_age:
            outfile_overlap = outfile_overlap.replace(".pdf","_1yronly.pdf")
            outfile_venn = outfile_venn.replace(".pdf","_1yronly.pdf")
        variable_importance_overlap(df_model_healthy, df_model_healthy_secondary, species_oi, outfile_overlap, outfile_venn, ds1=args.healthy_dataset, ds2=args.healthy_secondary, filter_age=args.restrict_age)
    if ("prevalence_plot" in args.graphs_do) | ("all" in args.graphs_do):
        prevalence_plot_compare(df_model_healthy, df_model_cf, species_oi, "Figure_Prevalence_{}_{}_Compare.pdf".format(args.healthy_dataset, args.cf_dataset))
        prevalence_plot_compare(df_model_healthy, df_model_cf, species_oi, "Figure_Prevalence_{}_{}_Compare_WithSecondary.pdf".format(args.healthy_dataset, args.cf_dataset), df_model_healthy_secondary)
        hayden_bins = False
        if "hayden" in args.healthy_secondary:
            hayden_bins = True
        prevalence_plot_compare(df_model_healthy, df_model_healthy_secondary, species_oi, "Figure_Prevalence_{}_{}_Compare.pdf".format(args.healthy_dataset, args.healthy_secondary), hayden_bins=hayden_bins)
        prevalence_plot(df_model_healthy, species_oi, "Figure_Prevalence_{}.pdf".format(args.healthy_dataset))
    if ("relative_age" in args.graphs_do) | ("all" in args.graphs_do):
        df_relative_age, df_relative_age_CF = calc_relative_age(df_model_healthy, df_model_cf)
        df_relative_age.to_csv("Relative_age_{}_{}_CF.csv".format(args.healthy_dataset, args.cf_dataset), index=False)
        df_relative_age.to_csv("Relative_age_{}_{}_Healhty.csv".format(args.healthy_dataset, args.cf_dataset), index=False)
        plot_relative_age(df_relative_age_CF, "Figure_Relative_Age_pred_{}_train_{}.pdf".format(args.cf_dataset, args.healthy_dataset))
        df_relative_age, df_relative_age_CF = calc_relative_age(df_model_healthy, df_model_cf, species_oi=species_oi)
        individual_relative_age(df_relative_age_CF, df_model_cf, 116, species_oi,"Figure_116_age_{}_train_{}_species.pdf".format(args.healthy_dataset,len(species_oi)))
        plot_relative_age(df_relative_age_CF, "Figure_Relative_Age_pred_{}_train_{}_SubSpecies.pdf".format(args.cf_dataset, args.healthy_dataset))
    if ("relative_age_predictions" in args.graphs_do) | ("all" in args.graphs_do):
        outfile_crossdataset_predictions = "Figure_PredictionsCrossDataset_{}_{}.pdf".format(args.healthy_dataset,args.healthy_secondary)
        outfile_simplified = "Figure_Agemodel_Simplified_{}.pdf".format(args.healthy_dataset)
        if args.restrict_age:
            outfile_crossdataset_predictions = outfile_crossdataset_predictions.replace(".pdf","_1yronly.pdf")
            outfile_simplified = outfile_simplified.replace(".pdf","_1yronly.pdf")

        age_model_simplified_species(df_model_healthy, species_oi, outfile_simplified, species_oi_alt)
        predictions_crossdataset(df_model_healthy, df_model_healthy_secondary, args.healthy_dataset,args.healthy_secondary,outfile_crossdataset_predictions, args.restrict_age)

    if ("DMM" in args.graphs_do) | ("all" in args.graphs_do):
        if args.healthy_dataset == "teddy":
            df_model_healthy = Figure_Data.load_TEDDY(use_counts=True)
        elif args.healthy_dataset == "diabimmune":
            df_model_healthy = Figure_Data.load_diabimmune_data(use_counts=True)
        else:
            raise Exception("not implement args.healthy_dataset")

        if args.cf_dataset == "dart":
            df_model_cf = Figure_Data.load_dartmouth(True)
            df_model_cf["Individual"] = df_model_cf["Individual"].astype(str)

        dmm_clustering(df_model_healthy, df_model_cf, species_oi, "Figure_DMM_Igraph_{}.pdf".format(args.healthy_dataset), "DMM_Assignments_{}.csv".format(args.healthy_dataset),
                                                                  "Figure_DMM_Igraph_{}_CF.pdf".format(args.healthy_dataset), "DMM_Assignments_{}_CF.csv".format(args.healthy_dataset),
                                                                  "Figure_DMM_heatmap_{}.pdf".format(args.healthy_dataset), "Figure_DMM_Selection_{}.pdf".format(args.healthy_dataset),
                                                                  "Figure_DMM_Selection_Repeat_{}.pdf".format(args.healthy_dataset),"Figure_DMM_CF_ClusterImportance_{}.pdf".format(args.healthy_dataset), 
                                                                  "model_{}.RData".format(args.healthy_dataset))



    if ("DMM_posthoc" in args.graphs_do):
        dmm_important_species(df_model_healthy, df_model_cf, species_oi, infile_assignments, infile_assignments_cf, infile_map, outfile_linegraph)
    if ("fp_strains" in args.graphs_do):
        df_cluster_assignments, df_cluster_assignments_cf = metagenomic_strains.load_cluster_assignments("DMM_Assignments_{}.csv".format(args.healthy_dataset), "DMM_Assignments_{}_CF.csv".format(args.healthy_dataset), "DMM_Assignments_{}_map.csv".format(args.healthy_dataset))
        fp_strains_barplot(df_model_healthy, df_model_cf, args.infile_strains_healthy, args.infile_strains_cf, "healthy", "Figure_Fp_strains_DMM_Healthy_{}.pdf".format(args.healthy_dataset))
        #fp_strains_barplot(df_model_healthy, df_model_cf,  args.infile_strains_healthy, args.infile_strains_cf, "middle", "Figure_Fp_strains_DMM_middle_{}.pdf".format(args.healthy_dataset))
        #fp_strains_barplot(df_model_healthy, df_model_cf,  args.infile_strains_healthy, args.infile_strains_cf, "late", "Figure_Fp_strains_DMM_late_{}.pdf".format(args.healthy_dataset))
        fp_strains_barplot(df_model_healthy, df_model_cf,  args.infile_strains_healthy, args.infile_strains_cf, "all", "Figure_Fp_strains_DMM_all_{}.pdf".format(args.healthy_dataset))
        write_fp_tree()