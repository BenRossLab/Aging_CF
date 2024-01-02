# Persistent delay in maturation of the developing gut microbiota in infants with cystic fibrosis

## Abstract:

The healthy human infant gut microbiome undergoes stereotypical changes in taxonomic composition between birth and maturation to an adult-like stable state. During this time, extensive communication between microbiota and the host immune system contributes to health status later in life. Although there are many reported associations between microbiota compositional alterations and disease in adults, less is known about how microbiome development is altered in pediatric diseases. One pediatric disease linked to altered gut microbiota composition is cystic fibrosis (CF), a multi-organ genetic disease involving impaired chloride secretion across epithelia and heightened inflammation both in the gut and at other body sites. Here, we use shotgun metagenomics to profile the strain-level composition and developmental dynamics of the infant fecal microbiota from several CF and non-CF longitudinal cohorts spanning from birth to greater than 36 months of life. We identify a set of keystone species whose prevalence and abundance reproducibly define microbiota development in early life in non-CF infants, but are missing or decreased in relative abundance in infants with CF. The consequences of these CF-specific differences in gut microbiota composition and dynamics are a delayed pattern of microbiota maturation, persistent entrenchment in a transitional developmental phase, and subsequent failure to attain an adult-like stable microbiota. We also detect the increased relative abundance of oral-derived bacteria and higher levels of fungi in CF, features that are associated with decreased gut bacterial density in inflammatory bowel diseases. Our results define key differences in the gut microbiota during ontogeny in CF and suggest the potential for directed therapies to overcome developmental delays in microbiota maturation.

## Main Figures

### [Figure 1:](Figures/Figure_1.pdf) Compositional differences in the CF Gut Microbiome

This figure shows the composition of the samples from each of our studies in an attempt to highlight significant differences not only between CF and healthy, but also between healthy cohorts.

A. First we show composition at the phylum level, organized by abundance of Proteobacteria for each cohort.

- Phylum plot was generated from `Lib/proteobacteria_plot_scripts.R`

B-E. PCoA plots of Dartmouth CF cohort vs. TEDDY and DIABIMMUNE cohorts

- Plots were generated from  `Lib/parse_metagenome.py`, `Lib/pcoa_analysis.py`, and `Lib/Phylo_PCoA.R`

### [Figure 2:](Figures/Figure_2.pdf) Age model results from supervised machine learning methods

A-B. The age model identifies a subset of species that are important for healthy development, and these species are reproducible across datasets. 

C. Some species fail to colonize the CF microbiome at the same times while other species successfully colonize the gut microbiome. 

D-E. Relative microbiota age is delayed in the Dartmouth CF cohort using both the TEDDY and DIABIMMUNE trained models 

- Plots generated from `Lib/InfantAgeModels.R`

### [Figure 3:](Figures/Figure_3.pdf) DMM clustering of TEDDY and DIABIMMUNE samples, compared to Dartmouth CF cohort

This figure shows the results of the unsupervised machine learning models created from the TEDDY and DIABIMMUNE cohorts, and how Dartmouth CF compares to these two models.

- Plot generated from `Lib/DMM_Clustering.R`

### [Figure 4:](Figures/Figure_4.pdf) Volcano plots of significant species 

This figure shows volcano plots of species abundance across phases (Developmental: 0-12 Months, Transitional: 13-18 Months, Stable: 19-36 Months) with the Dartmouth CF cohort separately compared with TEDDY and DIABIMMUNE.

- Plots were generated from `Lib/Humann_functions.R`

### [Figure 5:](Figures/Figure_5.pdf) Oral and Fungal abundance differences between datasets

A. First we show the oral microbe abundance across datasets and by phase.

- Oral microbe plot was generated from `Lib/Oral_Microbe_Functions.R`

B. Total fungal prevalence across datasets and phases

C. Fungal prevalence broken down by genera across datasets and phases

- Fungal plots were generated from `Lib/Fungal_Functions.R`

### [Figure 6:](Figures/Figure_6.pdf) *Faecalibacterium prausnitzii* strain diversity

A-B. *F. prauznitzii* strain prevalence across datasets

C. Phylogenetic tree showing *F. prausnitzii* strains detected in Dartmouth CF cohort clustering with publically available cat metagenomes

D. The number of Faecalibacterium clades present across datasets and by phase

- Clade figure was generated from `Lib/Faecalibacterium_functions.R`

## Supplemental Figures:

### [Supplemental Figure 1:](Figures/Figure_S1.pdf) Phylum level differences over time 

Supplementing Figure 1, we show phylum level changes over time in the Dartmouth CF, TEDDY and DIABIMMUNE cohorts. 

- Plots were generated from `Lib/Phylum_boxplots.R`

### [Supplemental Figure 2:](Figures/Figure_S2.pdf) Phylum level differences between CF cohorts

Here we compare the Dartmouth CF cohort to the Hayden CF cohort in the first 12 months of life. We also show PCoA results comparing our CF cohorts vs. TEDDY and DIABIMMUNE, and compare TEDDY and DIABIMMUNE to each other. 

A. Phylum level abundances of Dartmouth CF vs Hayden CF in the first 12 months

- Phylum plot was generated from `Lib/proteobacteria_plot_scripts.R`

- PCoA plots were generated from `Lib/parse_metagenome.py`, `Lib/pcoa_analysis.py`, and `Lib/Phylo_PCoA.R`

### [Supplemental Figure 3:](Figures/Figure_S3.pdf) Family level and diversity differences between cohorts

A. We further investigate our cohorts by time, and at the family level to begin isolating key differences within between our samples and healthy infant microbiomes.

- Family plot was generated from `Lib/Family_plot_functions.R`

B-E. We also explored alpha diversity between our Dartmouth CF, TEDDY, and DIABIMMUNE cohorts, and complemented these data with figures showing the number of species over time

- Plot was generated from `Lib/Alpha_Diversity_Functions.R`

### [Supplemental Figure 4:](Assembled_Figures/Figure_S4.pdf) Cross Dataset Validation

The age model is validated by comparing the models created with the TEDDY and DIABIMMUNE cohorts. 

- Plot generated from `Lib/InfantAgeModels.R`

### [Supplemental Figure 5:](Figures/Figure_S5,pdf) Heatmaps for DMM clusters

This figure shows the species abundances for each sample in each cluster resulting from the unsupervised machine learning methods used in this work.

- Plot generated from `Lib/DMM_Clustering.R`

### [Supplemental Figure 6:](Figures/Figure_S6.pdf) Cross dataset DMM cluster validation

Clustering DIABIMMUNE based on the TEDDY trained DMM model and vice versa in order to validate that both models recapitualte similar trends in non-CF cohorts. 

- Plot generated from `Lib/InfantAgeModels.R`

### [Supplemental Figure 7:](Figures/Figure_S7.pdf) DMM results of Leonard et al. Celiac Disease cohort

DMM clustering of Celiac Disease infants compared to the TEDDY and DIABIMMUNE cohorts

- Plot generated from `Lib/DMM_Clustering.R`

### [Supplemental Figure 8:](Figures/Figure_S8.pdf) HUMAnN metabolic analysis comparison of Dartmouth CF and TEDDY cohorts
Plots showing Dartmouth compared to TEDDY HUMAnN analysis in a heatmap, plus most important species based on a Wilcoxon ranked sum test. Then showing top species function distribution, followed by the top species based on average relative abundances for comparison. 

- Plot generated from `Lib/Humann_functions.R`

### [Supplemental Figure 9:](Figures/Figure_S9.pdf) HUMAnN metabolic analysis comparison of Dartmouth CF and DIABIMMUNE cohorts
Plots showing Dartmouth compared to DIABIMMUNE HUMAnN analysis in a heatmap, plus most important species based on a Wilcoxon ranked sum test. Then showing top species function distribution, followed by the top species based on average relative abundances for comparison. 

- Plot generated from `Lib/Humann_functions.R`

### [Supplemental Figure 10:](Figures/Figure_S10.pdf) Bar plots showing results of top most statistically significant species

A-B. Bar plots showing the distribution and abundances of KEGG pathways for the top 5 most statistically significant species across datasets and phases.

- Plots generated from `Lib/Humann_functions.R`

C-D. Boxplots showing relative abundance of top 5 most statistically significant species across datasets and phases 

- Plots generated from `Lib/Humann_functions.R`

### [Supplemental Figure 11:](Figures/Figure_S11.pdf) *Faecalibacterium prausnitzii* prevalence across datasets

A-B. *F. prausnitzii* prevalence across the TEDDY and DIABIMMUNE datasets, separated by developmental phase

- Plots generated from `Lib/metagenomic_strains.py`
