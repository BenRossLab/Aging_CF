## Code used to create Figures 1B-E and Supplemental Figures 2B-F
# PCoA data
library('ape')
library('phyloseq')
library('tidyverse')
library('vegan')
library('extrafont')
library('svglite')
library('ggplot2')
library('RColorBrewer')

# Running PCoA
tree <- read.tree(text=tree_str)
if (str_detect(tree$tip.label[1],fixed("|")))
  tree$tip.label <- sapply(tree$tip.label, function(x) str_split(x,fixed("|"))[[1]][7])

otu <- otu_table(as.matrix(df_sub), TRUE)
phylo <- phyloseq(otu, tree)
dists <- UniFrac(phylo)
pcoa.results <- pcoa(dists)

# Run significance testing
sig <- adonis(formula = t(dists) ~ Cohort, data = metadata)
pval <- sig$aov.tab$`Pr(>F)`[1]
print(pval)
# Setting up for initial plot
df_plot <- pcoa.results$vectors %>% data.frame %>% select(Axis.1, Axis.2)
df_plot$sample <- str_replace_all(rownames(df_plot), "^X","")
df_plot <- left_join(df_plot, metadata)

# calculate loading vector for plot
x <- pcoa.results
Y <- t(df_sub)
plot.axes <- c(1,2)

pr.coo <- x$vectors

n <- nrow(Y)
points.stand <- scale(pr.coo[, plot.axes])

S <- cov(Y, points.stand)
U <- S %*% diag((x$values$Eigenvalues[plot.axes]/(n - 
                                                    1))^(-0.5))

colnames(U) <- colnames(pr.coo[, plot.axes])
# Finding the top species, must calculate euclidean distance of the vectors
U_lengths <- data.frame(U)
# Finding euclidean distance and ordering by highest values 
U_lengths['Len'] <- sqrt(abs(U_lengths['Axis.1']*2 + U_lengths['Axis.2']*2))
lengths <- U_lengths[order(U_lengths['Len'],decreasing=TRUE),]
print(head(lengths))

## Plot by Cohort
p <- ggplot(df_plot, aes(x = Axis.1, y = Axis.2)) + geom_point(aes(fill=Cohort), pch=21, size=3, alpha = 0.7, color = "gray50") + 
  scale_fill_manual(values = c("#1F78B4", "#33A02C")) + theme_bw(16) + theme(axis.text=element_text(color="black"),
                                                                             panel.grid.major.x = element_line(linetype="dashed"),
                                                                             panel.grid.minor.x = element_blank(),
                                                                             panel.grid.minor.y = element_blank(),
                                                                             panel.grid.major.y = element_line(linetype="dashed")) + xlab("PCoA1") + ylab("PCoA2") +
  annotate(geom="text", x=min(df_plot$`Axis.1`), y=max(df_plot$`Axis.2`), label=paste0("P-value: ",pval), hjust=-0.05, vjust=1)

# Add top 3 vectors -- Change as necessary
p <- p + geom_segment(size=0.25, color="grey20", aes(x = 0, y = 0, xend = lengths["s__Bifidobacterium_breve",1], yend = lengths["s__Bifidobacterium_breve",2]),
                      arrow = arrow(length = unit(0.25, "cm"))) + annotate("text",x = lengths["s__Bifidobacterium_breve",1]+0.02, y = lengths["s__Bifidobacterium_breve",2]+0.01, label="B. breve") +
  geom_segment(size=0.25, color="grey20", aes(x = 0, y = 0, xend = lengths["s__Bifidobacterium_bifidum",1], yend = lengths["s__Bifidobacterium_bifidum",2]),
               arrow = arrow(length = unit(0.25, "cm"))) + annotate("text",x = lengths["s__Bifidobacterium_bifidum",1]+0.03, y = lengths["s__Bifidobacterium_bifidum",2]+0.01, label="B. bifidum") + 
  geom_segment(size=0.25, color="grey20", aes(x = 0, y = 0, xend = lengths["s__Escherichia_coli",1], yend = lengths["s__Escherichia_coli",2]),
               arrow = arrow(length = unit(0.25, "cm"))) + annotate("text",x = lengths["s__Escherichia_coli",1]+0.03, y = lengths["s__Escherichia_coli",2], label="E. coli")

## Plot by Age 
p <- ggplot(df_plot, aes(x = Axis.1, y = Axis.2)) + 
  geom_point(data = df_plot %>% filter(Cohort == "DIABIMMUNE"), aes(fill=Age), pch=21, size=3, alpha = 0.7, color = "gray50") + 
  scale_fill_fermenter(palette = "RdYlBu", breaks = c(3, 6, 9, 12, 18, 24, 30, 36), name = "Month") + 
  theme_bw(16) + theme(axis.text=element_text(color="black"),
                       panel.grid.major.x = element_line(linetype="dashed"),
                       panel.grid.minor.x = element_blank(),
                       panel.grid.minor.y = element_blank(),
                       panel.grid.major.y = element_line(linetype="dashed")) + xlab("PCoA1") + ylab("PCoA2") +
  annotate(geom="text", x=min(df_plot$`Axis.1`), y=max(df_plot$`Axis.2`), label=paste0("P-value: ",pval), hjust=-0.05, vjust=1)

p <- p + geom_point(data = df_plot %>% filter(Cohort == "Dartmouth"), fill = "gray50", pch=21, size=3, alpha = 0.7, color = "gray50")

# Add top species -- change as necessary 
p <- p + geom_segment(size=0.25, color="grey20", aes(x = 0, y = 0, xend = lengths["s__Bifidobacterium_breve",1], yend = lengths["s__Bifidobacterium_breve",2]),
                      arrow = arrow(length = unit(0.25, "cm"))) + annotate("text",x = lengths["s__Bifidobacterium_breve",1]+0.02, y = lengths["s__Bifidobacterium_breve",2]+0.01, label="B. breve") +
  geom_segment(size=0.25, color="grey20", aes(x = 0, y = 0, xend = lengths["s__Bifidobacterium_bifidum",1], yend = lengths["s__Bifidobacterium_bifidum",2]),
               arrow = arrow(length = unit(0.25, "cm"))) + annotate("text",x = lengths["s__Bifidobacterium_bifidum",1]+0.03, y = lengths["s__Bifidobacterium_bifidum",2]+0.01, label="B. bifidum") + 
  geom_segment(size=0.25, color="grey20", aes(x = 0, y = 0, xend = lengths["s__Escherichia_coli",1], yend = lengths["s__Escherichia_coli",2]),
               arrow = arrow(length = unit(0.25, "cm"))) + annotate("text",x = lengths["s__Escherichia_coli",1]+0.03, y = lengths["s__Escherichia_coli",2], label="E. coli")
